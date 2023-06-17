# Component-wise MLP of Tank, et. al. (2022)
# Relevant code adapted from https://github.com/iancovert/Neural-GC

import torch
import torch.nn as nn
import numpy as np
from copy import deepcopy

def activation_helper(activation, dim=None):
    if activation == 'sigmoid':
        act = nn.Sigmoid()
    elif activation == 'tanh':
        act = nn.Tanh()
    elif activation == 'relu':
        act = nn.ReLU()
    elif activation == 'leakyrelu':
        act = nn.LeakyReLU()
    elif activation is None:
        def act(x):
            return x
    else:
        raise ValueError('unsupported activation: %s' % activation)
    return act

class MLP(nn.Module):
    def __init__(self, num_series, lag, hidden, activation):
        super(MLP, self).__init__()
        self.activation = activation_helper(activation)

        # Set up network.
        layer = nn.Conv1d(num_series, hidden[0], lag)
        modules = [layer]

        for d_in, d_out in zip(hidden, hidden[1:] + [1]):
            layer = nn.Conv1d(d_in, d_out, 1)
            modules.append(layer)

        # Register parameters.
        self.layers = nn.ModuleList(modules)

    def forward(self, X):
        X = X.transpose(2, 1)
        for i, fc in enumerate(self.layers):
            if i != 0:
                X = self.activation(X)
            X = fc(X)

        return X.transpose(2, 1)

class cMLP(nn.Module):
    def __init__(self, num_series, lag, hidden, activation='relu'):
        '''
        cMLP model with one MLP per time series.
        Args:
          num_series: dimensionality of multivariate time series.
          lag: number of previous time points to use in prediction.
          hidden: list of number of hidden units per layer.
          activation: nonlinearity at each layer.
        '''
        super(cMLP, self).__init__()
        self.p = num_series
        self.lag = lag
        self.activation = activation_helper(activation)

        # Set up networks.
        self.networks = nn.ModuleList([
            MLP(num_series, lag, hidden, activation)
            for _ in range(num_series)])

    def forward(self, X):
        '''
        Perform forward pass.
        Args:
          X: torch tensor of shape (batch, T, p).
        '''
        return torch.cat([network(X) for network in self.networks], dim=2)

    def GC(self, threshold=True, ignore_lag=True):
        '''
        Extract learned Granger causality.
        Args:
          threshold: return norm of weights, or whether norm is nonzero.
          ignore_lag: if true, calculate norm of weights jointly for all lags.
        Returns:
          GC: (p x p) or (p x p x lag) matrix. In first case, entry (i, j)
            indicates whether variable j is Granger causal of variable i. In
            second case, entry (i, j, k) indicates whether it's Granger causal
            at lag k.
        '''
        if ignore_lag:
            GC = [torch.norm(net.layers[0].weight, dim=(0, 2))
                  for net in self.networks]
        else:
            GC = [torch.norm(net.layers[0].weight, dim=0)
                  for net in self.networks]
        GC = torch.stack(GC)
        if threshold:
            return (GC > 0).int()
        else:
            return GC
  
def prox_update(network, lam, lr, penalty):
    '''
    Perform in place proximal update on first layer weight matrix.
    Args:
      network: MLP network.
      lam: regularization parameter.
      lr: learning rate.
      penalty: one of GL (group lasso), GSGL (group sparse group lasso),
        H (hierarchical).
    '''
    W = network.layers[0].weight
    hidden, p, lag = W.shape
    if penalty == 'GL':
        norm = torch.norm(W, dim=(0, 2), keepdim=True)
        W.data = ((W / torch.clamp(norm, min=(lr * lam)))
                  * torch.clamp(norm - (lr * lam), min=0.0))
    elif penalty == 'GSGL':
        norm = torch.norm(W, dim=0, keepdim=True)
        W.data = ((W / torch.clamp(norm, min=(lr * lam)))
                  * torch.clamp(norm - (lr * lam), min=0.0))
        norm = torch.norm(W, dim=(0, 2), keepdim=True)
        W.data = ((W / torch.clamp(norm, min=(lr * lam)))
                  * torch.clamp(norm - (lr * lam), min=0.0))
    elif penalty == 'H':
        # Lowest indices along third axis touch most lagged values.
        for i in range(lag):
            norm = torch.norm(W[:, :, :(i + 1)], dim=(0, 2), keepdim=True)
            W.data[:, :, :(i+1)] = (
                (W.data[:, :, :(i+1)] / torch.clamp(norm, min=(lr * lam)))
                * torch.clamp(norm - (lr * lam), min=0.0))
    else:
        raise ValueError('unsupported penalty: %s' % penalty)
        
def regularize(network, lam, penalty):
     '''
     Calculate regularization term for first layer weight matrix.
     Args:
       network: MLP network.
       penalty: one of GL (group lasso), GSGL (group sparse group lasso),
         H (hierarchical).
     '''
     W = network.layers[0].weight
     hidden, p, lag = W.shape
     if penalty == 'GL':
         return lam * torch.sum(torch.norm(W, dim=(0, 2)))
     elif penalty == 'GSGL':
         return lam * (torch.sum(torch.norm(W, dim=(0, 2)))
                       + torch.sum(torch.norm(W, dim=0)))
     elif penalty == 'H':
         # Lowest indices along third axis touch most lagged values.
         return lam * sum([torch.sum(torch.norm(W[:, :, :(i+1)], dim=(0, 2)))
                           for i in range(lag)])
     else:
         raise ValueError('unsupported penalty: %s' % penalty)       
  
def ridge_regularize(network, lam):
    '''Apply ridge penalty at all subsequent layers.'''
    return lam * sum([torch.sum(fc.weight ** 2) for fc in network.layers[1:]])        
        
def restore_parameters(model, best_model):
    '''Move parameter values from best_model to model.'''
    for params, best_params in zip(model.parameters(), best_model.parameters()):
        params.data = best_params      
        
def train_model_ista(cmlp, X, lr, max_iter, lam=0, lam_ridge=0, penalty='H',
                     lookback=5, check_every=100, verbose=1):
    '''Train model with ISTA.'''
    lag = cmlp.lag
    p = X.shape[-1]
    loss_fn = nn.MSELoss(reduction='mean')
    train_loss_list = []

    # For early stopping.
    best_it = None
    best_loss = np.inf
    best_model = None

    # Calculate smooth error.
    loss = sum([loss_fn(cmlp.networks[i](X[:, :-1]), X[:, lag:, i:i+1])
                for i in range(p)])
    ridge = sum([ridge_regularize(net, lam_ridge) for net in cmlp.networks])
    smooth = loss + ridge

    for it in range(max_iter):
        # Take gradient step.
        smooth.backward()
        for param in cmlp.parameters():
            param.data = param - lr * param.grad

        # Take prox step.
        if lam > 0:
            for net in cmlp.networks:
                prox_update(net, lam, lr, penalty)

        cmlp.zero_grad()

        # Calculate loss for next iteration.
        loss = sum([loss_fn(cmlp.networks[i](X[:, :-1]), X[:, lag:, i:i+1])
                    for i in range(p)])
        ridge = sum([ridge_regularize(net, lam_ridge) for net in cmlp.networks])
        smooth = loss + ridge

        # Check progress.
        if (it + 1) % check_every == 0:
            # Add nonsmooth penalty.
            nonsmooth = sum([regularize(net, lam, penalty)
                             for net in cmlp.networks])
            mean_loss = (smooth + nonsmooth) / p
            train_loss_list.append(mean_loss.detach())

            if verbose > 0:
                print(('-' * 10 + 'Iter = %d' + '-' * 10) % (it + 1))
                print('Loss = %f' % mean_loss)
                print('Variable usage = %.2f%%'
                      % (100 * torch.mean(cmlp.GC().float())))

            # Check for early stopping.
            if mean_loss < best_loss:
                best_loss = mean_loss
                best_it = it
                best_model = deepcopy(cmlp)
            elif (it - best_it) == lookback * check_every:
                if verbose:
                    print('Stopping early')
                break

    # Restore best model.
    restore_parameters(cmlp, best_model)

    return train_loss_list


import scipy.io
mat = scipy.io.loadmat("dat.mat")['data_out']
dat_response = np.transpose(np.asmatrix(mat[((mat[:,0] >= 11.04) & (mat[:,0] < 14.7)), 4]))

# 30-40ms lag for stimuli is a 30ms initial selection with a 40 time step lag at 4kHz
dat_stimuli = np.asmatrix(mat[((mat[:,0] >= 11.01) & (mat[:,0] < 14.67)), 1:4])

data = np.concatenate((dat_response, dat_stimuli), axis = 1)

dat = torch.tensor(((data - np.mean(data, axis = 0)) / np.sqrt(np.var(data, axis = 0)))[np.newaxis], dtype=torch.float32)
dimension = data.shape[-1]

cmlp1 = cMLP(dimension, lag = 40, hidden = [250], activation = "tanh")
train_model_ista(cmlp1, dat, lam = 0.01, lam_ridge = 1e-2, lr = 1e-3, penalty = 'H', max_iter = 20000, check_every = 100)
grangerH1 = cmlp1.GC(threshold = False, ignore_lag = True).detach().numpy()
np.savetxt("grangerH1.csv", grangerH1, delimiter=",")

cmlp2 = cMLP(dimension, lag = 40, hidden = [250], activation = "tanh")
train_model_ista(cmlp2, dat, lam = 0.05, lam_ridge = 1e-2, lr = 1e-3, penalty = 'H', max_iter = 20000, check_every = 100)
grangerH2 = cmlp2.GC(threshold = False, ignore_lag = True).detach().numpy()
np.savetxt("grangerH2.csv", grangerH2, delimiter=",")

cmlp3 = cMLP(dimension, lag = 40, hidden = [250], activation = "tanh")
train_model_ista(cmlp3, dat, lam = 0.1, lam_ridge = 1e-2, lr = 1e-3, penalty = 'H', max_iter = 20000, check_every = 100)
grangerH3 = cmlp3.GC(threshold = False, ignore_lag = True).detach().numpy()
np.savetxt("grangerH3.csv", grangerH3, delimiter=",")

cmlp4 = cMLP(dimension, lag = 40, hidden = [250], activation = "tanh")
train_model_ista(cmlp4, dat, lam = 0.5, lam_ridge = 1e-2, lr = 1e-3, penalty = 'H', max_iter = 20000, check_every = 100)
grangerH4 = cmlp4.GC(threshold = False, ignore_lag = True).detach().numpy()
np.savetxt("grangerH4.csv", grangerH4, delimiter=",")













