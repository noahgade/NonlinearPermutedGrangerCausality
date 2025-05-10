#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Functions for Nonlinear Permuted Granger Causality

arma::mat Scale(arma::mat data) {
  // Standardizes data to mean zero and unit variance
  // Retains mean and scale vectors for inverse operation
  int dim = data.n_cols;
  arma::mat sdata = arma::mat(size(data)); 
  arma::vec means = arma::vec(dim);
  arma::vec scales = arma::vec(dim);
  for(int i = 0; i < dim; i++) {
    means(i) = arma::mean(data.col(i));
    scales(i) = arma::stddev(data.col(i));
    sdata.col(i) = (data.col(i) - means(i)) / scales(i);
  }
  return sdata;
}

double nrmse(arma::mat Mat, arma::mat RefMat) {
  // Computes a normalized root mean squared error to RefMat
  // Returns a value relative to the standard deviations of RefMat
  arma::vec DiffMat = arma::mean(arma::pow(Mat - RefMat, 2), 0);
  arma::vec NDiffMat = arma::vec(size(DiffMat));
  for(int c = 0; c < RefMat.n_cols; c++) {
    double scale = arma::stddev(RefMat.col(c));
    NDiffMat(c) = DiffMat(c) / scale;
  }
  double Out = mean(NDiffMat);
  return Out;
}

arma::mat Proximal(arma::mat Wi, double Threshold, int Penalty, arma::ivec VarKey, arma::ivec PenaltyKey) {
  // Computes proximal gradient for input matrix given a LASSO penalty
  // Penalty = 1 adds a LASSO penalty on all input weights Wi individually
  // Penalty = 2 adds a LASSO penalty on all input rows of Wi 
  // Penalty = 3 adds a group LASSO penalty by variable
  // Penalty = 4 adds a group sparse group LASSO penalty, combining 2 & 3
  // Penalty = 5 adds a hierarchical group LASSO penalty
  arma::mat Output = arma::mat(size(Wi));
  int NRows = Wi.n_rows;
  if(Penalty == 1) {
    for(int r = 0; r < NRows; r++) {
      for(int c = 0; c < Wi.n_cols; c++) {
        if(PenaltyKey(r) > 0) {
          double Value = 1 - Threshold / sqrt(pow(Wi(r, c), 2));
          if(Value < 0) {
            Value = 0;
          } 
          Output(r, c) = Value * Wi(r, c);
        } else {
          Output(r, c) = Wi(r, c);
        }
      }
    }
  } else if(Penalty == 2) {
    for(int r = 0; r < NRows; r++) {
      for(int c = 0; c < Wi.n_cols; c++) {
        if(PenaltyKey(r) > 0) {
          double Value = 1 - Threshold / sqrt(arma::accu(arma::pow(Wi.row(r), 2)));
          if(Value < 0) {
            Value = 0;
          } 
          Output(r, c) = Value * Wi(r, c);
        } else {
          Output(r, c) = Wi(r, c);
        }
      }
    }
  } else if(Penalty == 3) {
    for(int r = 0; r < NRows; r++) {
      for(int c = 0; c < Wi.n_cols; c++) {
        arma::uvec Group = find(VarKey == VarKey(r));
        if(PenaltyKey(r) > 0) {
          double Value = 1 - Threshold / sqrt(arma::accu(arma::pow(Wi.rows(Group), 2)));
          if(Value < 0) {
            Value = 0;
          } 
          Output(r, c) = Value * Wi(r, c);
        } else {
          Output(r, c) = Wi(r, c);
        }
      }
    }
  } else if(Penalty == 4) {
    for(int r = 0; r < Wi.n_rows; r++) {
      for(int c = 0; c < Wi.n_cols; c++) {
        arma::uvec Group = find(VarKey == VarKey(r));
        double enet = 0.5;
        if(PenaltyKey(r) > 0) {
          double Value = 1 - (1 - enet) * Threshold / sqrt(arma::accu(arma::pow(Wi.rows(Group), 2)));
          if(Value < 0) {
            Value = 0;
          } 
          double Value1 = 1 - enet * Threshold / sqrt(arma::accu(arma::pow(Wi.row(r), 2)));
          if(Value1 < 0) {
            Value1 = 0;
          }
          Output(r, c) = Value * Value1 * Wi(r, c);
        } else {
          Output(r, c) = Wi(r, c);
        }
      }
    }
  } else if(Penalty == 5) {
    for(int r = 0; r < Wi.n_rows; r++) {
      for(int c = 0; c < Wi.n_cols; c++) {
        arma::uvec Group = find(VarKey == VarKey(r));
        arma::uvec HGroup = Group.rows(find(Group >= r));
        if(PenaltyKey(r) > 0) {
          double Value = 1 - Threshold / sqrt(arma::accu(arma::pow(Wi.rows(HGroup), 2)));
          if(Value < 0) {
            Value = 0;
          }
          Output(r, c) = Value * Wi(r, c);
        } else {
          Output(r, c) = Wi(r, c);
        }
      }
    }
  } else {
    Output = Wi;
  }
  return Output;
}

List ForwardMLP1(arma::mat Input, arma::mat Wi, arma::mat Wo, int Activation) {
  // Forward propagation of a one-layer MLP
  // Activation Function
  // tanh = 0
  // ReLU = 1
  // Leaky ReLU = 2
  int L = Input.n_rows;
  int IDim = Input.n_cols;
  int HDim = Wi.n_cols;
  arma::mat H = arma::ones(L, 1) * Wi.row(0) + Input * Wi.rows(1, IDim);
  if(Activation == 1) {
    for(int r = 0; r < H.n_rows; r++) {
      for(int c = 0; c < H.n_cols; c++) {
        if(H(r, c) < 0) {
          H(r, c) = 0;
        }
      }
    }
  } else if(Activation == 2) {
    for(int r = 0; r < H.n_rows; r++) {
      for(int c = 0; c < H.n_cols; c++) {
        if(H(r, c) < 0) {
          H(r, c) = -0.01 * H(r, c);
        }
      }
    }
  } else {
    H = tanh(H);
  }
  arma::mat hOutput = arma::ones(L, 1) * Wo.row(0) + H * Wo.rows(1, HDim);
  return List::create(_["hOutput"] = hOutput, _["H"] = H);
}

List ForwardMLP2(arma::mat Input, arma::mat Wi, arma::mat Wh, arma::mat Wo, int Activation) {
  // Forward propagation of a two-layer MLP
  // Activation Function
  // tanh = 0
  // ReLU = 1
  // Leaky ReLU = 2
  int L = Input.n_rows;
  int IDim = Input.n_cols;
  int HDim = Wi.n_cols;
  arma::mat H1 = arma::ones(L, 1) * Wi.row(0) + Input * Wi.rows(1, IDim);
  if(Activation == 1) {
    for(int r = 0; r < H1.n_rows; r++) {
      for(int c = 0; c < H1.n_cols; c++) {
        if(H1(r, c) < 0) {
          H1(r, c) = 0;
        }
      }
    }
  } else if(Activation == 2) {
    for(int r = 0; r < H1.n_rows; r++) {
      for(int c = 0; c < H1.n_cols; c++) {
        if(H1(r, c) < 0) {
          H1(r, c) = -0.01 * H1(r, c);
        }
      }
    }
  } else {
    H1 = tanh(H1);
  }
  arma::mat H2 = arma::ones(L, 1) * Wh.row(0) + H1 * Wh.rows(1, HDim);
  if(Activation == 1) {
    for(int r = 0; r < H2.n_rows; r++) {
      for(int c = 0; c < H2.n_cols; c++) {
        if(H2(r, c) < 0) {
          H2(r, c) = 0;
        }
      }
    }
  } else if(Activation == 2) {
    for(int r = 0; r < H2.n_rows; r++) {
      for(int c = 0; c < H2.n_cols; c++) {
        if(H2(r, c) < 0) {
          H2(r, c) = -0.01 * H2(r, c);
        }
      }
    }
  } else {
    H2 = tanh(H2);
  }
  arma::mat hOutput = arma::ones(L, 1) * Wo.row(0) + H2 * Wo.rows(1, HDim);
  return List::create(_["hOutput"] = hOutput, _["H2"] = H2, _["H1"] = H1);
}

List TrainMLP(arma::mat trainInput, arma::mat trainOutput, arma::mat validInput, arma::mat validOutput, int NLayers, int HDim, int Penalty, double Lambda, arma::ivec VarKey, arma::ivec PenaltyKey, arma::mat Wi, arma::mat Wh, arma::mat Wo, int Activation) {
  // Trained MLP with depth NLayers and width HDim
  // Training with batch Adam proximal gradient descent
  // Penalty > 0 adds a LASSO penalty of Lambda to the training
  // Penalty = 1 adds a LASSO penalty on all input weights Wi individually
  // Penalty = 2 adds a LASSO penalty on all input rows of Wi 
  // Penalty = 3 adds a group LASSO penalty by variable
  // Penalty = 4 adds a group sparse group LASSO penalty, combining 2 & 3
  // Penalty = 5 adds a hierarchical group LASSO penalty
  int max_iter = 20000; // Maximum number of training iterations
  int patience = 500; // Continuation of training past minimum
  double Epsilon = 0.001; // Learning Rate
  double Rho1 = 0.9; // Decay, moment 1
  double Rho2 = 0.999; // Decay, moment 2
  double Delta = 1e-8; // Stabilization term
  double Ctol = 1e-4; // Tolerance for Condition goal
  double Ntol = 1e-2; // Tolerance for NRMSE goal
  double Condition = 999;
  double NRMSE = 999;
  int NTrain = trainInput.n_rows;
  int IDim = trainInput.n_cols;
  int ODim = trainOutput.n_cols;
  arma::mat LossMat = Condition * arma::ones(max_iter + 1, 3);
  arma::mat tWi = Wi;
  arma::mat tWh = Wh;
  arma::mat tWo = Wo;
  arma::mat SWi = arma::zeros(IDim + 1, HDim);
  arma::mat RWi = arma::zeros(IDim + 1, HDim); 
  arma::mat dWi = arma::zeros(IDim + 1, HDim); 
  arma::mat SWh = arma::zeros(HDim + 1, HDim);
  arma::mat RWh = arma::zeros(HDim + 1, HDim); 
  arma::mat dWh = arma::zeros(HDim + 1, HDim); 
  arma::mat SWo = arma::zeros(HDim + 1, ODim);
  arma::mat RWo = arma::zeros(HDim + 1, ODim); 
  arma::mat dWo = arma::zeros(HDim + 1, ODim);
  int iteration = 0;
  if(NLayers < 2) {
    List TrainPass = ForwardMLP1(trainInput, tWi, tWo, Activation);
    arma::mat hSOTrain = TrainPass[0];
    arma::mat HTrain = TrainPass[1];
    List ValidPass = ForwardMLP1(validInput, tWi, tWo, Activation);
    arma::mat hSOValid = ValidPass[0];
    arma::mat TrainLoss = arma::pow((trainOutput - hSOTrain), 2);
    arma::mat ValidLoss = arma::pow((validOutput - hSOValid), 2);
    LossMat(iteration, 0) = arma::accu(TrainLoss) / TrainLoss.n_elem;
    LossMat(iteration, 1) = arma::accu(ValidLoss) / ValidLoss.n_elem;
    LossMat(iteration, 2) = nrmse(hSOValid, validOutput);
    double MinLoss = LossMat(iteration, 1);
    while((Condition > Ctol) and (NRMSE > Ntol)) {
      arma::cube dWi_cube = arma::zeros(IDim + 1, HDim, NTrain);
      arma::cube dWo_cube = arma::zeros(HDim + 1, ODim, NTrain);
      arma::mat dLdhO = -2 * (trainOutput - hSOTrain);
      for(int time = 0; time < NTrain; time++) {
        arma::mat dLdH = dLdhO.row(time) * trans(tWo.rows(1, HDim));
        arma::mat dLdWo = trans(arma::join_rows(arma::ones(1, 1), HTrain.row(time))) * dLdhO.row(time);
        arma::mat dHTrain = arma::mat(size(HTrain.row(time)));
        if(Activation == 1) {
          for(int e = 0; e < HTrain.n_cols; e++) {
            if(HTrain(time, e) < 0) {
              dHTrain(e) = 0;
            } else {
              dHTrain(e) = 1;
            }
          }
        } else if(Activation == 2) {
          for(int e = 0; e < HTrain.n_cols; e++) {
            if(HTrain(time, e) < 0) {
              dHTrain(e) = -0.01;
            } else {
              dHTrain(e) = 1;
            }
          }
        } else {
          dHTrain = 1 - pow(HTrain.row(time), 2);
        }
        arma::mat dLdWi = trans(arma::join_rows(arma::ones(1, 1), trainInput.row(time))) * (dHTrain % dLdH); 
        dWi_cube.slice(time) = dLdWi;
        dWo_cube.slice(time) = dLdWo;
      }
      arma::mat gWi = arma::mean(dWi_cube, 2); 
      arma::mat gWo = arma::mean(dWo_cube, 2);
      iteration = iteration + 1;
      SWi = Rho1 * SWi + (1 - Rho1) * gWi;
      SWo = Rho1 * SWo + (1 - Rho1) * gWo;
      RWi = Rho2 * RWi + (1 - Rho2) * (gWi % gWi);
      RWo = Rho2 * RWo + (1 - Rho2) * (gWo % gWo);
      arma::mat ratioWi = ((1 - pow(Rho1, iteration)) * SWi) / ((1 - pow(Rho2, iteration)) * RWi + Delta);
      arma::mat ratioWo = ((1 - pow(Rho1, iteration)) * SWo) / ((1 - pow(Rho2, iteration)) * RWo + Delta);
      dWi = Epsilon * ratioWi.clamp(-1, 1);
      dWo = Epsilon * ratioWo.clamp(-1, 1);
      tWi = tWi - dWi;
      double Threshold = Lambda * Epsilon;
      tWi = Proximal(tWi, Threshold, Penalty, VarKey, PenaltyKey);
      tWo = tWo - dWo;
      TrainPass = ForwardMLP1(trainInput, tWi, tWo, Activation);
      arma::mat hSOT = TrainPass[0];
      hSOTrain = hSOT;
      arma::mat HT =  TrainPass[1];
      HTrain = HT;
      ValidPass = ForwardMLP1(validInput, tWi, tWo, Activation);
      arma::mat hSOV = ValidPass[0];
      hSOValid = hSOV;
      TrainLoss = arma::pow((trainOutput - hSOTrain), 2);
      ValidLoss = arma::pow((validOutput - hSOValid), 2);
      LossMat(iteration, 0) = arma::accu(TrainLoss) / TrainLoss.n_elem;
      LossMat(iteration, 1) = arma::accu(ValidLoss) / ValidLoss.n_elem;
      LossMat(iteration, 2) = nrmse(hSOValid, validOutput);
      if(LossMat(iteration, 1) < MinLoss) {
        MinLoss = LossMat(iteration, 1);
        Wi = tWi;
        Wo = tWo;
      }
      if(iteration > max_iter) {
        break;
      } else if(iteration >= patience) {
        Condition = LossMat(iteration - patience, 1) - LossMat(iteration, 1);
        NRMSE = LossMat(iteration, 2);
      }
    }
  } else {
    List TrainPass = ForwardMLP2(trainInput, tWi, tWh, tWo, Activation);
    arma::mat hSOTrain = TrainPass[0];
    arma::mat H2Train = TrainPass[1];
    arma::mat H1Train = TrainPass[2];
    List ValidPass = ForwardMLP2(validInput, tWi, tWh, tWo, Activation);
    arma::mat hSOValid = ValidPass[0];
    arma::mat TrainLoss = arma::pow((trainOutput - hSOTrain), 2);
    arma::mat ValidLoss = arma::pow((validOutput - hSOValid), 2);
    LossMat(iteration, 0) = arma::accu(TrainLoss) / TrainLoss.n_elem;
    LossMat(iteration, 1) = arma::accu(ValidLoss) / ValidLoss.n_elem;
    LossMat(iteration, 2) = nrmse(hSOValid, validOutput);
    double MinLoss = LossMat(iteration, 1);
    while((Condition > Ctol) and (NRMSE > Ntol)) {
      arma::cube dWi_cube = arma::zeros(IDim + 1, HDim, NTrain);
      arma::cube dWh_cube = arma::zeros(HDim + 1, HDim, NTrain);
      arma::cube dWo_cube = arma::zeros(HDim + 1, ODim, NTrain);
      arma::mat dLdhO = -2 * (trainOutput - hSOTrain);
      for(int time = 0; time < NTrain; time++) {
        arma::mat dLdH2 = dLdhO.row(time) * trans(tWo.rows(1, HDim));
        arma::mat dLdWo = trans(arma::join_rows(arma::ones(1, 1), H2Train.row(time))) * dLdhO.row(time);
        arma::mat dH2Train = arma::mat(size(H2Train.row(time)));
        arma::mat dH1Train = arma::mat(size(H1Train.row(time)));
        if(Activation == 1) {
          for(int e = 0; e < H2Train.n_cols; e++) {
            if(H2Train(time, e) < 0) {
              dH2Train(e) = 0;
            } else {
              dH2Train(e) = 1;
            }
          }
          for(int e = 0; e < H1Train.n_cols; e++) {
            if(H1Train(time, e) < 0) {
              dH1Train(e) = 0;
            } else {
              dH1Train(e) = 1;
            }
          }
        } else if(Activation == 2) {
          for(int e = 0; e < H2Train.n_cols; e++) {
            if(H2Train(time, e) < 0) {
              dH2Train(e) = -0.01;
            } else {
              dH2Train(e) = 1;
            }
          }
          for(int e = 0; e < H1Train.n_cols; e++) {
            if(H1Train(time, e) < 0) {
              dH1Train(e) = -0.01;
            } else {
              dH1Train(e) = 1;
            }
          }
        } else {
          dH2Train = 1 - pow(H2Train.row(time), 2);
          dH1Train = 1 - pow(H1Train.row(time), 2);
        }
        arma::mat dLdWh = trans(arma::join_rows(arma::ones(1, 1), H1Train.row(time))) * (dH2Train % dLdH2);
        arma::mat dLdH1 = (dH2Train % dLdH2) * trans(tWh.rows(1, HDim));
        arma::mat dLdWi = trans(arma::join_rows(arma::ones(1, 1), trainInput.row(time))) * (dH1Train % dLdH1);
        dWi_cube.slice(time) = dLdWi;
        dWh_cube.slice(time) = dLdWh;
        dWo_cube.slice(time) = dLdWo;
      }
      arma::mat gWi = arma::mean(dWi_cube, 2); 
      arma::mat gWh = arma::mean(dWh_cube, 2); 
      arma::mat gWo = arma::mean(dWo_cube, 2);
      iteration = iteration + 1;
      SWi = Rho1 * SWi + (1 - Rho1) * gWi;
      SWh = Rho1 * SWh + (1 - Rho1) * gWh;
      SWo = Rho1 * SWo + (1 - Rho1) * gWo;
      RWi = Rho2 * RWi + (1 - Rho2) * (gWi % gWi);
      RWh = Rho2 * RWh + (1 - Rho2) * (gWh % gWh);
      RWo = Rho2 * RWo + (1 - Rho2) * (gWo % gWo);
      arma::mat ratioWi = ((1 - pow(Rho1, iteration)) * SWi) / ((1 - pow(Rho2, iteration)) * RWi + Delta);
      arma::mat ratioWh = ((1 - pow(Rho1, iteration)) * SWh) / ((1 - pow(Rho2, iteration)) * RWh + Delta);
      arma::mat ratioWo = ((1 - pow(Rho1, iteration)) * SWo) / ((1 - pow(Rho2, iteration)) * RWo + Delta);
      dWi = Epsilon * ratioWi.clamp(-1, 1);
      dWh = Epsilon * ratioWh.clamp(-1, 1);
      dWo = Epsilon * ratioWo.clamp(-1, 1);
      tWi = tWi - dWi;
      double Threshold = Lambda * Epsilon;
      tWi = Proximal(tWi, Threshold, Penalty, VarKey, PenaltyKey);
      tWh = tWh - dWh;
      tWo = tWo - dWo;
      TrainPass = ForwardMLP2(trainInput, tWi, tWh, tWo, Activation);
      arma::mat hSOT = TrainPass[0];
      hSOTrain = hSOT;
      arma::mat H2T =  TrainPass[1];
      H2Train = H2T;
      arma::mat H1T =  TrainPass[2];
      H1Train = H1T;
      ValidPass = ForwardMLP2(validInput, tWi, tWh, tWo, Activation);
      arma::mat hSOV = ValidPass[0];
      hSOValid = hSOV;
      TrainLoss = arma::pow((trainOutput - hSOTrain), 2);
      ValidLoss = arma::pow((validOutput - hSOValid), 2);
      LossMat(iteration, 0) = arma::accu(TrainLoss) / TrainLoss.n_elem;
      LossMat(iteration, 1) = arma::accu(ValidLoss) / ValidLoss.n_elem;
      LossMat(iteration, 2) = nrmse(hSOValid, validOutput);
      if(LossMat(iteration, 1) < MinLoss) {
        MinLoss = LossMat(iteration, 1);
        Wi = tWi;
        Wh = tWh;
        Wo = tWo;
      }
      if(iteration > max_iter) {
        break;
      } else if(iteration >= patience) {
        Condition = LossMat(iteration - patience, 1) - LossMat(iteration, 1);
        NRMSE = LossMat(iteration, 2);
      }
    }
  }
  arma::mat finalLossMat = LossMat.rows(find(LossMat.col(0) < 999));
  return List::create(_["Wi"] = Wi, _["LossMat"] = finalLossMat, _["Wh"] = Wh, _["Wo"] = Wo, _["L"] = NLayers, _["N"] = HDim);
}

List TrainELM(arma::mat Input, arma::mat Output, int HDim, double Regular, arma::mat Wi, int Activation) {
  // Trained MLP with depth NLayers and width HDim
  // Using Ridge regression with regularization parameter Regular
  int L = Input.n_rows;
  int ODim = Output.n_cols;
  arma::mat Wo = arma::zeros(HDim + 1, ODim);
  List Forward = ForwardMLP1(Input, Wi, Wo, Activation);
  arma::mat H = Forward[1];
  arma::mat AugH =  arma::join_rows(arma::ones(L, 1), H);
  Wo = arma::inv(trans(AugH) * AugH + Regular * arma::eye(HDim + 1, HDim + 1)) * trans(AugH) * Output;
  arma::mat hSOData = arma::ones(L, 1) * Wo.row(0) + H * Wo.rows(1, HDim);
  arma::mat LossMat = arma::pow((Output - hSOData), 2); 
  double Loss = arma::accu(LossMat) / LossMat.n_elem;
  return List::create(_["Wi"] = Wi, _["Loss"] = Loss, _["Wo"] = Wo, _["N"] = HDim);
}

// [[Rcpp::export]]
arma::mat NPGC(arma::mat Y, arma::mat Ylag, arma::mat Z, arma::mat X, int Type, int Omega, int M, int K, int R, int HDim, int Activation) {
  // Nonlinear Permuted Granger Causality
  // Type = 0, ELM architecture
  // Type = 1, FNN architecture
  arma::mat SY = Scale(Y);
  arma::mat SYlag = Scale(Ylag);
  arma::mat SZ = Scale(Z);
  arma::mat SX = Scale(X);
  arma::mat Final = arma::zeros(M + 1, R);
  int L = SY.n_rows / Omega;
  int Penalty = 0;
  double Lambda = 0;
  arma::ivec VarKey = arma::join_cols(arma::zeros<arma::ivec>(1), arma::join_cols(1 * arma::ones<arma::ivec>(SYlag.n_cols), arma::join_cols(2 * arma::ones<arma::ivec>(SZ.n_cols), 3 * arma::ones<arma::ivec>(SX.n_cols))));
  arma::ivec PenaltyKey = arma::zeros<arma::ivec>(1 + SYlag.n_cols + SZ.n_cols + SX.n_cols);
  arma::uvec FoldKey = arma::shuffle(arma::floor(K * arma::regspace<arma::uvec>(0, L - 1) / L));
  arma::imat Folds = arma::zeros<arma::imat>(L * Omega, K);
  for(int k = 0; k < K; k++) {
    arma::uvec Key = arma::shuffle(arma::find(FoldKey != k));
    int NKey = Key.n_rows;
    arma::uvec TrainKey = Key.rows(0, ceil(4 * NKey / 5) - 1);
    arma::uvec ValidKey = Key.rows(ceil(4 * NKey / 5), NKey - 1);
    for(int rT : TrainKey) {
      Folds(rT, k) = 2;
    }
    for(int rV : ValidKey) {
      Folds(rV, k) = 1;
    }
  }
  arma::umat Permutations = arma::umat(L * Omega, M + 1);
  Permutations.rows(0, L - 1).col(0) = arma::regspace<arma::uvec>(0, L - 1);
  for(int m = 1; m <= M; m++) {
    Permutations.rows(0, L - 1).col(m) = arma::randperm(L);
  }
  if(Omega > 1) {
    for(int w = 1; w < Omega; w++) {
      Folds.rows(w * L, (w + 1) * L - 1) = Folds.rows(0, L - 1);
      Permutations.rows(w * L, (w + 1) * L - 1) = Permutations.rows(0, L - 1) + w * L;
    }
  }
  arma::cube Input = arma::join_slices(arma::join_rows(SYlag, arma::join_rows(SZ, SX)), arma::zeros(SX.n_rows, SYlag.n_cols + SZ.n_cols + SX.n_cols, M));
  for(int m = 1; m <= M; m++) {
    Input.slice(m) = arma::join_rows(SYlag, arma::join_rows(SZ, SX.rows(Permutations.col(m))));
  }
  int IDim = Input.n_cols;
  int ODim = SY.n_cols;
  if(Type == 0) {
    arma::cube Wi = arma::randn(IDim + 1, HDim, R);
    for(int r = 0; r < R; r++) {
      for(int m = 0; m <= M; m++) {
        for(int k = 0; k < K; k++) {
          arma::mat trainInput = Input.slice(m).rows(arma::find(Folds.col(k) != 0));
          arma::mat trainOutput = SY.rows(arma::find(Folds.col(k) != 0));
          arma::mat testInput = Input.slice(m).rows(arma::find(Folds.col(k) == 0));
          arma::mat testOutput = SY.rows(arma::find(Folds.col(k) == 0));
          List ELM = TrainELM(trainInput, trainOutput, HDim, 1e-8, Wi.slice(r), Activation);
          arma::mat fWi = ELM[0];
          arma::mat fWo = ELM[2];
          List Forward = ForwardMLP1(testInput, fWi, fWo, Activation);
          arma::mat Hat = Forward[0];
          double Residuals = arma::accu(arma::pow((testOutput - Hat), 2)) / testOutput.n_rows;
          Final(m, r) = Final(m, r) + Residuals / K;
        }
      }
    }
  } else {
    double Xavier = sqrt(2 / (float(IDim) + float(ODim)));
    arma::cube Wi = Xavier * arma::randn(IDim + 1, HDim, R);
    arma::cube Wh = Xavier * arma::randn(HDim + 1, HDim, R);
    arma::cube Wo = Xavier * arma::randn(HDim + 1, ODim, R);
    for(int r = 0; r < R; r++) {
      for(int m = 0; m <= M; m++) {
        for(int k = 0; k < K; k++) {
          arma::mat trainInput = Input.slice(m).rows(arma::find(Folds.col(k) == 2));
          arma::mat trainOutput = SY.rows(arma::find(Folds.col(k) == 2));
          arma::mat validInput = Input.slice(m).rows(arma::find(Folds.col(k) == 1));
          arma::mat validOutput = SY.rows(arma::find(Folds.col(k) == 1));
          arma::mat testInput = Input.slice(m).rows(arma::find(Folds.col(k) == 0));
          arma::mat testOutput = SY.rows(arma::find(Folds.col(k) == 0));
          List MLP = TrainMLP(trainInput, trainOutput, validInput, validOutput, Type, HDim, Penalty, Lambda, VarKey, PenaltyKey, Wi.slice(r), Wh.slice(r), Wo.slice(r), Activation);
          arma::mat fWi = MLP[0];
          arma::mat fWh = MLP[2];
          arma::mat fWo = MLP[3];
          List Forward = ForwardMLP1(testInput, fWi, fWo, Activation);
          arma::mat Hat = Forward[0];
          double Residuals = arma::accu(arma::pow((testOutput - Hat), 2)) / testOutput.n_rows;
          Final(m, r) = Final(m, r) + Residuals / K;
        }
      }
    }
  }
 return Final;
}

// [[Rcpp::export]]
arma::mat GAUSS(arma::mat Y, arma::mat Ylag, arma::mat Z, arma::mat X, int Type, int Omega, int M, int K, int R, int HDim, int Activation) {
  // Gaussian Substitution Granger Causality
  // Type = 0, ELM architecture
  // Type = 1, FNN architecture
  arma::mat SY = Scale(Y);
  arma::mat SYlag = Scale(Ylag);
  arma::mat SZ = Scale(Z);
  arma::mat SX = Scale(X);
  arma::mat Final = arma::zeros(M + 1, R);
  int L = SY.n_rows / Omega;
  int Penalty = 0;
  double Lambda = 0;
  arma::ivec VarKey = arma::join_cols(arma::zeros<arma::ivec>(1), arma::join_cols(1 * arma::ones<arma::ivec>(SYlag.n_cols), arma::join_cols(2 * arma::ones<arma::ivec>(SZ.n_cols), 3 * arma::ones<arma::ivec>(SX.n_cols))));
  arma::ivec PenaltyKey = arma::zeros<arma::ivec>(1 + SYlag.n_cols + SZ.n_cols + SX.n_cols);
  arma::uvec FoldKey = arma::shuffle(arma::floor(K * arma::regspace<arma::uvec>(0, L - 1) / L));
  arma::imat Folds = arma::zeros<arma::imat>(L * Omega, K);
  for(int k = 0; k < K; k++) {
    arma::uvec Key = arma::shuffle(arma::find(FoldKey != k));
    int NKey = Key.n_rows;
    arma::uvec TrainKey = Key.rows(0, ceil(4 * NKey / 5) - 1);
    arma::uvec ValidKey = Key.rows(ceil(4 * NKey / 5), NKey - 1);
    for(int rT : TrainKey) {
      Folds(rT, k) = 2;
    }
    for(int rV : ValidKey) {
      Folds(rV, k) = 1;
    }
  }
  if(Omega > 1) {
    for(int w = 1; w < Omega; w++) {
      Folds.rows(w * L, (w + 1) * L - 1) = Folds.rows(0, L - 1);
    }
  }
  arma::cube Input = arma::join_slices(arma::join_rows(SYlag, arma::join_rows(SZ, SX)), arma::cube(SX.n_rows, SYlag.n_cols + SZ.n_cols + SX.n_cols, M));
  for(int m = 1; m <= M; m++) {
    Input.slice(m) = arma::join_rows(SYlag, arma::join_rows(SZ, arma::randn(size(SX))));
  }
  int IDim = Input.n_cols;
  int ODim = SY.n_cols;
  if(Type == 0) {
    arma::cube Wi = arma::randn(IDim + 1, HDim, R);
    for(int r = 0; r < R; r++) {
      for(int m = 0; m <= M; m++) {
        for(int k = 0; k < K; k++) {
          arma::mat trainInput = Input.slice(m).rows(arma::find(Folds.col(k) != 0));
          arma::mat trainOutput = SY.rows(arma::find(Folds.col(k) != 0));
          arma::mat testInput = Input.slice(m).rows(arma::find(Folds.col(k) == 0));
          arma::mat testOutput = SY.rows(arma::find(Folds.col(k) == 0));
          List ELM = TrainELM(trainInput, trainOutput, HDim, 1e-8, Wi.slice(r), Activation);
          arma::mat fWi = ELM[0];
          arma::mat fWo = ELM[2];
          List Forward = ForwardMLP1(testInput, fWi, fWo, 0);
          arma::mat Hat = Forward[0];
          double Residuals = arma::accu(arma::pow((testOutput - Hat), 2)) / testOutput.n_rows;
          Final(m, r) = Final(m, r) + Residuals / K;
        }
      }
    }
  } else {
    double Xavier = sqrt(2 / (float(IDim) + float(ODim)));
    arma::cube Wi = Xavier * arma::randn(IDim + 1, HDim, R);
    arma::cube Wh = Xavier * arma::randn(HDim + 1, HDim, R);
    arma::cube Wo = Xavier * arma::randn(HDim + 1, ODim, R);
    for(int r = 0; r < R; r++) {
      for(int m = 0; m <= M; m++) {
        for(int k = 0; k < K; k++) {
          arma::mat trainInput = Input.slice(m).rows(arma::find(Folds.col(k) == 2));
          arma::mat trainOutput = SY.rows(arma::find(Folds.col(k) == 2));
          arma::mat validInput = Input.slice(m).rows(arma::find(Folds.col(k) == 1));
          arma::mat validOutput = SY.rows(arma::find(Folds.col(k) == 1));
          arma::mat testInput = Input.slice(m).rows(arma::find(Folds.col(k) == 0));
          arma::mat testOutput = SY.rows(arma::find(Folds.col(k) == 0));
          List MLP = TrainMLP(trainInput, trainOutput, validInput, validOutput, Type, HDim, Penalty, Lambda, VarKey, PenaltyKey, Wi.slice(r), Wh.slice(r), Wo.slice(r), Activation);
          arma::mat fWi = MLP[0];
          arma::mat fWh = MLP[2];
          arma::mat fWo = MLP[3];
          List Forward = ForwardMLP1(testInput, fWi, fWo, Activation);
          arma::mat Hat = Forward[0];
          double Residuals = arma::accu(arma::pow((testOutput - Hat), 2)) / testOutput.n_rows;
          Final(m, r) = Final(m, r) + Residuals / K;
        }
      }
    }
  }
  return Final;
}

// [[Rcpp::export]]
arma::mat ZERO(arma::mat Y, arma::mat Ylag, arma::mat Z, arma::mat X, int Type, int Omega, int K, int R, int HDim, int Activation) {
  // Zero Substitution Granger Causality
  // Type = 0, ELM architecture
  // Type = 1, FNN architecture
  arma::mat SY = Scale(Y);
  arma::mat SYlag = Scale(Ylag);
  arma::mat SZ = Scale(Z);
  arma::mat SX = Scale(X);
  arma::mat Final = arma::zeros(2, R);
  int L = SY.n_rows / Omega;
  int Penalty = 0;
  double Lambda = 0;
  arma::ivec VarKey = arma::join_cols(arma::zeros<arma::ivec>(1), arma::join_cols(1 * arma::ones<arma::ivec>(SYlag.n_cols), arma::join_cols(2 * arma::ones<arma::ivec>(SZ.n_cols), 3 * arma::ones<arma::ivec>(SX.n_cols))));
  arma::ivec PenaltyKey = arma::zeros<arma::ivec>(1 + SYlag.n_cols + SZ.n_cols + SX.n_cols);
  arma::uvec FoldKey = arma::shuffle(arma::floor(K * arma::regspace<arma::uvec>(0, L - 1) / L));
  arma::imat Folds = arma::zeros<arma::imat>(L * Omega, K);
  for(int k = 0; k < K; k++) {
    arma::uvec Key = arma::shuffle(arma::find(FoldKey != k));
    int NKey = Key.n_rows;
    arma::uvec TrainKey = Key.rows(0, ceil(4 * NKey / 5) - 1);
    arma::uvec ValidKey = Key.rows(ceil(4 * NKey / 5), NKey - 1);
    for(int rT : TrainKey) {
      Folds(rT, k) = 2;
    }
    for(int rV : ValidKey) {
      Folds(rV, k) = 1;
    }
  }
  if(Omega > 1) {
    for(int w = 1; w < Omega; w++) {
      Folds.rows(w * L, (w + 1) * L - 1) = Folds.rows(0, L - 1);
    }
  }
  arma::cube Input = arma::join_slices(arma::join_rows(SYlag, arma::join_rows(SZ, SX)), arma::join_rows(SYlag, arma::join_rows(SZ, arma::zeros(size(SX)))));
  int IDim = Input.n_cols;
  int ODim = SY.n_cols;
  if(Type == 0) {
    arma::cube Wi = arma::randn(IDim + 1, HDim, R);
    for(int r = 0; r < R; r++) {
      for(int m = 0; m < 2; m++) {
        for(int k = 0; k < K; k++) {
          arma::mat trainInput = Input.slice(m).rows(arma::find(Folds.col(k) != 0));
          arma::mat trainOutput = SY.rows(arma::find(Folds.col(k) != 0));
          arma::mat testInput = Input.slice(m).rows(arma::find(Folds.col(k) == 0));
          arma::mat testOutput = SY.rows(arma::find(Folds.col(k) == 0));
          List ELM = TrainELM(trainInput, trainOutput, HDim, 1e-8, Wi.slice(r), Activation);
          arma::mat fWi = ELM[0];
          arma::mat fWo = ELM[2];
          List Forward = ForwardMLP1(testInput, fWi, fWo, Activation);
          arma::mat Hat = Forward[0];
          double Residuals = arma::accu(arma::pow((testOutput - Hat), 2)) / testOutput.n_rows;
          Final(m, r) = Final(m, r) + Residuals / K;
        }
      }
    }
  } else {
    double Xavier = sqrt(2 / (float(IDim) + float(ODim)));
    arma::cube Wi = Xavier * arma::randn(IDim + 1, HDim, R);
    arma::cube Wh = Xavier * arma::randn(HDim + 1, HDim, R);
    arma::cube Wo = Xavier * arma::randn(HDim + 1, ODim, R);
    for(int r = 0; r < R; r++) {
      for(int m = 0; m < 2; m++) {
        for(int k = 0; k < K; k++) {
          arma::mat trainInput = Input.slice(m).rows(arma::find(Folds.col(k) == 2));
          arma::mat trainOutput = SY.rows(arma::find(Folds.col(k) == 2));
          arma::mat validInput = Input.slice(m).rows(arma::find(Folds.col(k) == 1));
          arma::mat validOutput = SY.rows(arma::find(Folds.col(k) == 1));
          arma::mat testInput = Input.slice(m).rows(arma::find(Folds.col(k) == 0));
          arma::mat testOutput = SY.rows(arma::find(Folds.col(k) == 0));
          List MLP = TrainMLP(trainInput, trainOutput, validInput, validOutput, Type, HDim, Penalty, Lambda, VarKey, PenaltyKey, Wi.slice(r), Wh.slice(r), Wo.slice(r), Activation);
          arma::mat fWi = MLP[0];
          arma::mat fWh = MLP[2];
          arma::mat fWo = MLP[3];
          List Forward = ForwardMLP1(testInput, fWi, fWo, Activation);
          arma::mat Hat = Forward[0];
          double Residuals = arma::accu(arma::pow((testOutput - Hat), 2)) / testOutput.n_rows;
          Final(m, r) = Final(m, r) + Residuals / K;
        }
      }
    }
  }
  return Final;
}

// [[Rcpp::export]]
List LASSO(arma::mat Y, arma::mat Ylag, arma::mat Z, arma::mat X, int Type, int Penalty, arma::vec Lambdas, int Omega, int K, int R, int HDim, bool ChooseLambda, int Activation) {
  // LASSO Granger Causality
  // Type = 1, FNN architecture
  // Type = 2, MLP architecture
  if(Type < 2) {
    if(ChooseLambda) {
      arma::mat SY = Scale(Y);
      arma::mat SYlag = Scale(Ylag);
      arma::mat SZ = Scale(Z);
      arma::mat SX = Scale(X);
      int L = SY.n_rows / Omega;
      arma::ivec VarKey = arma::join_cols(arma::zeros<arma::ivec>(1), arma::join_cols(1 * arma::ones<arma::ivec>(SYlag.n_cols), arma::join_cols(2 * arma::ones<arma::ivec>(SZ.n_cols), 3 * arma::ones<arma::ivec>(SX.n_cols))));
      arma::ivec PenaltyKey = arma::join_cols(arma::zeros<arma::ivec>(1 + SYlag.n_cols + SZ.n_cols), arma::ones<arma::ivec>(SX.n_cols));
      arma::uvec FoldKey = arma::shuffle(arma::floor(K * arma::regspace<arma::uvec>(0, L - 1) / L));
      arma::imat Folds = arma::zeros<arma::imat>(L * Omega, K);
      for(int k = 0; k < K; k++) {
        arma::uvec Key = arma::shuffle(arma::find(FoldKey != k));
        int NKey = Key.n_rows;
        arma::uvec TrainKey = Key.rows(0, ceil(4 * NKey / 5) - 1);
        arma::uvec ValidKey = Key.rows(ceil(4 * NKey / 5), NKey - 1);
        for(int rT : TrainKey) {
          Folds(rT, k) = 2;
        }
        for(int rV : ValidKey) {
          Folds(rV, k) = 1;
        }
      }
      if(Omega > 1) {
        for(int w = 1; w < Omega; w++) {
          Folds.rows(w * L, (w + 1) * L - 1) = Folds.rows(0, L - 1);
        }
      }
      arma::mat Input = arma::join_rows(SYlag, arma::join_rows(SZ, SX));
      int IDim = Input.n_cols;
      int ODim = SY.n_cols;
      double Xavier = sqrt(2 / (float(IDim) + float(ODim)));
      arma::cube Wi = Xavier * arma::randn(IDim + 1, HDim, R);
      arma::cube Wh = Xavier * arma::randn(HDim + 1, HDim, R);
      arma::cube Wo = Xavier * arma::randn(HDim + 1, ODim, R);
      arma::mat Final = arma::zeros(SX.n_cols, R);
      int NLambdas = Lambdas.n_rows;
      arma::vec LamOut = arma::zeros(R);
      for(int r = 0; r < R; r++) {
        double Lambda;
        if(NLambdas > 1) {
          arma::mat Estimates = arma::zeros(NLambdas);
          for(int lam = 0; lam < NLambdas; lam++) {
            for(int k = 0; k < K; k++) {
              arma::mat trainInput = Input.rows(arma::find(Folds.col(k) == 2));
              arma::mat trainOutput = SY.rows(arma::find(Folds.col(k) == 2));
              arma::mat validInput = Input.rows(arma::find(Folds.col(k) == 1));
              arma::mat validOutput = SY.rows(arma::find(Folds.col(k) == 1));
              arma::mat testInput = Input.rows(arma::find(Folds.col(k) == 0));
              arma::mat testOutput = SY.rows(arma::find(Folds.col(k) == 0));
              List MLP = TrainMLP(trainInput, trainOutput, validInput, validOutput, Type, HDim, Penalty, Lambdas(lam), VarKey, PenaltyKey, Wi.slice(r), Wh.slice(r), Wo.slice(r), Activation);
              arma::mat fWi = MLP[0];
              arma::mat fWh = MLP[2];
              arma::mat fWo = MLP[3];
              List Forward = ForwardMLP1(testInput, fWi, fWo, Activation);
              arma::mat Hat = Forward[0];
              double Residuals = arma::accu(arma::pow((testOutput - Hat), 2)) / testOutput.n_rows;
              Estimates(lam) = Estimates(lam) + Residuals / K;
            }
          }
          arma::uword Est = Estimates.index_min();
          Lambda = Lambdas(Est);
        } else {
          Lambda = Lambdas(0);
        }
        LamOut(r) = Lambda;
        List fMLP = TrainMLP(Input, SY, Input, SY, Type, HDim, Penalty, Lambda, VarKey, PenaltyKey, Wi.slice(r), Wh.slice(r), Wo.slice(r), Activation);
        arma::mat ffWi = fMLP[0];
        Final.col(r) = arma::sqrt(arma::sum(arma::pow(ffWi.rows(arma::find(PenaltyKey == 1)), 2), 1));
      }
      return List::create(_["W0"]= Final, _["Lambdas"] = LamOut);
    } else {
      arma::mat SY = Scale(Y);
      arma::mat SYlag = Scale(Ylag);
      arma::mat SZ = Scale(Z);
      arma::mat SX = Scale(X);
      arma::ivec VarKey = arma::join_cols(arma::zeros<arma::ivec>(1), arma::join_cols(1 * arma::ones<arma::ivec>(SYlag.n_cols), arma::join_cols(2 * arma::ones<arma::ivec>(SZ.n_cols), 3 * arma::ones<arma::ivec>(SX.n_cols))));
      arma::ivec PenaltyKey = arma::join_cols(arma::zeros<arma::ivec>(1 + SYlag.n_cols + SZ.n_cols), arma::ones<arma::ivec>(SX.n_cols));
      arma::mat Input = arma::join_rows(SYlag, arma::join_rows(SZ, SX));
      int IDim = Input.n_cols;
      int ODim = SY.n_cols;
      double Xavier = sqrt(2 / (float(IDim) + float(ODim)));
      arma::cube Wi = Xavier * arma::randn(IDim + 1, HDim, R);
      arma::cube Wh = Xavier * arma::randn(HDim + 1, HDim, R);
      arma::cube Wo = Xavier * arma::randn(HDim + 1, ODim, R);
      int NLambdas = Lambdas.n_rows;
      arma::cube Final = arma::zeros(SX.n_cols, R, NLambdas);
      for(int lam = 0; lam < NLambdas; lam++) {
        for(int r = 0; r < R; r++) {
          List fMLP = TrainMLP(Input, SY, Input, SY, Type, HDim, Penalty, Lambdas(lam), VarKey, PenaltyKey, Wi.slice(r), Wh.slice(r), Wo.slice(r), Activation);
          arma::mat ffWi = fMLP[0];
          Final.slice(lam).col(r) = arma::sqrt(arma::sum(arma::pow(ffWi.rows(arma::find(PenaltyKey == 1)), 2), 1));
        }
      }
      return List::create(_["W0"]= Final, _["Lambdas"] = Lambdas);
    }
  } else {
    if(ChooseLambda) {
      arma::mat SY = Scale(Y);
      arma::mat SYlag = Scale(Ylag);
      arma::mat SZ = Scale(Z);
      arma::mat SX = Scale(X);
      int L = SY.n_rows / Omega;
      arma::ivec VarKey = arma::join_cols(arma::zeros<arma::ivec>(1), arma::join_cols(1 * arma::ones<arma::ivec>(SYlag.n_cols), arma::join_cols(2 * arma::ones<arma::ivec>(SZ.n_cols), 3 * arma::ones<arma::ivec>(SX.n_cols))));
      arma::ivec PenaltyKey = arma::join_cols(arma::zeros<arma::ivec>(1 + SYlag.n_cols + SZ.n_cols), arma::ones<arma::ivec>(SX.n_cols));
      arma::uvec FoldKey = arma::shuffle(arma::floor(K * arma::regspace<arma::uvec>(0, L - 1) / L));
      arma::imat Folds = arma::zeros<arma::imat>(L * Omega, K);
      for(int k = 0; k < K; k++) {
        arma::uvec Key = arma::shuffle(arma::find(FoldKey != k));
        int NKey = Key.n_rows;
        arma::uvec TrainKey = Key.rows(0, ceil(4 * NKey / 5) - 1);
        arma::uvec ValidKey = Key.rows(ceil(4 * NKey / 5), NKey - 1);
        for(int rT : TrainKey) {
          Folds(rT, k) = 2;
        }
        for(int rV : ValidKey) {
          Folds(rV, k) = 1;
        }
      }
      if(Omega > 1) {
        for(int w = 1; w < Omega; w++) {
          Folds.rows(w * L, (w + 1) * L - 1) = Folds.rows(0, L - 1);
        }
      }
      arma::mat Input = arma::join_rows(SYlag, arma::join_rows(SZ, SX));
      int IDim = Input.n_cols;
      int ODim = SY.n_cols;
      double Xavier = sqrt(2 / (float(IDim) + float(ODim)));
      arma::cube Wi = Xavier * arma::randn(IDim + 1, HDim, R);
      arma::cube Wh = Xavier * arma::randn(HDim + 1, HDim, R);
      arma::cube Wo = Xavier * arma::randn(HDim + 1, ODim, R);
      arma::mat Final = arma::zeros(SX.n_cols, R);
      int NLambdas = Lambdas.n_rows;
      arma::vec LamOut = arma::zeros(R);
      for(int r = 0; r < R; r++) {
        double Lambda;
        if(NLambdas > 1) {
          arma::mat Estimates = arma::zeros(NLambdas);
          for(int lam = 0; lam < NLambdas; lam++) {
            for(int k = 0; k < K; k++) {
              arma::mat trainInput = Input.rows(arma::find(Folds.col(k) == 2));
              arma::mat trainOutput = SY.rows(arma::find(Folds.col(k) == 2));
              arma::mat validInput = Input.rows(arma::find(Folds.col(k) == 1));
              arma::mat validOutput = SY.rows(arma::find(Folds.col(k) == 1));
              arma::mat testInput = Input.rows(arma::find(Folds.col(k) == 0));
              arma::mat testOutput = SY.rows(arma::find(Folds.col(k) == 0));
              List MLP = TrainMLP(trainInput, trainOutput, validInput, validOutput, Type, HDim, Penalty, Lambdas(lam), VarKey, PenaltyKey, Wi.slice(r), Wh.slice(r), Wo.slice(r), Activation);
              arma::mat fWi = MLP[0];
              arma::mat fWh = MLP[2];
              arma::mat fWo = MLP[3];
              List Forward = ForwardMLP2(testInput, fWi, fWh, fWo, Activation);
              arma::mat Hat = Forward[0];
              double Residuals = arma::accu(arma::pow((testOutput - Hat), 2)) / testOutput.n_rows;
              Estimates(lam) = Estimates(lam) + Residuals / K;
            }
          }
          arma::uword Est = Estimates.index_min();
          Lambda = Lambdas(Est);
        } else {
          Lambda = Lambdas(0);
        }
        LamOut(r) = Lambda;
        List fMLP = TrainMLP(Input, SY, Input, SY, Type, HDim, Penalty, Lambda, VarKey, PenaltyKey, Wi.slice(r), Wh.slice(r), Wo.slice(r), Activation);
        arma::mat ffWi = fMLP[0];
        Final.col(r) = arma::sqrt(arma::sum(arma::pow(ffWi.rows(arma::find(PenaltyKey == 1)), 2), 1));
      }
      return List::create(_["W0"]= Final, _["Lambdas"] = LamOut);
    } else {
      arma::mat SY = Scale(Y);
      arma::mat SYlag = Scale(Ylag);
      arma::mat SZ = Scale(Z);
      arma::mat SX = Scale(X);
      arma::ivec VarKey = arma::join_cols(arma::zeros<arma::ivec>(1), arma::join_cols(1 * arma::ones<arma::ivec>(SYlag.n_cols), arma::join_cols(2 * arma::ones<arma::ivec>(SZ.n_cols), 3 * arma::ones<arma::ivec>(SX.n_cols))));
      arma::ivec PenaltyKey = arma::join_cols(arma::zeros<arma::ivec>(1 + SYlag.n_cols + SZ.n_cols), arma::ones<arma::ivec>(SX.n_cols));
      arma::mat Input = arma::join_rows(SYlag, arma::join_rows(SZ, SX));
      int IDim = Input.n_cols;
      int ODim = SY.n_cols;
      double Xavier = sqrt(2 / (float(IDim) + float(ODim)));
      arma::cube Wi = Xavier * arma::randn(IDim + 1, HDim, R);
      arma::cube Wh = Xavier * arma::randn(HDim + 1, HDim, R);
      arma::cube Wo = Xavier * arma::randn(HDim + 1, ODim, R);
      int NLambdas = Lambdas.n_rows;
      arma::cube Final = arma::zeros(SX.n_cols, R, NLambdas);
      for(int lam = 0; lam < NLambdas; lam++) {
        for(int r = 0; r < R; r++) {
          List fMLP = TrainMLP(Input, SY, Input, SY, Type, HDim, Penalty, Lambdas(lam), VarKey, PenaltyKey, Wi.slice(r), Wh.slice(r), Wo.slice(r), Activation);
          arma::mat ffWi = fMLP[0];
          Final.slice(lam).col(r) = arma::sqrt(arma::sum(arma::pow(ffWi.rows(arma::find(PenaltyKey == 1)), 2), 1));
        }
      }
      return List::create(_["W0"]= Final, _["Lambdas"] = Lambdas);
    }
  }
}