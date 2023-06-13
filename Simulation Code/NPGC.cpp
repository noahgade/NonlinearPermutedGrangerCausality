#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat scaleCpp(arma::mat input) {
  float dim = input.n_cols;
  for(int i = 0; i < dim; i++) {
    input.col(i) = (input.col(i) - arma::mean(input.col(i))) / arma::stddev(input.col(i));
  }
  return input;
}

// [[Rcpp::export]]
arma::mat lagmatCpp(arma::mat input, float lag) {
  float dim = input.n_cols;
  float out_length = input.n_rows - lag;
  arma::mat output = arma::mat(out_length, dim * lag);
  for(int i = 0; i < dim; i++) {
    for(int j = 0; j < lag; j++) {
      output.col(j + i * lag) = input.submat(lag - 1 - j, i, lag - 2 - j + out_length, i);
    }
  }
  return output;
}

// [[Rcpp::export]]
arma::cube generateW(float nrows, float ncols, float n_initializations) {
  arma::cube W = arma::randn(nrows, ncols, n_initializations);
  return W;
}

// [[Rcpp::export]]
arma::cube npgc(arma::mat Y, arma::mat X, arma::mat Ylag, arma::mat Z, float n_permutations, float n_folds, float n_initializations, float feature_dimension) {
  float L = Y.n_rows;
  arma::mat A = scaleCpp(arma::join_horiz(Ylag, Z));
  arma::mat A0 = arma::join_horiz(arma::mat(L, 1, arma::fill::ones), A);
  arma::mat X0 = scaleCpp(X);
  
  arma::cube Xp = arma::cube(X0.n_rows, X0.n_cols, n_permutations);
  Xp.slice(0) = X0;
  for(int permutation = 1; permutation < n_permutations; permutation++) {
    arma::uvec pshuffle = arma::randperm(floor(L));
    Xp.slice(permutation) = X0.rows(pshuffle);
  }
  
  arma::vec folds = ceil(arma::linspace(0, Y.n_rows, n_folds + 1));
  arma::uvec fshuffle = arma::randperm(floor(L));
  
  arma::cube output = arma::cube(n_folds, n_initializations, n_permutations); 
  
  for(int fold = 0; fold < n_folds; fold++) {
    arma::uvec test_index = fshuffle.rows(folds(fold), folds(fold + 1) - 1);
    arma::uvec train_index = arma::uvec(L - test_index.n_rows);
    
    if(fold == 0) {
      train_index = fshuffle.rows(folds(1), folds(n_folds) - 1);
    } else if(fold == (n_folds - 1)) {
      train_index = fshuffle.rows(folds(0), folds(n_folds - 1) - 1);
    } else {
      train_index = arma::join_vert(fshuffle.rows(folds(0), folds(fold) - 1), fshuffle.rows(folds(fold + 1), folds(n_folds) - 1));
    }
    
    arma::cube W = generateW(X0.n_cols + A0.n_cols, feature_dimension, n_initializations);
    
    for(int permutation = 0; permutation < n_permutations; permutation++) {
      for(int initialization = 0; initialization < n_initializations; initialization++) {
        arma::mat Htest = tanh(arma::join_horiz(A0.rows(test_index), Xp.slice(permutation).rows(test_index)) * W.slice(initialization));
        arma::mat Htrain = tanh(arma::join_horiz(A0.rows(train_index), Xp.slice(permutation).rows(train_index)) * W.slice(initialization));
        arma::mat Resid = Y.rows(test_index) - Htest * (arma::inv(trans(Htrain) * Htrain) * trans(Htrain) * Y.rows(train_index));
        output(fold, initialization, permutation) = arma::trace(trans(Resid) * Resid / test_index.n_rows);
      }
    }
  }
  
  return output;
}

// [[Rcpp::export]]
arma::mat naive(arma::mat Y, arma::mat X, arma::mat Ylag, arma::mat Z, float n_initializations, float feature_dimension) {
  float L = Y.n_rows;
  arma::mat A = scaleCpp(arma::join_horiz(Ylag, Z));
  arma::mat A0 = arma::join_horiz(arma::mat(L, 1, arma::fill::ones), A);
  arma::mat X0 = scaleCpp(X);
  
  arma::mat Xres = arma::mat(X0.n_rows, X0.n_cols, arma::fill::zeros);
  
  arma::cube W = generateW(X0.n_cols + A0.n_cols, feature_dimension, n_initializations);
  arma::mat output = arma::mat(3, n_initializations);
    
  for(int initialization = 0; initialization < n_initializations; initialization++) {
    arma::mat Xsub = scaleCpp(arma::randn(X0.n_rows, X0.n_cols));
    arma::mat Hunres = tanh(arma::join_horiz(A0, X0) * W.slice(initialization));
    arma::mat Hres = tanh(arma::join_horiz(A0, Xres) * W.slice(initialization));
    arma::mat Hsub = tanh(arma::join_horiz(A0, Xsub) * W.slice(initialization));
    
    arma::mat Runres = Y - Hunres * (arma::inv(trans(Hunres) * Hunres) * trans(Hunres) * Y);
    arma::mat Rres = Y - Hres * (arma::inv(trans(Hres) * Hres) * trans(Hres) * Y);
    arma::mat Rsub = Y - Hsub * (arma::inv(trans(Hsub) * Hsub) * trans(Hsub) * Y);
    output(0, initialization) = arma::trace(trans(Runres) * Runres / Runres.n_rows);
    output(1, initialization) = arma::trace(trans(Rres) * Rres / Rres.n_rows);
    output(2, initialization) = arma::trace(trans(Rsub) * Rsub / Rsub.n_rows);
  }
  
  return output;
}

// [[Rcpp::export]]
arma::cube application_npgc(arma::mat Y, arma::mat X, arma::mat Ylag, arma::mat Z, float n_permutations, float n_initializations, float feature_dimension, arma::cube W1, arma::cube W2, arma::cube W3, arma::cube W4, arma::cube W5) {
  float L = Y.n_rows;
  float n_folds = 5;
  arma::mat A = scaleCpp(arma::join_horiz(Ylag, Z));
  arma::mat A0 = arma::join_horiz(arma::mat(L, 1, arma::fill::ones), A);
  arma::mat X0 = scaleCpp(X);
  
  arma::cube Xp = arma::cube(X0.n_rows, X0.n_cols, n_permutations);
  Xp.slice(0) = X0;
  for(int permutation = 1; permutation < n_permutations; permutation++) {
    arma::uvec pshuffle = arma::randperm(floor(L));
    Xp.slice(permutation) = X0.rows(pshuffle);
  }
  
  arma::vec folds = ceil(arma::linspace(0, Y.n_rows, n_folds + 1));
  arma::uvec fshuffle = arma::randperm(floor(L));
  
  arma::cube output = arma::cube(n_folds, n_initializations, n_permutations); 
  arma::cube W = arma::cube(X0.n_cols + A0.n_cols, feature_dimension, n_initializations);
  
  for(int fold = 0; fold < n_folds; fold++) {
    arma::uvec test_index = fshuffle.rows(folds(fold), folds(fold + 1) - 1);
    arma::uvec train_index = arma::uvec(L - test_index.n_rows);
    
    if(fold == 0) {
      train_index = fshuffle.rows(folds(1), folds(n_folds) - 1);
    } else if(fold == (n_folds - 1)) {
      train_index = fshuffle.rows(folds(0), folds(n_folds - 1) - 1);
    } else {
      train_index = arma::join_vert(fshuffle.rows(folds(0), folds(fold) - 1), fshuffle.rows(folds(fold + 1), folds(n_folds) - 1));
    }
    
    if(fold == 0) {
      W = W1;
    } else if(fold == 1) {
      W = W2;
    } else if(fold == 2) {
      W = W3;
    } else if(fold == 3) {
      W = W4;
    } else {
      W = W5;
    }
    
    for(int permutation = 0; permutation < n_permutations; permutation++) {
      for(int initialization = 0; initialization < n_initializations; initialization++) {
        arma::mat Htest = tanh(arma::join_horiz(A0.rows(test_index), Xp.slice(permutation).rows(test_index)) * W.slice(initialization));
        arma::mat Htrain = tanh(arma::join_horiz(A0.rows(train_index), Xp.slice(permutation).rows(train_index)) * W.slice(initialization));
        arma::mat Resid = Y.rows(test_index) - Htest * (arma::inv(trans(Htrain) * Htrain) * trans(Htrain) * Y.rows(train_index));
        output(fold, initialization, permutation) = arma::trace(trans(Resid) * Resid / test_index.n_rows);
      }
    }
  }
  
  return output;
}