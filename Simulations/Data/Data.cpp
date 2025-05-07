#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Functions to generate simulated VAR data

arma::mat VARMat(int Dim, int Lag, arma::mat Causal) {
  // Generates VAR transition matrices
  // Dim: dimension of the data
  // Lag: lags to include, a VAR(lag) model
  // Causal: Dim x (Dim * Lag) cube of ones and zeros for each regime
  double MaxSpecRad = 0.8;
  double MinCausal = 0.2;
  arma::mat A = arma::zeros(Dim * Lag, Dim * Lag);
  bool OuterGoal = false;
  int OuterIter = 0;
  while(not OuterGoal) {
    bool StartOver = false;
    A.rows(0, Dim - 1) = (2 * arma::randu(Dim, Dim * Lag) - 1.0) % Causal;
    if(Lag > 1) {
      A.rows(Dim, Dim * Lag - 1).cols(0, Dim * (Lag - 1) - 1) = arma::eye(Dim * (Lag - 1), Dim * (Lag - 1));
    }
    bool InnerGoal0 = false;
    int InnerIter0 = 0;
    while(not InnerGoal0) {
      double Eigen = sort(abs(eig_gen(A))).max();
      A.rows(0, Dim - 1) = MaxSpecRad * A.rows(0, Dim - 1) / Eigen;
      for(int i = 0; i < Dim; i++) {
        for(int j = 0; j < Dim * Lag; j++) {
          if(Causal(i, j) != 0) {
            if((A(i, j) > 0) and (A(i, j) < MinCausal)) {
              A(i, j) = MinCausal;
            }
            if((A(i, j) < 0) and (A(i, j) > -1 * MinCausal)) {
              A(i, j) = -1 * MinCausal;
            }
          }
        }
      }
      double Eigen0 = sort(abs(eig_gen(A))).max();
      InnerIter0 = InnerIter0 + 1;
      double ME0 = abs(Eigen0 - MaxSpecRad);
      if((ME0 < 0.01) and (Eigen0 < 1)) {
        InnerGoal0 = true;
      }
      if((Eigen0 == Eigen) or (InnerIter0 > 100000)) {
        StartOver = true;
        break;
      }
    }
    if(not StartOver) {
      OuterGoal = true;
    }
    OuterIter = OuterIter + 1;
    if(OuterIter > 100) {
      MinCausal = MinCausal / 10;
    }
  }
  return A;
}

// [[Rcpp::export]]
List VAR(int L, int Dim, int Lag, arma::mat Causal) {
  // Vector autoregressive data
  // L: number of observations in the data
  // Dim: number of variables in the data
  // Lag: lags to include, a VAR(lag) model
  // Causal: Dim x (Dim * Lag) matrix of ones and zeros
  // MaxSpecRad: Maximum spectral radius of the VAR matrix
  // MinCausal: Minimum value of "Causal" element
  int BurnIn = 100;
  arma::mat tOutput = arma::zeros(L + BurnIn, Dim * Lag);
  arma::mat Innovations = arma::join_rows(arma::randn(L + BurnIn, Dim), arma::zeros(L + BurnIn, Dim * (Lag - 1)));
  arma::mat A = VARMat(Dim, Lag, Causal);
  for(int time = 1; time < (L + BurnIn); time++) {
    tOutput.row(time) = tOutput.row(time - 1) * trans(A) + Innovations.row(time);
  }
  arma::mat Output = tOutput.rows(BurnIn, L + BurnIn - 1).cols(0, Dim - 1);
  return::List::create(_["data"] = Output, _["VAR"] = A, _["Causal"] = Causal);
}

// Functions to generate simulated TAR data

List TARMat(int Dim, int Lag, arma::cube Causal) {
  // Generates TAR transition matrices
  // Dim: dimension of the data
  // Lag: lags to include, a VAR(lag) model
  // Causal: Dim x (Dim * Lag) cube of ones and zeros for each regime
  double MaxSpecRad = 0.8;
  double MinCausal = 0.2;
  arma::cube A = arma::zeros(Dim * Lag, Dim * Lag, 2);
  bool OuterGoal = false;
  int OuterIter = 0;
  while(not OuterGoal) {
    bool StartOver = false;
    A.slice(0).rows(0, Dim - 1) = (2 * arma::randu(Dim, Dim * Lag) - 1.0) % Causal.slice(0);
    A.slice(1).rows(0, Dim - 1) = (2 * arma::randu(Dim, Dim * Lag) - 1.0) % Causal.slice(1);
    if(Lag > 1) {
      A.slice(0).rows(Dim, Dim * Lag - 1).cols(0, Dim * (Lag - 1) - 1) = arma::eye(Dim * (Lag - 1), Dim * (Lag - 1));
      A.slice(1).rows(Dim, Dim * Lag - 1).cols(0, Dim * (Lag - 1) - 1) = arma::eye(Dim * (Lag - 1), Dim * (Lag - 1));
    }
    bool InnerGoal0 = false;
    int InnerIter0 = 0;
    while(not InnerGoal0) {
      double Eigen = sort(abs(eig_gen(A.slice(0)))).max();
      A.slice(0).rows(0, Dim - 1) = MaxSpecRad * A.slice(0).rows(0, Dim - 1) / Eigen;
      for(int i = 0; i < Dim; i++) {
        for(int j = 0; j < Dim * Lag; j++) {
          if(Causal.slice(0)(i, j) != 0) {
            if((A.slice(0)(i, j) > 0) and (A.slice(0)(i, j) < MinCausal)) {
              A.slice(0)(i, j) = MinCausal;
            }
            if((A.slice(0)(i, j) < 0) and (A.slice(0)(i, j) > -1 * MinCausal)) {
              A.slice(0)(i, j) = -1 * MinCausal;
            }
          }
        }
      }
      double Eigen0 = sort(abs(eig_gen(A.slice(0)))).max();
      InnerIter0 = InnerIter0 + 1;
      double ME0 = abs(Eigen0 - MaxSpecRad);
      if((ME0 < 0.01) and (Eigen0 < 1)) {
        InnerGoal0 = true;
      }
      if((Eigen0 == Eigen) or (InnerIter0 > 100000)) {
        StartOver = true;
        break;
      }
    }
    bool InnerGoal1 = false;
    int InnerIter1 = 0;
    while(not InnerGoal1) {
      double Eigen = sort(abs(eig_gen(A.slice(1)))).max();
      A.slice(1).rows(0, Dim - 1) = MaxSpecRad * A.slice(1).rows(0, Dim - 1) / Eigen;
      for(int i = 0; i < Dim; i++) {
        for(int j = 0; j < Dim * Lag; j++) {
          if(Causal.slice(1)(i, j) != 0) {
            if((A.slice(1)(i, j) > 0) and (A.slice(1)(i, j) < MinCausal)) {
              A.slice(1)(i, j) = MinCausal;
            }
            if((A.slice(1)(i, j) < 0) and (A.slice(1)(i, j) > -1 * MinCausal)) {
              A.slice(1)(i, j) = -1 * MinCausal;
            }
          }
        }
      }
      double Eigen0 = sort(abs(eig_gen(A.slice(1)))).max();
      InnerIter1 = InnerIter1 + 1;
      double ME0 = abs(Eigen0 - MaxSpecRad);
      if((ME0 < 0.01) and (Eigen0 < 1)) {
        InnerGoal1 = true;
      }
      if((Eigen0 == Eigen) or (InnerIter1 > 100000)) {
        StartOver = true;
        break;
      }
    }
    double Cond = abs(eig_gen(0.5 * arma::kron(A.slice(0), A.slice(0)) + 0.5 * arma::kron(A.slice(1), A.slice(1)))).max();
    if(Cond >= 1) {
      StartOver = true;
    }
    if(not StartOver) {
      OuterGoal = true;
    }
    OuterIter = OuterIter + 1;
    if(OuterIter > 100) {
      MinCausal = MinCausal / 10;
    }
  }
  arma::cube TAR0 = arma::cube(Dim * Lag, Dim * Lag, 1);
  TAR0.slice(0) = A.slice(0);
  arma::cube TAR1 = arma::cube(Dim * Lag, Dim * Lag, 1);
  TAR1.slice(0) = A.slice(1);
  return List::create(_["TAR0"] = TAR0, _["TAR1"] = TAR1);
}

// [[Rcpp::export]]
List TAR(int L, int Dim, int Lag, int RegimeL, arma::cube Causal) {
  // Threshold autoregressive data
  // L: number of observations in the data
  // Dim: number of variables in the data
  // Lag: lags to include, a VAR(lag) model
  // RegimeL: Length of regime dependence, sum of previous RegimeL observations of Dim = 0, < 0 yields Regime 0, otherwise Regime 1
  // Causal: Dim x (Dim * Lag) matrix of ones and zeros
  // MaxSpecRad: Maximum spectral radius of the VAR matrix
  // MinCausal: Minimum value of "Causal" element
  int BurnIn = 100;
  arma::vec Regimes = arma::zeros(L + BurnIn);
  arma::mat tOutput = arma::zeros(L + BurnIn + RegimeL, Dim * Lag);
  arma::mat Innovations = arma::join_rows(arma::randn(L + BurnIn, Dim), arma::zeros(L + BurnIn, Dim * (Lag - 1)));
  List TARMats = TARMat(Dim, Lag, Causal);
  arma::cube A = TARMats[0];
  arma::cube B = TARMats[1];
  for(int time = 0; time < (L + BurnIn); time++) {
    double RegimeSum = arma::accu(tOutput.rows(time, time + RegimeL - 1).col(0));
    arma::mat T = A.slice(0);
    int Regime = 0;
    if(RegimeSum >= 0) {
      Regime = 1;
      T = B.slice(0);
    }
    Regimes(time) = Regime;
    tOutput.row(time + RegimeL) = tOutput.row(time + RegimeL - 1) * trans(T) + Innovations.row(time);
  }
  arma::mat Output = tOutput.rows(RegimeL + BurnIn, L + RegimeL + BurnIn - 1).cols(0, Dim - 1);
  Regimes = Regimes.rows(BurnIn, L + BurnIn - 1);
  return::List::create(_["data"] = Output, _["TAR0"] = A, _["TAR1"] = B, _["Causal"] = Causal, _["Regimes"] = Regimes);
}