///////////////////////////////////////////////////////////////////////
// Author: K Gilje
// Date: October 19, 2016
// TODO:
//   Implement error level
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TTHROWMCTOY_HH
#define TTHROWMCTOY_HH


// Include standard headers
#include <iostream>
#include <algorithm>
#include <vector>

// Include ROOT headers
#include <TMath.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <TDecompChol.h>
#include <TRandom3.h>
#include <TF1.h>

class TThrowMCToy {
private:
  
  ///////////////////////////////////////////////////////////////////////
  // private member attributes
  ///////////////////////////////////////////////////////////////////////
  
  
  /// The number of parameters
  int nParameters;

  /// A vector of the output parameters
  TVectorD *v_Parameters;

  /// The covariance matrix
  TMatrixDSym *m_Covariance;

  /// The Cholesky Decomposed matrix
  TMatrixD *m_CholeskyDecomp;

  /// The random number generator
  TRandom3 randomGen;

  /// A boolean related to the success of the decomposition
  bool status = true;
  
  ///////////////////////////////////////////////////////////////////////
  //private member functions
  ///////////////////////////////////////////////////////////////////////

  /// A function to decompose a matrix (m_Covariance)
  void CholeskyDecompose(TMatrixD &m_CholDec);
  
  /// A function to create a random pair of numbers using the Box-Muller polar form
  void RandomGaussian(double *z);

public:
  /// Default Constructor
  TThrowMCToy();
  
  /// Constructor: Takes a vector of parameters (LE-bin values) and a covariance matrix as inputs
  TThrowMCToy(TVectorD &v_inputParameters, TMatrixDSym &m_inputCovariance);

  /// Destructor
  ~TThrowMCToy();
  
  ///////////////////////////////////////////////////////////////////////
  // public member functions
  ///////////////////////////////////////////////////////////////////////

  /// Set the seed of the random generator
  void SetSeed(int seed = 1867) {randomGen.SetSeed(seed);};
  
  /// Get the status of the decomposition
  bool GetStatus() {return status;};

  /// Get the number of input parameters
  int GetSize() {return nParameters;};
  
  /// Set the number of input parameters
  void SetSize(int size) {nParameters = size;};
  
  /// Get the parameters
  TVectorD* GetParameters() {return v_Parameters;};
  
  /// Set the parameters
  void SetParameters(TVectorD &parameters) {v_Parameters = new TVectorD(parameters);};

  /// Get the covariance matrix
  TMatrixDSym* GetCovariance() {return m_Covariance;};
  
  /// Set the covariance matrix
  void SetCovariance(TMatrixDSym &covariance) {m_Covariance = new TMatrixDSym(covariance);};
  
  /// Call to throw a single experiment
  void ThrowExperiment(std::vector<double> &thrownParameters);

};

#endif
