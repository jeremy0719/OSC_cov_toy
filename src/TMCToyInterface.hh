///////////////////////////////////////////////////////////////////////
// Author: K Gilje
// Date: October 20, 2016
// TODO:
//   Implement error level
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TMCTOYINTERFACE_HH
#define TMCTOYINTERFACE_HH


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

// Include package headers
#include "TExperiment.hh"
#include "TCovarianceMatrixInterface.hh"
#include "TThrowMCToy.hh"

class TMCToyInterface {
  
private:
  ///////////////////////////////////////////////////////////////////////
  // private member attributes
  ///////////////////////////////////////////////////////////////////////
  
  /// Instance of TExperiment
  TExperiment &fExperiment;
  
  // TCovarianceMatrixInterface object as a datamember
  TCovarianceMatrixInterface fCovMatInt;
  
  /// Vector of reactor on nominal signal
  TVectorD *v_ReactorOnSignal;
  
  /// Vector of reactor on nominal background
  TVectorD *v_ReactorOnBackground;
  
  /// Vector of reactor off nominal background
  TVectorD *v_ReactorOffBackground;
  
  /// The MC Toy thrower for reactor on signal
  TThrowMCToy *toyReactorOnSignal;
  
  /// The MC Toy thrower for reactor on background
  TThrowMCToy *toyReactorOnBackground;
  
  /// The MC Toy thrower for reactor off background
  TThrowMCToy *toyReactorOffBackground;
  
  /// The MC Toy vector for reactor on signal
  std::vector<double> v_RxOnSignalToy;
  
  /// The MC Toy vector for reactor on background
  std::vector<double> v_RxOnBackgroundToy;
  
  /// The MC toy vector for reactor off background
  std::vector<double> v_RxOffBackgroundToy;
  
  /// Reactor On L vs E histograms
  std::map<int, TH2D*> hReactorOn;
  
  /// Reactor Off L vs E histograms
  std::map<int, TH2D*> hReactorOff;
  
  /// Map of reduced covariance matrices from files for correlated signal
  std::map<int, TH2D*> hRCovFCSigMatrix;
  
  /// Map of reduced covariance matrices from files for uncorrelated signal
  std::map<int, TH2D*> hRCovUCSigMatrix;
  
  /// Map of reduced covariance matrices from files for correlated background
  std::map<int, TH2D*> hRCovFCBkgMatrix;
  
  /// Map of reduced covariance matrices from files for uncorrelated background
  std::map<int, TH2D*> hRCovUCBkgMatrix;
  
  ///////////////////////////////////////////////////////////////////////
  // private member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// A function to get one set of toys
  void GetToy();
  
  /// A function to make the toy output appear like data
  void CreateFakeData();
  
  /// A function to take in a map of histograms for histogram bounds and a vector of new values to fill said map with.
  std::map<int, TH2D*> FillOutputHistogram(std::map<int, TH2D*>, TVectorD, TString);
  
  /// A function to take the input map of TH2D histograms and turns them into TVectorD
  TVectorD* ExtractVector(std::map<int, TH2D*>);
  
  /// A function to pull the covariance matrices from files
  //void SetupCovarianceMatrices(int);
  
  /// Build up the total covariance matrix
  ///void BuildCovarianceMatrices();
  
public:
  /// Default Constructor
  TMCToyInterface(TExperiment &);
  
  /// Default destructor
  ~TMCToyInterface();
  
  ///////////////////////////////////////////////////////////////////////
  // public member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// The public function that will construct a single toy that looks like data.
  void ThrowToy();
  
  /// A function to write the output reactor on and reactor off histograms to an output file.
  void WriteHistograms(TFile &);
  
  /// A function that sets the random seed to use with the toy thrower
  void SetSeed(int);
  
  const std::map<int, TH2D*> & GetReactorOnHists() const{
    return hReactorOn;
  }
  
  const std::map<int, TH2D*> & GetReactorOffHists() const{
    return hReactorOff;
  }
  
};
#endif
