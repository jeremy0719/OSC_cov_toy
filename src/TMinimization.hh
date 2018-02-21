///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: June 23, 2016
// TODO: Implement the ability to create number of bins on the inputs
//   Implement error level
//   In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TMINIMIZATION_HH
#define TMINIMIZATION_HH

//Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

//Include ROOT headers
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"

//Include Osc_Sens related headers
//#include "TOscillationSimulator.hh"
#include "TOscillationModelBuilder.hh"
#include "TExperiment.hh"
#include "TCovarianceMatrixInterface.hh"
#include "TOscSensUtils.hh"

#include <time.h>

/// Minimization parent class (Sensitivity estimation and Oscillation fitting classes are derived from this class)
/*
 In future there could be a daugther class TCovarianceMinimization that would be a minimizer but a covariance matrix based minimizer.
 */
class TMinimization
{
public:

  /// Default minimization constructor
  //TMinimization(){}
  
  /// Minimization constructor, takes an OscillationModelBuilder as an input
  TMinimization(TOscillationModelBuilder &);
  
  /// Default
  virtual ~TMinimization();
  
  ///////////////////////////////////////////////////////////////////////
  //public attributes
  ///////////////////////////////////////////////////////////////////////
  
  /// TOscillationModelBuilder object as data member
  TOscillationModelBuilder &fOscillationModelBuilder;
  
  /// TExperiment object as a data member
  TExperiment &fExperiment;
  
  // TCovarianceMatrixInterface object as a datamember
  TCovarianceMatrixInterface fCovMatInt;
  
  /// Number of deltam2 bins
  double fNDeltam2;
  /// Number of sin22theta bins
  double fNSinSq2Theta;

  /// Number of energy bins, common for all the detectors.
  int fEneBins=0;
  
  /// Reference histograms to be compared for calculating chi2
  std::map<int, TH2D*> hReferenceLvsE;
  
  // The keys for all the following maps
  /// Map of null-oscillated signal LvsE
  std::map<int,TH2D*> hLvsENull;
  
  /// The final output chisquare map, separate for each detector
  std::map<int,TH2D*> hChiSquareMap;
  
  /// The final output chisquare map for the combination of the detectors
  TH2D* hCumulativeChiSquareMap;
  
  /// Deltam2 histograms
  std::map<int,TH2D*> hLvsEDeltam2;
  
  ///////////////////////////////////////////////////////////////////////
  //public member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// Function to setup minimizer
  /*
   Setup covariance matrices and chisquare maps
   */
  virtual void SetupMinimizer();
  
  /* Minimize, all the functinality is in <BuildChiSquareMap>"()" and <BuildCovarianceMatrices>"()"*/
  virtual void Minimize();
  
  /// Write the minimization related historgams to the file
  void WriteHistograms(TFile &);
  
  const TH2D* GetCumulativeChiSquareMap() const{return hCumulativeChiSquareMap;}
  
  void SetCumulativeChiSquareMapName(TString CumulativeChiSquareMapName) {hCumulativeChiSquareMap->SetName(CumulativeChiSquareMapName.Data());}
  
private:
  
  ///////////////////////////////////////////////////////////////////////
  //private attributes
  ///////////////////////////////////////////////////////////////////////
  
  // Detector on and off fractions
//  double detOnFraction = 0.0;
//  double detOffFraction = 0.0;
  // Boolean to tell if the chi2 should be calculated using realtive histograms or absolute historgams
  bool doRelativeMinimization=false;
  
  ///////////////////////////////////////////////////////////////////////
  //private member functions
  ///////////////////////////////////////////////////////////////////////
    
  /// Calculate chi2 for one particular value of deltam2 and sin22theta, returns the calculated value
  double CalculateChiSquareValue(TDetector&, double, double);
  
  /// Build Chisquare maps after minimization.
  void BuildChiSquareMap();
  
  /// Clear all the minimizer histograms
  void ClearChi2Histograms();

};

#endif
