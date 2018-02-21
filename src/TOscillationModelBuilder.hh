///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Oct, 2016
// Built on top of macros created by K. Gilje
// TODO: Implement the ability to create number of bins on the inputs
//   Implement error level
//   Implement a way to include information from more than 1 detector location
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TOSCILLATIONMODELBUILDER_HH
#define TOSCILLATIONMODELBUILDER_HH

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

//Include Osc_Sens related headers
#include "TExperiment.hh"
#include "TOscillationSimulator.hh"
#include "TDetectorResponseInterface.hh"

/// Oscillation model builder Class
/**
 * Class that corresponds to the oscillation models built for comparison to data/toys/models
 */
class TOscillationModelBuilder
{
public:
  
  /// Default Constructor
  TOscillationModelBuilder(TExperiment &);
  
  /// Default destructor
  ~TOscillationModelBuilder();
  
  ///////////////////////////////////////////////////////////////////////
  //public attributes
  ///////////////////////////////////////////////////////////////////////
  
  /// TExperiment object this class uses for some information
  TExperiment &fExperiment;
  
  /// Oscillation simulator for simulating oscillations
  TOscillationSimulator fOscillationSimulator;
  
  /// Number of bins to be used for deltam2
  int fNDeltam2;
  
  /// Number of bins to be used for sin22theta
  int fNSinSq2Theta;
  
  /// Position tree populated in TDetector
  std::map<int, TTree*> fPositionTree;
  
  // vectors and arrays for sin22theta and deltam2, there got to be a better way to do this.
  /// Vector of deltam2
  std::vector<double> fDeltam2;
  /// Vector of Sin22theta
  std::vector<double> fSinSq2Theta;
  /// deltam2 bins
  double* fDeltam2Bins;
  /// Sin22theta bins
  double* fSinSq2ThetaBins;
  
  /// Simulated neutrino energy
  std::map<int, TH1D*> hSimulatedNueEnergy;
  
  /// LvsE histogram for all the values of deltam2  corresponding to #fDeltam2
  std::map<int, TH2D*> hLvsEDeltam2;
  
  /// LvsE histogram for all the values of deltam2  corresponding to #fDeltam2
  std::map<int, TH2D*> hLvsEDeltam2Relative;
  
  ///////////////////////////////////////////////////////////////////////
  //public member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// Obtain sin2theta22 and deltam2 bin numbers
  int EstimateBinNumbers(double, double, int&, int&);
  /// Setup the parameters(deltam2 and sin22theta) for oscillations
  /** The input has to be the number deltam2 and sin22theta values to be used for oscillation in that particular order
   */
  void SetUpModelOscillation(int nDeltam2=57,int nSinSq2Theta=50);
  
  /// Simulates oscillation.
  /**
   Simulates oscillation for the given experiment, detector, the histogram to be filled a given value deltam2.
   */
  void SimulateModelOscillation(const TExperiment&,TDetector&,TH2D &, double);
  
  /// Simulates oscillation.
  /**
   Simulates oscillation for the given experiment and the given values of deltam2.
   */
  void SimulateModelOscillation();
  
  /// Write histograms corresponding to oscillation
  void WriteHistograms(TFile & );
  
private:
  
  //EDIT: Remove deltam2bin and sin22thetabin dynamic arrays
  /// Construct bins, used in creation of chi2 histograms
  void ConstructDeltam2Bins();
  
  /// Construct the values of deltam2 based on #fDeltam2
  void ConstructDeltam2();
  
  /// Construct bins, used in creation of chi2 histograms
  void ConstructSinSq2ThetaBins();
  
  /// Construct the values of sin22theta based on #fSinSq2ThetaBins
  void ConstructSinSq2Theta();
  
  /// Construct deltam2 and sin22theta bins from the input file
  /*
   The input file should either be a .txt file or .root file
   */
  void ConstructBinsFromFile(TString);
  
  /// Setup the histograms corresponding to each value of the oscillation parameter #fDeltam2
  void SetUpHistograms();
};

#endif
