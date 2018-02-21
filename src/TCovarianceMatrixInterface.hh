///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Nov, 2016
// TODO:
//   Implement error level
//   Ability to build covriance matrices from list supplied in the macro file
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TCOVARIANCEMATRIXINTERFACE_HH
#define TCOVARIANCEMATRIXINTERFACE_HH

//Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <bitset>

//Include ROOT headers
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

//Include Osc_Sens related headers
#include "TExperiment.hh"
#include "TOscSensUtils.hh"


 ///Interface class to make osc_Sens-ready covariance TMatrix objects from the already existing covariance matrix ROOT files
 /*Note: This class is not a singleton class, i.e., a new object is created everytime you call the constructor. This is to make sure that the class and it's functionality can be used in multiple locations (for example for minimization and for creating toys)  ans multiple times.
 */
class TCovarianceMatrixInterface
{
public:
  
  /// Default constructor
  TCovarianceMatrixInterface(TExperiment &, bool);
  
  /// Default destructor
  ~TCovarianceMatrixInterface();
  
  ///////////////////////////////////////////////////////////////////////
  //public member functions
  ///////////////////////////////////////////////////////////////////////
  /// Write histograms
  void WriteHistograms(TFile&);
  
  
  /// Function to setup covariance matrices
  /*
   Called once for each detector used in the experiment, takes the detector code as input argument. Obtain covariance matrices from the input file.
   */
  void SetupCovarianceMatrices(int);
  
  /// Function to build covariance matrices.
  /*This would be an overriden function of the base function from TCovarianceMinimiation if such class is used in future. */
  void BuildCovarianceMatrices();
  
  ///////////////////////////////////////////////////////////////////////
  //public attributes
  ///////////////////////////////////////////////////////////////////////
  /// TExperiment object as a data member
  TExperiment &fExperiment;
  
  /// Bitset to store what matrices are to be used for minimizatiom/toy generation
  std::map<int,std::bitset<50>> CovMatList;
  
  /// List of all the covariance matrices possible
  // The code WILL NOT function properly if you manually assign a value to any other but the first enum
  enum CovarianceMatrices {
    SSTAT=1, // Signal statisical covariance matrices
    BSTAT, // Background statistical covariance matrix
    SNORM, // Signal normalization matrix
    BNORM, // Background normalization matrix
    SSHAPE, // Signal shape covariance matrix
    BSHAPE, // Background shape covariance matrix
    SESCALE, // Signal energy scale covariance matrix
    SB2B, // Signal bin to bin covariance matrix
  };
  
  /// Number of energy bins, common for all the detectors.
  int fEneBins=0;
  
  /// Number of covariance bins, separate for each detector
  std::map<int,int> fDetCovMatBins;
  
  /// Number of cumulative covariance bins, common for all the detectors.
  int fCumulativeCovMatBins=0;
  
  /// The uncorrelated signal covariance matrix histogram, separate for each detector
  std::map<int,TH2D*> hRCovUCSigMatrix;
  
  /// The fully-correlated signal covariance matrix histogram, same for each detector,
  // EDIT: Need to have only one histogram for each detector
  std::map<int,TH2D*> hRCovFCSigMatrix;
  
  /// The uncorrelated background covariance matrix histogram, separate for each detector
  std::map<int,TH2D*> hRCovUCBkgMatrix;
  
  /// The fully-correlated background covariance matrix histogram, same for each detector
  // EDIT: Need to have only one histogram for each detector
  std::map<int,TH2D*> hRCovFCBkgMatrix;
  
  /// The fully-correlated signal covariance matrix, same for each detector
  // EDIT: Need to have only one matrix for each detector
  std::map<int,TMatrixD> fCovMatrixSigFC;
  
  /// The uncorrelated signal covariance matrix, separate for each detector
  std::map<int,TMatrixD> fCovMatrixSigUC;
  
  
  /// The fully-correlated background covariance matrix, same for each detector
  // EDIT: Need to have only one matrix for each detector
  std::map<int,TMatrixD> fCovMatrixBkgOnFC;
  
  /// The uncorrelated background covariance matrix, separate for each detector
  std::map<int,TMatrixD> fCovMatrixBkgOnUC;
  
  /// The fully-correlated background covariance matrix, same for each detector
  // EDIT: Need to have only one matrix for each detector
  std::map<int,TMatrixD> fCovMatrixBkgOffFC;
  
  /// The uncorrelated background covariance matrix, separate for each detector
  std::map<int,TMatrixD> fCovMatrixBkgOffUC;
  
  /// The correlated signal statical covariance matrix
  std::map<int,TMatrixD> fCovMatrixSigStat;
  
  /// The correlated background (corresponsing to  reactor-on) statical covariance matrix
  std::map<int,TMatrixD> fCovMatrixBkgOnStat;
  
  /// The correlated background (corresponsing to  reactor-off) statical covariance matrix
  std::map<int,TMatrixD> fCovMatrixBkgOffStat;
  
  /// The correlated background statical covariance matrix
  /* Since the background is "measured" when the reactor is off,
  the statistical error on the background should be determined
  by the statistics of the reactor off period.*/
  std::map<int,TMatrixD> fCovMatrixBkgStat;
  
  /// The total covariance matrix, separate for each detector
  std::map<int,TMatrixD> fCovMatrixTot;
  
  /// The total covariance matrix, separate for each detector, this one is not inverted to save to histogram
  std::map<int,TMatrixD> fCovMatrixTotHist;
  
  /// The total Relative covariance matrix, separate for each detector
  std::map<int,TMatrixD> fCovMatrixTotRelative;
  
  /// The total cumulative covariance matrix
  TMatrixD fCumulativeCovMatrix;
  
  /// The total cumulative covariance matrix, this one is not inverted to save to histogram
  TMatrixD fCumulativeCovMatrixHist;
  
  /// The total Relative cumulative covariance matrix
  TMatrixD fCumulativeCovMatrixRelative;
  
  /// Covariance Matrix for reactor on signal
  TMatrixDSym fReactorOnSigCovMatrix;
  
  /// Covariance Matrix for reactor on background
  TMatrixDSym fReactorOnBkgCovMatrix;
  
  /// Covariance Matrix for reactor off background
  TMatrixDSym fReactorOffBkgCovMatrix;
  
  // Final Covariance matrices that are written to the file.
  std::map<int,TH2D*> hSigCovMatrixFC;
  std::map<int,TH2D*> hSigCovMatrixUC;
  std::map<int,TH2D*> hBkgOnCovMatrixFC;
  std::map<int,TH2D*> hBkgOnCovMatrixUC;
  std::map<int,TH2D*> hBkgOffCovMatrixFC;
  std::map<int,TH2D*> hBkgOffCovMatrixUC;
  std::map<int,TH2D*> hSigStatCovMatrix;
  std::map<int,TH2D*> hBkgOnStatCovMatrix;
  std::map<int,TH2D*> hBkgOffStatCovMatrix;
  std::map<int,TH2D*> hBkgStatCovMatrix;
  std::map<int,TH2D*> hTotalCovMatrix;
  std::map<int,TH2D*> hTotalCovMatrixRelative;
  
  TH2D* hCumulativeCovMatrix;
  TH2D* hCumulativeCovMatrixRelative;
  TH2D* hReactorOnSigCovMatrix;
  TH2D* hReactorOnBkgCovMatrix;
  TH2D* hReactorOffBkgCovMatrix;

private:
  ///////////////////////////////////////////////////////////////////////
  //private member functions
  ///////////////////////////////////////////////////////////////////////
  /// Function to setup covariance matrices for a particular detector
  
  void displayMacroError(TString key)
  {
    printf("No macro input with the key %s\n",key.Data());
    printf("Check you macro file\n");
    exit(1);
  }
  
  void displayFileError(TString filename)
  {
    printf("The specified covariance matrix %s file has not been found\n",filename.Data());
    exit(1);
  }
  
  TString GetString(CovarianceMatrices cov)
  {
    switch (cov) {
      case SSTAT: return "SSTAT";
      case BSTAT: return "BSTAT";
      case SNORM: return "SNORM";
      case BNORM: return "BNORM";
      case SSHAPE: return "SSHAPE";
      case BSHAPE: return "BSHAPE";
      case SESCALE: return "SESCALE";
      case SB2B: return "SB2B";
    }
  }
  
  /// Add to the TH2 histogram the covariance matrix obtained from the provided type of covariance matrix
  void ExtractCovarianceMatrix(int,int, TH2D*);
  
  /// Extact the scaling factors for the provided detector number to scale the covariance matrices
  void ExtractScalingVectors(int);
  
  /// Extact the scaling factor to scale the covariance matrices
  void ExtractScalingVector();
  
  /// Setup the covariance matrices for the dtector number provided
  void SetupDetCovarianceMatrices(int);
  
  /// Extract covariance matrices for given detector number
  // Will devolve into something more complicated when covariance matrices need to be defined based on macro files
  void ExtractCovarianceMatrices(int);
  
  /// Function to define and initialize covariance matrix histograms for a given detector number
  void SetupCovarianceMatrixHists(int);
  
  /// Function to fill covariance histograms
  void FillCovHistograms();
  
  /// Function to track the list of covariance matrices to be used for minimization or toy
  void TrackCovMatrixList(int);
  
  
  /// Implement me!!!
  // A function to zero all bins in a covariance matrix based on how the detector is fiducialized
  void FiducializeCovMatrix();
  ///////////////////////////////////////////////////////////////////////
  //private attributes
  ///////////////////////////////////////////////////////////////////////

  /// Defined whether the interface is for toys
  bool isToyInterface=false;
  
  // Detector on and off fractions
  double detOnFraction = 0.0;
  double detOffFraction = 0.0;
  
  /// Individual detector Relative matrix to be used for covariance matrix generation
  /*
   This is the matrix that contains the m_ij=Signal_i*signal_j. This when multiplied element-by element with the total covariance matrix will give the covariance matrix for Relative histograms
   */
  std::map<int,TMatrixD> SignalScalingMatrix;
  /// Cumulative Relative matrix to be used for covariance matrix generation
  TMatrixD TotSignalScalingMatrix;

  /// A temporary histogram to store covariance matrices and store before adding it to the total cov matrix
  TH2D* hTempMat;
};
#endif
