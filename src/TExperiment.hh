///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Oct 2016
// Built on top of macros created by K. Gilje
// TODO:
//   Implement error level
//   Implement ability to create an object of this class with an input data file.
//   Implement the ability to choose to save histograms only when needed
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TEXPERIMENT_HH
#define TEXPERIMENT_HH

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
#include "TH3.h"
#include "TFile.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TF1.h"

//Include Osc_Sens related headers
#include "TReactor.hh"
#include "TDetector.hh"
#include "TMacroInterface.hh"
#include "TOscillationSimulator.hh"
#include "TDetectorResponseInterface.hh"
#include "TOscSensUtils.hh"

///Class for experiment
/**
 * Information regarding the experiment. Includes information pertaining to the Detector, reactor, spectrum from the reactor and everything else (cross-section, branching ration etc ) related to the experiment.
 */
class TExperiment
{
public:
  
  /// Default constructor
  TExperiment();
  
  /// Default constructor
  TExperiment(double, double);
  
  /// Default destructor
  ~TExperiment();
  
  ///////////////////////////////////////////////////////////////////////
  //public attributes
  ///////////////////////////////////////////////////////////////////////
  
  /// Detector objects as members
  /*
   Since each experiment might have more than one detector, a vector of detector objects were used
   */
  std::vector<TDetector*> fDetectors;
  
  /// Event tree
  /**
   * Tree storing information of neutrino origin, detection, baselines, weights etc
   * In future, might need to evolve into something more complicated.
   */
  std::map<int,TTree*> fEventTree;
  
  /// Reactor object, all the reactor-related information can be obtained from here
  /*
   May need to implement a vector if an experiment will have more than one reactor, I do not foresee any need for that now.
   */
  TReactor fReactor;
  
  /// The number of generated fake events
  int fNEvents = 1000000;
  
  /// The span of time (in number of seconds) used for detector data taking during reactor-on for each detector
  std::map<int, double> fRxOnExposure;
  
  /// The span of time (in number of seconds) used for detector data taking during reactor-off for each detector
  std::map<int, double> fRxOffExposure;
  
  std::map<int,double> fRxOnOffRatio;
  
  /// Conversion of Years to Seconds
  // EDIT: Change to static const
  const double fYearToSeconds = 3.15569e7;
  const double fDayToSeconds = 86400;
  // Warning:fCycleToSeconds might be different for different cycles
  const double fCycleToSeconds = 2134080;
  
  /// Oscillation parameters for the reference model
  double fReferenceDeltam2=0.0;
  double fReferenceSin22theta=0.0;
  
  // Oscillation simulator object
  TOscillationSimulator fOscillationSimulator;
  
  /// Simulated antineutrino spectrum from the reactor.
  /*
   Might need a separate class TAntiNueSpectrum
   */
  TH1D *hAntiNuSpectrum;
  TH1D *hIBDCrossSection;
  TH1D *hObsSpectrum;
  
  /// Maps the measured baseline bins to the segment
  std::map<int,TH2D*> hPositionToSegment;
  
  /// Maps the measured baseline bins to the fiducial segment
  /**
   * This cuts out the position bins that refer to non-fiducial segments
   */
  std::map<int,TH2D*> hPositionToSegmentFid;
  
  /// Maps the measured baseline bins to the segment
  std::map<int,TH2D*> hBaselineToSegment;
  
  /// Pseudo detector, represents all the simulated events in the detector volume
  std::map<int,TH3D*> hPseudoDetector;

  /// Simulated neutrino energy
  /* This is different from the #hAntiNuSpectrum in that it is properly binned and
   each detector gets it's own spectrum
   */
  std::map<int,TH1D*> hSimulatedNueEnergy;
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  // Beginning of all hists relevant to minimization
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////
  // All hists needed for the MC side of the spectrum
  ////////////////////////////////////////////////////////////////////////////////
  
  /// Unoscillated simulated absolute neutrino energy spectrum
  std::map<int,TH1D*> hETrue;
  
  /// Unoscillated simulated absolute deposited energy spectrum
  std::map<int,TH1D*> hENull;
  
  /// Unoscillated simulated neutrino energy spectrum as a function of segment
  std::map<int,TH2D*> hSegvsE;

  /// Unoscillated simulated neutrino energy spectrum as a function of segment
  std::map<int,TH2D*> hSegvsETrue;
  
  /// Unoscillated simulated deposited energy spectrum as a function of segment
  std::map<int,TH2D*> hSegvsENull;
  
  /// Unoscillated simulated deposited energy spectrum as a function of baseline
  std::map<int,TH2D*> hLvsENull;
  
  /// Unoscillated LvsE deposited spectrum that has been subtracted from the scaled absolute deposited energy spectrum #hENull
  std::map<int,TH2D*> hLvsERelative;
  
  ////////////////////////////////////////////////////////////////////////////////
  // All hists needed for the reference side of the analysis.
  // In case of sensitivity estimate, these are the same hists as the null versions,
  // in case of oscillation-fiiting, these are the data files and
  // in case of fake oscillation-fitting (i.e, fitting to MC oscillated hists),
  // these are the oscillated versions of the null hists.
  ////////////////////////////////////////////////////////////////////////////////
  
  /// Oscillated absolute neutrino spectrum.
  /*
   In case of sensitivity estimate, this is same as #hETrue,
   in case of oscillation-fitting, this irrelavant and hence not used
   in case of fake oscillation-fitting (i.e, fitting to MC oscillated hists),
   this is the oscillated version of #hETrue with the oscillation parameters #fReferenceDeltam2 and #fReferenceSin22theta.
   */
  std::map<int,TH1D*> hEOsc;
  
  /// Absolute reconstructed spectrum
  /*
   In case of sensitivity estimate, this is same as #hENull,
   in case of oscillation-fitting, this is background subtracted detected IBD spectrum,
   in case of fake oscillation-fitting (i.e, fitting to MC oscillated hists),
   this is the oscillated version of #hENull with the oscillation parameters #fReferenceDeltam2 and #fReferenceSin22theta.
   */
  std::map<int,TH1D*> hERef;

  /// Absolute oscillated neutrino spectrum as a fucntion of segment
  /*
   In case of sensitivity estimate, this is same as #hSegvsETrue,
   in case of oscillation-fitting, this is not relevant and so not used,
   in case of fake oscillation-fitting (i.e, fitting to MC oscillated hists),
   this is the oscillated version of #hSegvsENull with the oscillation parameters #fReferenceDeltam2 and #fReferenceSin22theta.
   */
  std::map<int,TH2D*> hSegvsEOsc;
  
  /// Absolute reconstructed spectrum as a fucntuon of segment
  /*
   In case of sensitivity estimate, this is same as #hSegvsENull,
   in case of oscillation-fitting, this is background subtracted detected IBD spectrum in each segment,
   in case of fake oscillation-fitting (i.e, fitting to MC oscillated hists),
   this is the oscillated version of #hSegvsENull with the oscillation parameters #fReferenceDeltam2 and #fReferenceSin22theta.
   */
  std::map<int,TH2D*> hSegvsERef;
  
  /// Absolute reconstructed spectrum as a fucntion of segment
  /*
   In case of sensitivity estimate, this is same as #hLvsENull,
   in case of oscillation-fitting, this is background subtracted detected IBD spectrum in each segment,
   in case of fake oscillation-fitting (i.e, fitting to MC oscillated hists),
   this is the oscillated version of #hLvsENull with the oscillation parameters #fReferenceDeltam2 and #fReferenceSin22theta.
   */
  std::map<int,TH2D*> hLvsERef;
  

  /// LvsE reconstructed spectrum that has been subtracted from the Relative absolute reconstructed energy spectrum #hERef
  std::map<int,TH2D*> hLvsERelativeRef;
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  // End of all hists relevant to minimization
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  /// Reactor on background 1D histogram
  /*
   For now identical for all the segments, might need to use a vector of histograms one for each segment.
   */
  std::map<int,TH1D*> hBkgRxon1D;
  
  /// Reactor on backgrounds as a function of segment
  std::map<int,TH2D*> hBkgRxonSeg;
  
  /// Reactor on backgrounds as a function of baseline
  std::map<int,TH2D*> hBkgRxOnLvsE;
  
  /// Reactor off backgrounds as a function of baseline
  std::map<int,TH2D*> hBkgRxOffLvsE;
  
  /// Segment to segment signal to background histogram, currently populated only if the reference is data
  TH2D* hStoB;
  
  /// Segment to segment reactor-on total rates, currently populated only if the reference is data
  TH2D* hRxOn;
  
  /// Segment to segment reactor-off total rates, currently populated only if the reference is data
  TH2D* hRxOff;
  
  ///////////////////////////////////////////////////////////////////////
  //public member functions
  ///////////////////////////////////////////////////////////////////////
  
  // Method to Append a detector to the experiment. The inputs are detector type (near, far etc) and detector location and number of cycles respectively. The default is 7 cycles.
  void AddDetector(int, int, double nCycles=7.0);
  
  /// Method to read energy spectrum information and fill the energy, flux and cross-section vectors
  void ReadEnergySpectrum(TString energySpectrumFileName="./inputs/EnergyTableHeu.txt");
  
  /// Create trees
  void CreateTrees();
  
  /// Filling both event and position trees
  void FillTrees();
  
  /// Write all the trees related to the experiment
  void WriteTrees(TFile&);
  
  /// Setup experiment
  void SetupExperiment();
  
  /// Fill and write all the histograms
  void WriteHistograms(TFile&);
  
  /// Print information pertaining to the experiment
  void Print();
  
  /// Function to Populate the details in TExperiment
  void PopulateExperimentFromData();
  
protected:
  
  ///////////////////////////////////////////////////////////////////////
  //protected attributes
  ///////////////////////////////////////////////////////////////////////
  
  /// This tells us if the reference histogram comes from the data (TDataExtractor) if it is a model
  bool isReferenceData=false;
  
  /// Raw reactor coordinates
  TVector3* reactorPos = new TVector3;
  
  /// Raw detector coordinates
  TVector3* detectorPos = new TVector3;
  
  /// The segment id of the event
  int segmentId = -1;
  
  /// The X segment id of the event
  int segmentIdX = -1;
  
  /// The Z segment id of the event
  int segmentIdZ = -1;
  
  /// Flag if event is in the fiducial volume
  bool isFiducial = true;
  
  /// Distance from Reactor-Pixel to Detector-Pixel
  /// Also known as the true baseline for the oscillation calculation.
  double trueBaseline = 0.0;
  
  /// Weight of event from baseline
  double weight = 0.0;
  
private:
  
  ///////////////////////////////////////////////////////////////////////
  //private attributes
  ///////////////////////////////////////////////////////////////////////
  
  // Reconstructed energy bins
  const int nTrueEBins = 67;
  const int nRecoEBins = 32;
  const double trueEBinWidth=0.1;
  const double recoEBinWidth=0.2;
  const double trueMinEEnergy = 1.8;
  const double recoMinEEnergy = 0.8;
  
  // Number of fiducial bins, key is the detnumber and value is the no of fidbins
  std::map<int,int> nBinsFid;
  
  /// Flag to check whether the experiment is defined or not
  bool isExperimentDefined=false;
  
  /// Flag to check whether trees are created
  bool areTreesCreated = false;
  
  /// Vector with antinu energy
  std::vector<double> fAntiNuEnergy;
  
  /// Vector with the energy dependent incoming antineutrino flux
  std::vector<double> fAntiNuFlux;
  
  /// Vector with the energy dependent IBD cross section
  std::vector<double> fIBDCrossSection;
  
  /// Data location
  TString dataFileName="ExtractedData.root";
  
  ///////////////////////////////////////////////////////////////////////
  //private methods
  ///////////////////////////////////////////////////////////////////////
  // EventTree related functions
  
  /// Returns the number of fiducial bins
  int GetNFidBins(const TDetector &);
  
  /// Function to generate fiducial binning using the non-fiducialised histogram for a specific detector.
  void GenerateFidBinning(const TDetector &,double*,int);
  
  /// Function to generate segment binning for a specific detector.
  void GenerateSegBinning(const TDetector &,double*,int);
  
  /// Function to generate true energy binning
  void GenerateTrueEBinning(double*,int);
  
  /// Function to generate reconstructed energy binning
  void GenerateRecoEBinning(double*,int);
  
  /// Method to check if the experiment is defined and if it is not exit
  void CheckIfExperimentDefined();
  
  /// Create the event tree
  void CreateEventTree();
  
  /// This function resets all the flat tree variables to default ones
  void ResetEventTreeValues();
  
  /// This function fills the tree
  void FillEventTree();
  
  /// Fills the event tree
  /**
   * The reactor positions of the event tree is created from an input file (either .root or textfile) whose path is provided by the TString argument
   */
  // TODO: Implement this, this will end up having to access TReactorSpectrum class
  void FillEventTree(TString);
  
  /// Write energy spectrum
  void WriteEnergySpectrum(TFile&);
  
  /// Write the event tree
  void WriteEventTree(TFile&);
  
  /// Setup energy spectrum histograms
  void SetupEnergySpectrumHists();
  
  /// Setup all the histograms
  void SetupSignalHists();
  
  /// Setup background histograms
  void SetupBackgroundHistograms(TString);
  
  /// Function to extract inputs from macro corresponsing to the Experiment Class
  void ReadMacroInputs(int);
  
};
#endif
