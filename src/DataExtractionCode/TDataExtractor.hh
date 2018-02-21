///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: June 23, 2016
// TODO:
// Take 2 input files, one for signal another for background
//   Implement error level
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TDATAEXTRACTOR_HH
#define TDATAEXTRACTOR_HH

//Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//Include ROOT headers
#include "TChain.h"
#include "TSystemDirectory.h"
#include "TH2D.h"
#include "TCollection.h"
#include "TFile.h"

// OscSens specific includes
#include "TDetector.hh"
#include "TExperiment.hh"

///Class for extracting PROSPECT data and converting into usable format.
class TDataExtractor
{
public:
  
  /// Constructor where the detector baselines are taken from TDetector
  TDataExtractor(TString extName,TExperiment &,TDetector &);
  
  /// Function to add a ttree, the argument has to be the parent directory in which the root files and subdirectories sit
  void AddTTree(TString);
 
  /// TChain for IBD ttrees
  TChain *fInputDataTChain;
  
  /// TTree for IBD-like events
  TTree *fInputDataTree;
  
  TBranch *bLineBranch;
  
  /// Segment vs E histograms
  std::map<int,TH1D*> hSegE;
  
  /// Segment vs E histograms
  TH2D *hSegvsE;

  
  /// Baseline vs E histograms
  TH2D *hLvsE;
  
  /// Baseline vs E histograms
  TH1D *hE;
  
  /// Time of the run, i.e., time differnce between the first and the last event in the run
  double fT;
  
  /* 
   Function to extract data from the root file given as the input in the AddTTree 
   method and fill the histograms LvsE, SegvsE and other relavant histograms
   */
  void ExtractData(); 
  
  void WriteHistograms(TFile &);
  /// Power-weighted time, normalized to 85 MW
  //std::map<int,double> fTPower;
  
private:
  ///////////////////////////////////////////////////////////////////////
  //private member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// Default constructor
  //TDataExtractor(){}
  
  /// Create histograms used in TExperiment
  void CreateHistograms();
  
  /// Fianlize the addition of the TTrees to the TChain, past this point there should not be any addition of TTrees
  void FinalizeTTreeAddition();
  
  /// Read data tree
  void ReadDataTree();
  
  ///////////////////////////////////////////////////////////////////////
  //private attributes
  ///////////////////////////////////////////////////////////////////////
  
  bool isFinal=false;
  
  /// Number of tress used to make the tchain
  int nFiles=0;
  
  /// TExperiment object this class uses for some information
  TExperiment &fExperiment;
  
  /// TDetector object that the dataextractor is used for.
  TDetector &fDetector;
  
  /// Event number
  Long64_t evt = -1;
  
  /// Absolute time of the event (seconds)
  double t_abs=-1;
  
  /// Energy of the event (MeV)
  float E = -1;
  
  /// Max segment
  int maxseg = -1;
  
  /// array of events location in X, Y and Z
  float xyz[3];
  
  /// Energy of the maximum segment (MeV)
  float E_maxseg=-1;
  
  /// Energy of the adjacent segment (MeV)
  float E_adjacent=-1;
  
  /// Segment multiplicity
  int segmult=-1;
  
  /// dimater of the event
  float diameter=-1;
  
  /// Time spread (between the PMTs ?) of the the event (ns?)
  float tspread = -1;
  
  /// Neutron capture dt(μs)
  float ncapt_dt = -1;
  
  /// Segment number of neutron
  int n_seg=-1;
  
  /// array of floats
  float n_xyz[3];
  
  /// time to nearest veto event [ns]
  float veto_t=-1;
  
  /// Cut number of the events
  int cut;
  
  /// Geometry flags corresponing to the position of the detector
  int detgeom=-1;
  
  /// Reactor power
  float rxpwr = -1 ;
  
  //Baseline lengths, loaded from input in a macro file

  double baselineX;
  double baselineZ;
 
  /// Reconstructed baseline of the event, will be used by the TExperiment class
  double bLine;
  
  TString extractorName;
  
  // Detector code being used to generate the baselines
  int fDetCode;
  
  bool isDataExtracted=false;
  
  bool isTDetectorDefined=false;
  
};
#endif
