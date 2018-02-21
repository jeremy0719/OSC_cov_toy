///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: June 04, 2016
// Built on top of macros created by K. Gilje
// TODO:
//   Implement error level
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TDETECTOR_HH
#define TDETECTOR_HH

//Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//Include ROOT headers
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"

// OscSens related headers
#include "TDetectorLocation.hh"

/// Detector class, contains all the details pertaining to PROSPECT detector
class TDetector
{
public:
  
  /// number of detectors used, helpful in TExperiment class
  /*
   The member is updted everytime there is a call to a TDetector constructor
   */
  static int DetectorCount;
  
  /// Default constructor
  TDetector(){}
  
  /// constructor that takes detector number and location respectively as the arguments
  TDetector(int, int);
  
  /// Destructor
  virtual ~TDetector(){}
  
  ///////////////////////////////////////////////////////////////////////
  //public attributes
  ///////////////////////////////////////////////////////////////////////
  
  /// An object for the location of the detector
  TDetectorLocation fDetLocation;
  
  /// One integer defining the detector type and position, is defined as 100*detType+detPosition.
  /* The value ranges from 101 to 999, 101 being the near detector in near position.
   */
  int fDetectorCode=0;
  
  /// Detector name
  // Edit: Change this to enum 
  TString fDetectorName;
  
  /// Number of segments in X direction
  int fNSegX;
  
  /// Number of segments in Z direction
  int fNSegZ;
  
  /// Total number of segments
  int fNSeg;
  
  /// Number of fiducial segments
  int fNFiducialSegments;
  
  /// Number of Events to create
  int fNEvents;
  
  /// Position tree
  TTree* fPositionTree;
  
  /// Detector Event Tree
  TTree* fDetectorEventTree;
  
  //--------Things most likely to be common among all (or the near and the far) detectors---------//
  /// Segment width (in meters) along x direction
  /** Unless otherwise specified, segment width in X and Z direction should not be different\n
   *  Here X and Z are defined as the directions along the short length of the segment \n
   *  and Y is defined as the direction along the long length of the segment.
   */
  double fSegWidthX = 0.1461;
  
  /// Segment width (in meters) along z direction
  double fSegWidthZ = 0.1461;
  
  // generic segment width
  /*
   highly unlikely for the segment cross-section to be changed from square to rectangle but just in case, this needs to be changed.
   */
  double fSegWidth = fSegWidthX;
  
  /// Segment length (in meters) in y direction
  double fSegLength = 1.17475;
  
  /// Proton density (m^-3) of the liquid scintillator
  /**TODO: May need to have a separate class for the scintillator properties*/
  double fProtonDensity = 5.46e28; // EJ-309. Detector Protons per volume (m^-3)
  
  /// Efficiency of the detector
  /**TODO: May need to implement segment-based and/or energy-based efficiency*/
  double fDetectorEfficiency = 0.436;
  
  /// Predicted detector position resolution (in meters) in x direction
  double fPosResX = fSegWidthX;
  
  /// Predicted detector position resolution (in meters) in y direction
  double fPosResY = 0.07;
  
  /// Predicted detector position resolution (in meters) in z direction
  double fPosResZ = fSegWidthZ;
  
  /// Predicted energy resolution /sqrt(E)
  double fEnergyResolution = 0.45;
  
  /// Minimum baseline value of the detector
  double fMinBaseline;
  
  /// Maximum baseline of the detector
  double fMaxBaseline;
  
  /// Number of segments along the baseline
  int fNBinsL=0;
  
  /// Number of outer layers, most often, it will be the same in all directions
  int fNOuterXLeft = 1;
  int fNOuterXRight = 1;
  int fNOuterZTop = 1;
  int fNOuterZBottom = 1;
  
  ///////////////////////////////////////////////////////////////////////
  //public member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// Define the detector, called when the detector is instantiated as a default detector
  void DefineDetector(int,int);
  
  /// Move a detector event from Detector Local coordinates to Global Experiment Coordinates.
  void LocalToGlobal();
  
  /// Change a position from Detector Local coordinates to Global Experiment Coordinates.
  void LocalToGlobal(double&, double&, double&, double, double, double);
  
  /// Create Position TTree
  void CreatePositionTree();
  
  /// Create position TTree using the arguments as the
  void CreatePositionTree(int&, int&, int&,double&);
  
  /// This function resets all the flat tree variables to default ones
  void ResetPositionTreeValues();
  
  /// This function fills the tree
  void FillPositionTree();
  
  /// Read position tree into the variables supplied as arguments
  void ReadPositionTree(int& ,int& ,int& ,double&);
  
  /// This function gives access to reading the position tree
  void ReadPositionTree(TString);
  
  /// Write position tree.
  void WritePositionTree(TFile&);
  
  /// Create a tree of events in the detector.
  void CreateDetectorEventTree();
  
  /// This function resets all the tree variables to default in the detector event tree.
  void ResetDetectorEventTree();
  
  /// Create detector events and fill a tree with values.
  void FillDetectorEventTree();

  /// Access to reading the detector event tree, with access to segment ids
  void ReadDetectorEventTree(int&, int&, int&, bool&);
  
  /// Access to reading the detector event tree, with access to global position
  void ReadDetectorEventTree(double&, double&, double&);
  
  /// Access to reading the detector event tree
  void ReadDetectorEventTree(TString);
  
  /// Write the Detector Event Tree
  void WriteDetectorEventTree(TFile&);
  
  /// Print detector info
  void Print();
  
private:
  
  ///////////////////////////////////////////////////////////////////////
  //private attributes
  ///////////////////////////////////////////////////////////////////////
  // EDIT: Remove these variables from data members. Makes more sense to have a separate segment class at a later point.
  /// The index of the segment
  int segmentIndex = 0;
  
  /// The x index of the segment
  int segmentIndexX = 0;
  
  /// The z index of the segment
  int segmentIndexZ = 0;
  
  /// The "measured" baseline to the segment
  double baseline = 0.0;
  
  /// Check whether position tree is created
  bool isPositionTreeCreated = false;
  
  /// The index of the segment of an event
  int eventSegmentIndex = 0;
  
  /// The x index of a segment of an event
  int eventSegmentIndexX = 0;
  
  /// The y index of a segment of an event
  int eventSegmentIndexZ = 0;
  
  /// If an event is in the fiducial volume
  bool isEventFiducial = true;
  
  /// The x position of an event in local detector coordinates
  double eventLocalX = 0.;
  
  /// The y position of an event in local detector coordinates
  double eventLocalY = 0.;
  
  /// The z position of an event in local detector coordinates
  double eventLocalZ = 0.;
  
  /// The x position of an event in global detector coordinates
  double eventGlobalX = 0.;
  
  /// The y position of an event in global detector coordinates
  double eventGlobalY = 0.;
  
  /// The z position of an event in global detector coordinates
  double eventGlobalZ = 0.;
  
  /// Check whether detector event tree is created
  bool isDetectorEventTreeCreated = false;
  
  ///////////////////////////////////////////////////////////////////////
  //private member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// Check whether detector is defined, done by checking if the name of the detector is defined
  void CheckIfDetectorDefined();
  
  /// Compute the paramters of the detector based on the location and type of the detector
  void ComputeDetectorParameters(int,int);
  
  /// Function to extract inputs from macro corresponsing to the Detector Class
  void ReadMacroInputs();
  
};

#endif
