///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: June 05, 2016
// Built on top of macros created by K. Gilje
// TODO:
//   Implement error level
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TDETECTORLOCATION_HH
#define TDETECTORLOCATION_HH

//Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//Include ROOT headers
#include "TString.h"
#include "TMath.h"

//Include Osc_Sens related headers
#include "TMacroInterface.hh"

/// Struct for the detector location of PROSPECT.
/**
 * Should be used by the detectors both near and far to identify the location and the distances of the segments, might need a complex implementation later on in future
 */
class TDetectorLocation
{
public:
  /// Default constructor
  TDetectorLocation();
  
  /// Constructor, takes detector type and location number as the input
  TDetectorLocation(int, int);
  
  /// Destructor
  virtual ~TDetectorLocation(){}
  
  ///////////////////////////////////////////////////////////////////////
  //public attributes
  ///////////////////////////////////////////////////////////////////////
  
  /// Identifier for detector type, 1 for PROSPECT Near and 2 for PROSPECT Far
  int fDetectorType;
  
  /// Identifier for the detector position
  int fDetectorPosition;
  
  /// Detector X distance from reactor
  /*
   Distance from the reactor center to the front of the detector in the X direction
   Note: Here X is defined at the radial distance, Y is defined as the distance along the length of the segments, and Z is defined in vertical direction
   */
  double fDistX; // meters
  
  /// Detector Y distance from reactor
  double fDistY; // meters
  
  /// Detector Z distance from reactor
  double fDistZ; // meters
  
  /// Detector angular rotation with respect to radial axis of reactor.
  double fAngleOffset=0.0; // degrees
  
  ///////////////////////////////////////////////////////////////////////
  //public member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// Assign the detector type and position in if the object was initialized by default constuctor
  void DefineDetectorLocation(int, int);
   
  /// Print detector info
  void Print();
  
private:
  
  
  ///////////////////////////////////////////////////////////////////////
  //private member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// Method to compute distances of the detector based on the type and position
  void ComputeDistances();
  
  /// Function to read inputs from macro file
  void ReadMacroInputs();
};

#endif
