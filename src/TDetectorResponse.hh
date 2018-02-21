///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Aug 2017
// TODO:
//   Implement error level
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TDETECTORRESPONSE_HH
#define TDETECTORRESPONSE_HH

//Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cassert>

//Include ROOT headers
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TVector3.h"

///Class for Detector reponse
/**
 Apply detector response to the provided histograms.
 The smearing is
 */
class TDetectorResponse
{
public:
  
  /// Default constructor
  TDetectorResponse(){}
  
  /// Default destructor
  ~TDetectorResponse(){}
  
  ///////////////////////////////////////////////////////////////////////
  //public attributes
  ///////////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////////
  //public member functions
  ///////////////////////////////////////////////////////////////////////
  
  
  /*
   Apply gaussian smearing of a given value to the input histogram.
   The following energy-based smearing is used σ=smearing*sqrt(E)
   Resolution smearing from Xianyi Zhang, modified by Pranava
   Assumes that the binning is regular. Doesn't work for variable binning
   Works only for the histograms with minimum bin value is > 0
   */
  void ApplyResponse(const TH1 &,TH1*,double);
  
  /**
   Apply gaussian smearing of a given value to the input histogram.
   This function extracts 1D histograms corresponding to each Y bin and smears that particular histogram by calling on void ApplyResponse"("TH1D*,double")"
   */
  void ApplyResponse(const TH2 &,TH2*, double);
  
  
  /// Apply response using a detector response matrix
  void ApplyResponse(const TH1 &, TH1*,TH2*);
  
  /** Apply response using a detector response matrix
   This function extracts 1D histograms corresponding to each Y bin and applies detector response to that particular histogram by calling on void ApplyResponse"("TH1&, TH1*,TH2*")"
   */
  void ApplyResponse(const TH2 &, TH2*,TH2*);
  
  //Apply response to a set of histograms using detector response matrices
  /**
   Apply response to a set of histograms using detector response matrices that are input as a map of histograms. Each key in the map corresponds to the hisotgram number i.e., if you are applying response to a bunch of PROSPECT segments, the key corresponds to the segment number
   */
  void ApplyResponse(const TH2 &, TH2*,const std::map<int,TH2*> &);
  
  // by xiaobin lu 2018.02.15 modification
  ///////////////////////////////////////////////////////////////////////
  // overloading each segmented detector response function
  // int i = 1 for displaying overloading info
  // void ApplyResponse(const TH2 &, TH2* , TH2* , int  );

  

















  ///////////////////////////////////////////////////////////////////////

protected:
  
  ///////////////////////////////////////////////////////////////////////
  //private attributes
  ///////////////////////////////////////////////////////////////////////
private:
  
  ///////////////////////////////////////////////////////////////////////
  //private attributes
  ///////////////////////////////////////////////////////////////////////
};
#endif
