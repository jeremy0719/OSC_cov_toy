///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Sep 2017
// TODO:
//   Implement error level
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TDETECTORRESPONSEINT_HH
#define TDETECTORRESPONSEINT_HH

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

// OscSens specific headers
#include "TDetectorResponse.hh"
#include "TMacroInterface.hh"

///Class for interfacing between TDetectorResponseInterface and the rest of the classes.
/*
 This is the class where all the inputs histograms for detector reponse is extracted from
 and given as input to the TDetectorReponse class.
 */
class TDetectorResponseInterface
{
public:
  
  /// Default constructor
  TDetectorResponseInterface();
  
  /// Constructor that takes the detectorcode as input
  TDetectorResponseInterface(int);
  
  /// Default destructor
  ~TDetectorResponseInterface();
  
  ///////////////////////////////////////////////////////////////////////
  //public attributes
  ///////////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////////
  //public member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// Apply response using a detector response matrix
  /*
   The detector response matrix is extracted before inputting to ApplyResponse function in TDetectorResponse.
   */
  void ApplyDetectorResponse(const TH1&,TH1*);
  
  /// Apply response using a detector response matrix.
  /**This function extracts 1D histograms corresponding to each Y bin and applies detector response to that particular histogram by calling on void ApplyResponse"("TH1*,std::map<int,TH2*>")" of the #detResp
   */
  void ApplyDetectorResponse(const TH2&, TH2*);
  
  /*
   Apply gaussian smearing (smearing/sart(E)) of a given value to the input histogram.
   The following energy-based smearing is used σ=smearing*sqrt(E)
   */
  void ApplyDetectorResponse(const TH1 &inputHist,TH1* outputHist, double smearing){detResp.ApplyResponse(inputHist,outputHist,smearing);}
  
  /**
   Apply gaussian smearing of a given value to the input histogram.
   This function extracts 1D histograms corresponding to each Y bin and smears that particular histogram by calling on void ApplyDetectorResponse"("TH1D*,double")"
   */
  
  void ApplyDetectorResponse(const TH2 &inputHist,TH2* outputHist, double smearing){detResp.ApplyResponse(inputHist,outputHist,smearing);}

  /*
  modification by xiaobin lu 2018.02.20 
  */
  void ApplyDetectorResponse(const TH2 &inputHist, TH2* outputHist, TH1* outputHistFull);



  
protected:
  
  ///////////////////////////////////////////////////////////////////////
  //private attributes
  ///////////////////////////////////////////////////////////////////////
private:

  ///////////////////////////////////////////////////////////////////////
  //private attributes
  ///////////////////////////////////////////////////////////////////////
  
  TDetectorResponse detResp;
  
  TString detRespFileName;
  
  /// Code for the the detector for extracting detector response matrix.
  /* By default the detector at near position is taken.
   */
  int detectorCode = 101;
  
  /// Boolean to check if the constructor called to create the object is default or not.
  bool isDefaultConstructor=true;
  
  /// Get energy resolution for the detector.
  double fEnergyResolution=0.045;
  
  /// Detector response matrix file that contains detector response for each individual file and the full detector response matrix.
  TFile *detRespFile;
};
#endif
