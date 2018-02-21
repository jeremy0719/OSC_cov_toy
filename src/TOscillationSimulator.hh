///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Oct 2016
// Built on top of macros created by K. Gilje
// TODO: Implement error level
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TOSCILLATIONSIMULATOR_HH
#define TOSCILLATIONSIMULATOR_HH

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
#include "TMath.h"

/// Oscillation Class
/**
 * Utility class for obtaining oscillated plots from unoscillated plots
 */
class TOscillationSimulator
{
public:
  
  /// Default Constructor
  TOscillationSimulator();
  
  /// Default destructor
  ~TOscillationSimulator();
  
  ///////////////////////////////////////////////////////////////////////
  //public member functions
  ///////////////////////////////////////////////////////////////////////
  
  /// Obtain the m2 term in the antineutrino disappearance/appearance probability
  /**
   Return value of the function is \f$ \sin^{2}(\frac{1.27 \delta m^{2} L }{E})\f$
   */
  double GetDeltam2Term(double, double, double);
  
  /// Compute \f$ \bar\nu_e \f$ disappearance probability or \f$ P_{e \rightarrow \alpha} \f$
  /**
   Return value of the function is \f$ \sin^{2}2\theta \sin^{2}(\frac{1.27 \delta m^{2} L }{E})\f$
   */
  double ComputeNueDisProb(double, double , double, double);
  
  /// Compute  \f$ \bar\nu_e \f$ survival probability \f$ P_{e \rightarrow e}\f$
  /**
   Return value of the function is \f$ 1- \sin^{2}2\theta \sin^{2}(\frac{1.27 \delta m^{2} L }{E})\f$
   */
  double ComputeNueSurvProb(double , double, double, double);
  
  /// Oscillate an LvsE histogram based on the deltam2 and Sin2theta values
  TH2D* Oscillate2DHistogram(const TH2D*,double, double);
  
private:
};

#endif
