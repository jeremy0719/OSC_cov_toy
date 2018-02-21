///////////////////////////////////////////////////////////////////////////
///// Name : TOscSensUtils.hh
///// Desc : Utility functions for Oscillation Sensitivity framework
///// Author : P T Surukuchi
///// Last edited : 2017
///// TODO :
/////
//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifndef TOSCSENSUTILITIES_HH
#define TOSCSENSUTILITIES_HH

// ROOT headers
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include "TVectorD.h"
// cpp headers
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>

//OscSens specific includes
#include "TExperiment.hh"

/// Obtain number of bins in a 2D histogram
int Get2DHistBins(const TH2D*);

/// A function to take an input map of TH2D histograms and turns them into TVectorD
void ExtractVector(TH2D*, TVectorD&);

/// A function to take a TH2D histograms and turns them into TVectorD
void ExtractVector(std::map<int, TH2D*>,TVectorD&);

/// A function to take an input map of TH2D histograms and turns them into TVectorD
void ExtractInvertedVector(TH2D*, TVectorD&);

/// A function to take a TH2D histograms and turns them into TVectorD
void ExtractInvertedVector(std::map<int, TH2D*>,TVectorD&);

/// Relativize the histogram

/* Takes an input LvsE histogram and spits out a histgram that has been rescaled to the simulated LvsE histogram.
 If inputLvsE and inputE are input LvsE(L,E) and E(E) histograms respectively where LvsE is a funtion of baseline and energy and E is function of energy. The output LvsE histogram is given by LvsERelative(L,E):
 \f[
 LvsERelative(L,E)=\frac{inputLvsE(L,E)/inputE(E)}{LvsENull(L,E)/ENull(E)}
 \f]
 */
void RelativizeHistogram(const TH2D&,const TH1D&,const TH2D&, const TH1D&, TH2D&);

/// Check for consistecy of 1D histograms
/*
 Checks if the axes ranges and bin centers are same for both histograms
 */
bool CheckConsistency(const TH1 &,const TH1 &);

/// Check for consistecy of 2D histograms
/*
 Checks if the axes ranges and bin centers are same for both histograms
 */
bool CheckConsistency(const TH2 &,const TH2 &);
#endif
