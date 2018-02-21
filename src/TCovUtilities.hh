///////////////////////////////////////////////////////////////////////////
///// Name : CovUtilities.hh
///// Desc : Utility functions for Covariance Analysis
///// Author : P T Surukuchi
///// Last edited : 2016
///// TODO :
/////
//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifndef TCOVUTILITIES_HH
#define TCOVUTILITIES_HH

// ROOT headers
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
// cpp headers
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>

TH2D* rebinTH2(TH2*);

std::vector<TH1D*> CreateThrowHistograms( int, double);

std::vector<TH1D*> CreateToyHistograms(std::vector<TH1D*>);

#endif
