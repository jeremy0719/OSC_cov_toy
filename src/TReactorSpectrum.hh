///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: June 23, 2016
// TODO:
//   Implement error level
//   In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TREACTORSPECTRUM_HH
#define TREACTORSPECTRUM_HH

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

//Include Osc_Sens related headers
#include "TReactor.hh"

/// Dummy class for working with spetrum from the reactor, need to be implemented for using spectra from different sources (e.g HM, DL etc )
class TReactorSpectrum
{
public:
  TReactorSpectrum();
};

#endif