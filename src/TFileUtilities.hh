///////////////////////////////////////////////////////////////////////////
///// Name : TFileUtilities.hh
///// Desc : Utlilities pertaining to dealing with other format of files
///// Author : P T Surukuchi
///// Last edited : 2017/10/16
//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifndef TFILEUTILITIES_HH
#define TFILEUTILITIES_HH

#include <cassert>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
//#include "RioStream.h"

#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>


void usage();

// Extract Data from a 1D histogram with a given name, takes name of the ROOT file and name of the histogram in the ROOT file in that order
void hDataExtractor(TString, TString);

// Extract Data from a 2D histogram with a given name, takes name of the ROOT file and name of the histogram in the ROOT file in that order
void tdhDataExtractor(TString,TString);

TH1D* makeTXTTH1(TString,int,int);

TH2D* makeTXTTH2(TString,int,int);

void makeTXTTGraph(TString,TGraph *);

#endif
