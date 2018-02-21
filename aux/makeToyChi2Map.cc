///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Mar, 2017
// Executable to take a minimization file as input and produce
// 1D chi2 distribution 
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <map>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TList.h"
#include "TKey.h"
#include "PROSPECTStyle.hh"

void usage(int nInputs)
{
  printf("Are you serious!, you think this program could run with %i inputs\n",nInputs);
  printf("Usage: ./makeToyChi2Map inputfile\n");
  exit(0);
}

// Currently the file structure is very rigid, could be change later
void catchFilenameError(TString inputFileName)
{
	printf("The input file name is %s not right\n",inputFileName.Data());
	printf("It should start with Minimization and have sin22 and deltam2 bins in that order\n");
	printf("Example: Minimization_0.100_3.02.root");
	exit(0);
}

// Obtaining the bin numbers based on the name of the input files 
void GetAxesBinNumbers(TFile *inputFile, int &xBin, int &yBin)
{
	TString binValues(inputFile->GetName());
	TObjArray *objArr = binValues.Tokenize("/");
	binValues	= ((TObjString*)objArr->Last())->GetString();
	binValues = binValues.ReplaceAll("Minimization","");
	binValues = binValues.ReplaceAll(".root","");
	objArr = binValues.Tokenize("_");
	if(objArr->GetEntries() != 2)catchFilenameError(inputFile->GetName());
	double s22binValue=(((TObjString*)objArr->First())->GetString()).Atof();
	double dm2binValue=(((TObjString*)objArr->Last())->GetString()).Atof();
	if(!(inputFile->GetListOfKeys()->Contains("CumulativeChiMap0"))){
		printf("CumulativeChiMap0 doesn't exist");
		exit(1);
	}
	TH2D* CumulativeHistogram = (TH2D*)inputFile->Get("CumulativeChiMap0");
	xBin = CumulativeHistogram->GetXaxis()->FindBin(s22binValue); 
	yBin = CumulativeHistogram->GetYaxis()->FindBin(dm2binValue); 
}

int main(int argc, char** argv)
{
  if(argc !=2)
  {
    usage(argc);
    return -1;
  }
  
  TString inputFileName(argv[1]);
  TFile *inputFile = TFile::Open(inputFileName);
  TString outputFileName(argv[1]);
	if(!(outputFileName.Contains("Minimization"))) catchFilenameError(outputFileName);
  outputFileName.ReplaceAll("Minimization","Chi2Values");
  TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");

	int actualXBin=-1;
	int actualYBin =-1;
	GetAxesBinNumbers(inputFile,actualXBin,actualYBin);
	printf("ΔΧ2 values are calculated w.r.t the original points corresponding to the bins for sin2θ=%i and Δm2=%i",actualXBin,actualYBin);
	
	//Load style
  SetupProspectStyle();
  
  TH1D *chi2Values = new TH1D("Chi2Values","Chi2Values",2000,0,5000);
  
  TList *listOfKeys=inputFile->GetListOfKeys();
  for (int i = 0; i < listOfKeys->GetSize(); i++) {
    TString keyName=((TKey*)listOfKeys->At(i))->GetName();
    if(!(keyName.Contains("CumulativeChiMap"))) continue;
    TH2D* CumulativeHistogram;
    CumulativeHistogram = (TH2D*)inputFile->Get(keyName.Data());
    int x,y,z;
    CumulativeHistogram->GetMinimumBin(x,y,z);
    double bestBinCont=CumulativeHistogram->GetBinContent(x,y);
		double actualBinCont=CumulativeHistogram->GetBinContent(actualXBin,actualYBin);
    chi2Values->Fill(actualBinCont-bestBinCont);
		std::cout << actualBinCont << "   " << actualXBin << " "<< actualYBin<<std::endl;
		std::cout << bestBinCont << "   " << x << " " << y<<std::endl;
		std::cout <<" " <<std::endl;
	}
  outputFile->cd();
  chi2Values->Write();
  outputFile->Close();
}
