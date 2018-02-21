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
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TKey.h"
#include "PROSPECTStyle.hh"

void usage(int nInputs)
{
  printf("Are you serious!, you think this program could run with %i inputs\n",nInputs);
  printf("Usage: ./generateCriticalChi2Plot input_directory Reference_CumulativeChi2Plot Confidence Level\n");
  exit(0);
}

void catchFilenameError(TString inputFileName)
{
	printf("The input file name is %s not right\n",inputFileName.Data());
	printf("It should start with Chi2Values and have sin22 and deltam2 bins seperated by '_' in that order\n");
	printf("Example: Minimization_0.100_3.02.root\n");
	exit(0);
}

void GetAxesValues(TFile *inputFile, double &xValue, double &yValue)
{
	TString binValues(inputFile->GetName());
	TObjArray *objArr = binValues.Tokenize("/");
	binValues	= ((TObjString*)objArr->Last())->GetString();
	binValues = binValues.ReplaceAll("Chi2Values","");
	binValues = binValues.ReplaceAll(".root","");
	objArr = binValues.Tokenize("_");
	if(objArr->GetEntries() != 2) catchFilenameError(inputFile->GetName());
	xValue=(((TObjString*)objArr->First())->GetString()).Atof();
	yValue=(((TObjString*)objArr->Last())->GetString()).Atof();
}

int main(int argc, char** argv)
{
  if(argc !=4)
  {
    usage(argc);
    return -1;
  }
	
		
  TString inputFileLocation(argv[1]);
	TString inputReferenceFile(argv[2]);
	TFile *inputFile=TFile::Open(inputReferenceFile.Data());
	if(!(inputFile->IsOpen()) || inputFile->IsZombie()){
		printf("Check the reference file supplied");
		exit(1);
	}
	
	if(!(inputFile->GetListOfKeys()->Contains("CumulativeChiMap0"))){
		printf("CumulativeChiMap0 doesn't exist\n");
		exit(1);
	}
	TH2D *hRef=(TH2D*)inputFile->Get("CumulativeChiMap0");

  double confidenceLevel=std::atof(argv[3]);
	TString outputFileName(inputFileLocation);
	outputFileName.Append("/CriticalChi2Map");
	outputFileName.Append(argv[3]);
	outputFileName.Append(".root");
	std::cout << outputFileName.Data() <<std::endl;
	TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
	TSystemDirectory *inputDirectory=new TSystemDirectory("Chi2ValueDirectory",inputFileLocation.Data());
  TList *fileList =inputDirectory->GetListOfFiles();
	TIter it(fileList);
	
	TH2D* hCriticalChi2Map;
	hCriticalChi2Map=(TH2D*)hRef->Clone("hCriticalChi2Map");
	hCriticalChi2Map->Reset();

	double xq[1];
	double yq[1];
	xq[0] = confidenceLevel;
	for(int i=0;i<fileList->GetSize();i++){
		TSystemFile *file = (TSystemFile *)fileList->At(i);
		if(file->IsDirectory()) continue;
		TString fileName=file->GetName();
		if(!(fileName.Contains(".root"))) continue;
		if(!(fileName.Contains("Chi2Values_"))) continue;

		fileName.Prepend(inputFileLocation.Data());
		TFile *inputFile = TFile::Open(fileName);
		TH1D* chi2Map =(TH1D*)inputFile->Get("Chi2Values");
		double xValue, yValue;
		GetAxesValues(inputFile,xValue,yValue);
		int xBin=(hRef->GetXaxis())->FindBin(xValue);
		int yBin=(hRef->GetYaxis())->FindBin(yValue);
		chi2Map->GetQuantiles(1,yq,xq);
		hCriticalChi2Map->SetBinContent(xBin,yBin,yq[0]);
		inputFile->Close();
	}
  outputFile->cd();
  hCriticalChi2Map->Write();
  outputFile->Close();
}
