/// A file to run simple tests on reliability of the fitter.
/// Author: P. T. Surukuchi
/// Dec 2017

#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <map>

#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TCollection.h"
#include "TList.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TH2D.h"
#include "TMath.h"

void usage(int arguments)
{
  printf("Cannot use %i arguments\n",arguments);
  printf("Usage:\n");
  printf("$ ./testOscillationFits \n");
  exit(0);
}

int main(int argc,char** argv)
{
  TString inputFile("");
  if(argc!=2) usage(argc);
  TString inputDirPath(argv[1]);
  
  TList *lsFiles;
  //Open tests system directory
  TSystemDirectory testsDir("TestsDir",inputDirPath);
  //Make sure it is an actual system directory
  if(!testsDir.IsDirectory()) {
    printf("%s is not a direcory\n",testsDir.GetName());
    exit(1);
  }
  // Get the list of files contained in the directory
  lsFiles=(TList*)testsDir.GetListOfFiles();
  if(lsFiles){// Make sure there is actually a list
    TIter next(lsFiles);// Instantiate an iterator
    while(TSystemFile *file=(TSystemFile*)next()){// Iterate
      TString fileName=file->GetName();
      //Only get the files that e with ext extension
      if(!(file->IsDirectory()) && (fileName.EndsWith("root"))){
        
        //Get input deltam2 and sin22theta values from the name of the ROOT file
        TObjArray *strArray=fileName.Tokenize("_"); // Tokenize by '-'
        // Get the substring containing dm2 term
        TString dm2Str=((TObjString*)strArray->At(strArray->GetEntries()-2))->GetString();
        double inputMinY=dm2Str.Atof();// convert to floating point integer
        // Get the substring containing s22 term
        TString s22Str=((TObjString*)strArray->Last())->GetString();
        s22Str.ReplaceAll(".root","");
        double inputMinX=s22Str.Atof();// convert to floating point integer
        
        fileName.Prepend("/");
        fileName.Prepend(inputDirPath);// Make sure to prepend the directory name or it will look in the directory from which the program is run.
        TFile *rootFile =(TFile*)TFile::Open(fileName.Data());
        TH2D* chi2Hist=(TH2D*)rootFile->Get("CumulativeChiSquareMap");// Get chi2 histogram
        int x,y,z;
        chi2Hist->GetMinimumBin(x,y,z);
        // Get the x (s22t) value of minimum bin
        double fitMinX=chi2Hist->GetXaxis()->GetBinCenter(x);
        // Get the x (dm2) value of minimum bin
        double fitMinY=chi2Hist->GetYaxis()->GetBinCenter(y);
        // Check if the difference between the input sin22t and dm2 are more than one bin width away, if yes the test has failed.
        double fitGap=TMath::Sqrt(TMath::Power(fitMinX-inputMinX,2)*TMath::Power(fitMinY-inputMinY,2));
        double inputGap=chi2Hist->GetXaxis()->GetBinWidth(x)*chi2Hist->GetYaxis()->GetBinWidth(y);
        if(fitGap<inputGap) {
          printf("PASS: Test passed for bins corresponding to x=%f and y=%f\n",inputMinX,inputMinY);
        }
        else {
          printf("xxxx Test failed for bins corresponding to x=%f and y=%f \n",inputMinX,inputMinY);
          printf("  ---> The values of input and fits for Deltam2 and sin22theta are %f!=%f and %f!=%f respectively\n",inputMinY,fitMinY,inputMinX,fitMinX);
        }
      }
    }
  }
  delete lsFiles;
  return 0;
}

