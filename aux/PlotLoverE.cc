///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Oct 2017
// Built on top of macros created by K. Gilje
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

void usage(int nInputs)
{
  // Use bins 34 and 27 for Kopp best fit.
  printf("Are you serious!, you think this program could run with %i inputs\n",nInputs);
  printf("Usage: ./PlotLoverE massValue s22tValue inputFileName \n");
  exit(1);
}
/*
// Number of bins to be used for sin22theta
// Make sure to have the right number of bins
int fNSinSq2Theta=49;
//std::vector<double> fSinSq2Theta;// Vector of Sin22theta
double* fSinSq2ThetaBins;// Sin22theta bins

void ConstructSinSq2ThetaBins()
{
  fSinSq2ThetaBins = new double[fNSinSq2Theta+1];
  for (int i = 1; i <= fNSinSq2Theta; i++) {
    fSinSq2ThetaBins[i-1] = fSinSq2Theta.at(i-1)*TMath::Power(10, -0.025);
  }
  fSinSq2ThetaBins[fNSinSq2Theta] = fSinSq2Theta.at(fNSinSq2Theta-1)*TMath::Power(10, 0.025);
}

void ConstructSinSq2Theta()
{
  fSinSq2Theta.clear();
  fSinSq2Theta.push_back(0.0);
  double sinSq2ThetaStepSize = 2.45/fNSinSq2Theta;
  for (int i = 1; i <= fNSinSq2Theta; i++) {
    fSinSq2Theta.push_back(TMath::Power(10, -2.4+ sinSq2ThetaStepSize*(i)));
  }
  ConstructSinSq2ThetaBins();
}
*/
int main(int argc,char** argv) {
  
  if(argc!=4) usage(argc);
  
  double massValue = std::stof(argv[1]);
  double s22tValue = std::stof(argv[2]);
  TString inputFileName = argv[3];
  
  TFile *inputFile = TFile::Open(inputFileName.Data(),"READ");
  
  if (!(inputFile) && !(inputFile->IsZombie())) {
    printf("File not open\n");
    exit(1);
  }
  
  printf("Using deltam2= %f and Sin22theta=%f\n",massValue,s22tValue);
//  ConstructSinSq2Theta();
  
  std::vector<int> detNumbers;
  // Currently we have 2 detectors and at 3 positions, so maxi =2 and maxj=3
  for(int i=0;i<2;i++) {// Iterate through Detector
    int j=0;
    for(j=0;j<3;j++){//Itearate through positions
      int detNumber = 100*(i+1)+j+1;
      TString detName;
      detName.Form("LvsENull%i",detNumber);
      //Check if the detector exists in the supplied ROOT file
      if(inputFile->FindObjectAny(detName.Data()))
      {
        detName.ReplaceAll("LvsENull","");
        detNumbers.push_back(detName.Atoi());
      }
    }
  }
  
  TString histName;
  // Create Output File Name
  TString outputName(inputFileName);
  outputName.ReplaceAll(".root","_LoverE.root");
  
  TFile* outputFile = new TFile(outputName, "RECREATE");
  
  histName.Form("ChiSquareMap%i",detNumbers.at(0));
  TH2D* chi2Hist=(TH2D*)inputFile->Get(histName);
  
  // Get the bin number corresponding to the value provided as the input
  int massBin = chi2Hist->GetYaxis()->FindBin(massValue);
  
  histName.Form("LvsENull%i",detNumbers.at(0));
  TH2D* NullOscillationCumulative=new TH2D(*((TH2D*)inputFile->Get(histName.Data())));
  NullOscillationCumulative->Scale(0);
  NullOscillationCumulative->SetName("NullOscillation");
  
  
  double minRange = 0;
  double maxRange = 9;
  int nBins = 72;
  
  TH1D* LoverEOscCumulative = new TH1D("LvsEOscCumulative",";L/E (m/MeV); Entries", nBins, minRange, maxRange);
  TH1D* LoverEBkgCumulative = new TH1D("LoverEBkgCumulative",";L/E (m/MeV); Entries", nBins, minRange, maxRange);
  TH1D* LoverENullCumulative = new TH1D("LvsENullCumulative",";L/E (m/MeV); Entries", nBins, minRange, maxRange);
  TH1D* LoverERatioStatCumulative = new TH1D("LoverERatioStatCumulative",";L/E (m/MeV); Osc/Null", nBins, minRange, maxRange);
  TH1D* LoverERatioCumulative = new TH1D("LoverERatioCumulative",";L/E (m/MeV); Osc/Null", nBins, minRange, maxRange);
  TH1D* LoverEErrorCumulative = new TH1D("LoverEErrorCumulative",";L/E (m/MeV); Entries", nBins, minRange, maxRange);
  
  const double threshold=0.1;
  
  for(const int& i : detNumbers)
  {
    histName.Form("LvsENull%i",i);
    
    TH2D* NullOscillation = (TH2D*)inputFile->Get(histName.Data());
//    NullOscillationCumulative->Add(NullOscillation);
    
    histName.Form("LvsEDeltam2_%i_%i",i,massBin);
    TH2D* Deltam2Hist = (TH2D*)inputFile->Get(histName.Data());
    
    histName.Form("BkgRxOnLvsE%i",i);
    TH2D* bkgHist = (TH2D*)inputFile->Get(histName.Data());
    
    histName.Form("MinTotalCovMatrix%i",i);
    TH2D* covMatrix=(TH2D*)inputFile->Get(histName);
    
    // Get number of bins
    int nBinsE = NullOscillation->GetXaxis()->GetNbins();
    int nBinsP = NullOscillation->GetYaxis()->GetNbins();
    
    histName.Form("NullOscillation%i",i);
    // Fill the null oscillation histogram
    NullOscillation->SetName(histName);
    NullOscillation->SetTitle("Null Oscillation; Reconstructed Energy (MeV); Distance from Reactor (m)");
    NullOscillation->Write();
    
    Deltam2Hist->Write();
    
    histName.Form("Oscillation%i",i);
    TH2D* Oscillation = (TH2D*)NullOscillation->Clone(histName);
    // Subtract Oscillation
    Oscillation->Add(Deltam2Hist, -1.0*s22tValue);
    Oscillation->Write();
    
    
    histName.Form("LoverEOsc_%i",i);
    TH1D* LoverEOsc = new TH1D(histName.Data(),";L/E (m/MeV); Entries", nBins, minRange, maxRange);
    histName.Form("LoverENull_%i",i);
    TH1D* LoverENull = new TH1D(histName.Data(),";L/E (m/MeV); Entries", nBins, minRange, maxRange);
    histName.Form("LoverEBkg_%i",i);
    TH1D* LoverEBkg = new TH1D(histName.Data(),";L/E (m/MeV); Entries", nBins, minRange, maxRange);
    histName.Form("LoverEError_%i",i);
    TH1D* LoverEError = new TH1D(histName.Data(),";L/E (m/MeV); Entries", nBins, minRange, maxRange);
    
    double maxPbin = Deltam2Hist->FindLastBinAbove(10,2);
    
//    double LoverERatioTotalErrors[nBins];
    for (int j = 0; j < nBinsE; j++) {
      for (int k = 0; k < maxPbin; k++) {
        // Calculate True Antineutrino Energy of bin
        double EVal = Oscillation->GetXaxis()->GetBinCenter(j+1)+0.8;
        // Calculate the position of the bin
        double PVal = Oscillation->GetYaxis()->GetBinCenter(k+1);
        int PEBin=j*nBinsP+k;
        
        LoverEOsc->Fill(PVal/EVal, Oscillation->GetBinContent(j+1,k+1));
        LoverENull->Fill(PVal/EVal, NullOscillation->GetBinContent(j+1,k+1));
        LoverEBkg->Fill(PVal/EVal,bkgHist->GetBinContent(j+1,k+1));
        LoverEError->Fill(PVal/EVal,covMatrix->GetBinContent(PEBin,PEBin));
        
        LoverEOscCumulative->Fill(PVal/EVal, Oscillation->GetBinContent(j+1,k+1));
        LoverEBkgCumulative->Fill(PVal/EVal,bkgHist->GetBinContent(j+1,k+1));
        LoverENullCumulative->Fill(PVal/EVal, NullOscillation->GetBinContent(j+1,k+1));
        LoverEErrorCumulative->Fill(PVal/EVal,covMatrix->GetBinContent(PEBin,PEBin));
      }
    }
    
    LoverEOsc->Sumw2();
    LoverENull->Sumw2();
    
    histName.Form("LoverERatio_%i",i);
    TH1D* LoverERatio = new TH1D(histName.Data(),";L/E (m/MeV); Osc/Null", nBins, minRange, maxRange);
    histName.Form("LoverERatioStat_%i",i);
    TH1D* LoverERatioStat = new TH1D(histName.Data(),";L/E (m/MeV); Osc/Null", nBins, minRange, maxRange);
//    histName.Form("LoverERatioSyst_%i",i);
//    TH1D* LoverERatioSyst = new TH1D(histName.Data(),";L/E (m/MeV); Osc/Null", nBins, minRange, maxRange);
    
    for (int i = 0; i < LoverEOsc->GetNbinsX(); i++) {
      if (LoverENull->GetBinContent(i+1) < threshold) continue;
      LoverERatio->SetBinContent(i+1,LoverEOsc->GetBinContent(i+1)/LoverENull->GetBinContent(i+1));
      LoverERatioStat->SetBinContent(i+1,LoverEOsc->GetBinContent(i+1)/LoverENull->GetBinContent(i+1));
//      LoverERatioSyst->SetBinContent(i+1,LoverEOsc->GetBinContent(i+1)/LoverENull->GetBinContent(i+1));
      LoverERatioStat->SetBinError(i+1, TMath::Sqrt(LoverEOsc->GetBinContent(i+1)+LoverEBkg->GetBinContent(i+1))/LoverENull->GetBinContent(i+1));
      LoverERatio->SetBinError(i+1, TMath::Sqrt(LoverEError->GetBinContent(i+1))/LoverENull->GetBinContent(i+1));
    }
    LoverENull->Write();
    LoverEOsc->Write();
    LoverERatioStat->Write();
    LoverERatio->Write();
    LoverEBkg->Write();
    LoverEError->Write();
  }
  NullOscillationCumulative->Write();
  
  LoverEOscCumulative->Sumw2();
  LoverENullCumulative->Sumw2();
  
  for (int i = 0; i < LoverEOscCumulative->GetNbinsX(); i++) {
    double error=TMath::Sqrt(LoverEOscCumulative->GetBinContent(i+1)+LoverEBkgCumulative->GetBinContent(i+1))/LoverENullCumulative->GetBinContent(i+1);
    // check to make sure the error is not too big
    if (error > threshold) continue;
    // Check to make sure that the null histogram doesn't have 0.
    if (LoverENullCumulative->GetBinContent(i+1) < threshold) continue;
    LoverERatioCumulative->SetBinContent(i+1,LoverEOscCumulative->GetBinContent(i+1)/LoverENullCumulative->GetBinContent(i+1));
    LoverERatioStatCumulative->SetBinContent(i+1,LoverEOscCumulative->GetBinContent(i+1)/LoverENullCumulative->GetBinContent(i+1));
    
    LoverERatioStatCumulative->SetBinError(i+1,error);
    LoverERatioCumulative->SetBinError(i+1,TMath::Sqrt(LoverEErrorCumulative->GetBinContent(i+1))/LoverENullCumulative->GetBinContent(i+1));
  }
  
  LoverEOscCumulative->Write();
  LoverENullCumulative->Write();
  LoverEBkgCumulative->Write();
  
  LoverERatioCumulative->Write();
  LoverERatioStatCumulative->Write();
  
  outputFile->Close();
}
