///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Jan 2017
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
  printf("Are you serious!, you think this program could run with %i inputs\n",nInputs);
  printf("Usage: ./PlotLoverE massBin s22tBin inputFileName \n");
  exit(1);
}

// Number of bins to be used for sin22theta
// Make sure to have the right number of bins
int fNSinSq2Theta=49;

// Vector of Sin22theta
std::vector<double> fSinSq2Theta;
// Sin22theta bins
double* fSinSq2ThetaBins;

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


int main(int argc,char** argv) {
  
  if(argc!=4) usage(argc);
  
  int massBin = std::stoi(argv[1]);
  int s22tBin = std::stoi(argv[2]);
  TString inputFileName = argv[3];
  
  TFile *inputFile = TFile::Open(inputFileName.Data(),"READ");
  
  if (!(inputFile) && !(inputFile->IsZombie())) {
    printf("File not open\n");
    exit(1);
  }
  
  ConstructSinSq2Theta();
  
  std::vector<int> detNumbers;
  for(int i=0;i<2;i++) {
    int j=0;
    for(j=0;j<3;j++){
      int detNumber = 100*(i+1)+j+1;
      TString detName;
      detName.Form("LvsENull%i",detNumber);
      if(inputFile->FindObjectAny(detName.Data()))
      {
        detName.ReplaceAll("LvsENull","");
        detNumbers.push_back(detName.Atoi());
      }
    }
  }
  
  
  
  TString detName;
  detName.Form("LvsENull%i",detNumbers.at(0));
  // Create Output File Name
  TString outputName(inputFileName);
  outputName.ReplaceAll(".root","_LoverE.root");
  
  TFile* outputFile = new TFile(outputName, "RECREATE");
  
  TH2D* NullOscillationCumulative=new TH2D(*((TH2D*)inputFile->Get(detName.Data())));
  NullOscillationCumulative->Scale(0);
  NullOscillationCumulative->SetName("NullOscillation");
  
  
  double minRange = 0;
  double maxRange = 9;
  int nBins = 72;
//  TH2D* Deltam2Hist = (TH2D*)inputFile->Get(detName.Data());
//  TH2D* Oscillation = (TH2D*)NullOscillation->Clone(detName);
  TH1D* LoverEOscCumulative = new TH1D("LvsEOscCumulative",";L/E (m/MeV); Entries", nBins, minRange, maxRange);
  TH1D* LoverENullCumulative = new TH1D("LvsENullCumulative",";L/E (m/MeV); Entries", nBins, minRange, maxRange);
  // Background L/E histogram
  TH1D* LoverEBKG = new TH1D("LoverEBKG",";L/E (m/MeV); Entries", nBins, minRange, maxRange);
  TH1D* LoverERatioCumulative = new TH1D("LoverERatioCumulative",";L/E (m/MeV); Osc/Null", nBins, minRange, maxRange);
  TH1D* LoverENullErr = new TH1D("LoverERatioErr",";L/E (m/MeV); Osc/Null", nBins, minRange, maxRange);

  
  const double threshold=100.0;
  
  for(const int& i : detNumbers)
  {
    detName.Form("LvsENull%i",i);
    
    TH2D* NullOscillation = (TH2D*)inputFile->Get(detName.Data());
//    NullOscillationCumulative->Add(NullOscillation);
    
    detName.Form("LvsEDeltam2_%i_%i",i,massBin);
    TH2D* Deltam2Hist = (TH2D*)inputFile->Get(detName.Data());
    
    detName.Form("BkgRxOffLvsE%i",i);
    TH2D* BkgRxOffLvsE = (TH2D*)inputFile->Get(detName.Data());
    
    // Get number of bins
    int nBinsE = NullOscillation->GetXaxis()->GetNbins();
    int nBinsP = NullOscillation->GetYaxis()->GetNbins();
    
    detName.Form("NullOscillation%i",i);
    // Fill the null oscillation histogram
    NullOscillation->SetName(detName);
    NullOscillation->SetTitle("Null Oscillation; Reconstructed Energy (MeV); Distance from Reactor (m)");
    NullOscillation->Write();
    
    Deltam2Hist->Write();
    
    
    detName.Form("LvsERef_%i",i);
    TH2D* Oscillation = (TH2D*)inputFile->Get(detName.Data());
    Oscillation->Write();
    
    
    detName.Form("LoverEOsc_%i",i);
    TH1D* LoverEOsc = new TH1D(detName.Data(),";L/E (m/MeV); Entries", nBins, minRange, maxRange);
    detName.Form("LoverENull_%i",i);
    TH1D* LoverENull = new TH1D(detName.Data(),";L/E (m/MeV); Entries", nBins, minRange, maxRange);
    
    double maxPbin = Oscillation->FindLastBinAbove(50,2);
    for (int j = 0; j < nBinsE; j++) {
      for (int k = 0; k < maxPbin; k++) {
        std::cout << nBinsE << "   "  <<nBinsP <<std::endl;
        // Calculate True Antineutrino Energy of bin
        double EVal = Oscillation->GetXaxis()->GetBinCenter(j+1)+0.8;
        // Calculate the position of the bin
        double PVal = Oscillation->GetYaxis()->GetBinCenter(k+1);
        LoverEOsc->Fill(PVal/EVal, Oscillation->GetBinContent(j+1,k+1));
        LoverENull->Fill(PVal/EVal, NullOscillation->GetBinContent(j+1,k+1));
        
        LoverEOscCumulative->Fill(PVal/EVal, Oscillation->GetBinContent(j+1,k+1));
        LoverENullCumulative->Fill(PVal/EVal, NullOscillation->GetBinContent(j+1,k+1));
        LoverEBKG->Fill(PVal/EVal,BkgRxOffLvsE->GetBinContent(j+1,k+1));
        
        double oscErr=TMath::Power(Oscillation->GetBinContent(j+1,k+1),0.5);
        LoverENullErr->Fill(PVal/EVal,oscErr);
      }
    }
    
    LoverEOsc->Sumw2();
    LoverENull->Sumw2();
    
    detName.Form("LoverERatio_%i",i);
    TH1D* LoverERatio = new TH1D(detName.Data(),";L/E (m/MeV); Osc/Null", nBins, minRange, maxRange);
    
    for (int i = 0; i < LoverEOsc->GetNbinsX(); i++) {
      if (LoverENull->GetBinContent(i+1) < threshold) continue;
      if (LoverEOsc->GetBinContent(i+1) < threshold) continue;
      LoverERatio->SetBinContent(i+1,LoverEOsc->GetBinContent(i+1)/LoverENull->GetBinContent(i+1));
      LoverERatio->SetBinError(i+1, TMath::Sqrt(LoverEOsc->GetBinContent(i+1))/LoverENull->GetBinContent(i+1));
    }
    LoverENull->Write();
    LoverEOsc->Write();
    LoverERatio->Write();
  }
  NullOscillationCumulative->Write();
  
  LoverEOscCumulative->Sumw2();
  LoverENullCumulative->Sumw2();
  
  for (int i = 0; i < LoverEOscCumulative->GetNbinsX(); i++) {
    if (LoverENullCumulative->GetBinContent(i+1) < threshold) continue;
    if (LoverEOscCumulative->GetBinContent(i+1) < threshold) continue;
    LoverERatioCumulative->SetBinContent(i+1,LoverEOscCumulative->GetBinContent(i+1)/LoverENullCumulative->GetBinContent(i+1));
    
    double err=TMath::Sqrt(LoverEOscCumulative->GetBinContent(i+1)+LoverENullCumulative->GetBinContent(i+1)+LoverEBKG->GetBinContent(i+1));
    err=err/LoverENullCumulative->GetBinContent(i+1);
    std::cout << err <<std::endl;
    LoverERatioCumulative->SetBinError(i+1, err);
  }
  
  LoverEOscCumulative->Write();
  LoverENullCumulative->Write();
  
  LoverERatioCumulative->Write();
  
  outputFile->Close();
}
