///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: June 2017
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
#include "TStyle.h"


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
  
  gStyle->SetOptStat("");
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
  TString outputName("./untracked/LoverE_Plots/");
  inputFileName.ReplaceAll(".root","");
  outputName.Append(inputFileName);
  outputName.Append("_LoverE.root");
  
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
  TH1D* LNullCumulative = new TH1D("LNullCumulative",";L/E (m/MeV); Entries", 31, 6.715596,11.21505);
  TH1D*  LoverENullCumulative= new TH1D("LvsENullCumulative",";L/E (m/MeV); Entries", nBins, minRange, maxRange);
  
  const double threshold=100.0;
  
  for(const int& i : detNumbers)
  {
    detName.Form("LvsENull%i",i);
    
    std::cout << detName.Data() <<std::endl;
    TH2D* NullOscillation = (TH2D*)inputFile->Get(detName.Data());
    
    detName.Form("LvsEDeltam2_%i_%i",i,massBin);
    TH2D* Deltam2Hist = (TH2D*)inputFile->Get(detName.Data());
    
    // Get number of bins
    int nBinsE = NullOscillation->GetXaxis()->GetNbins();
    int nBinsP = NullOscillation->GetYaxis()->GetNbins();
    
    detName.Form("NullOscillation%i",i);
    // Fill the null oscillation histogram
    NullOscillation->SetName(detName);
    NullOscillation->SetTitle("Null Oscillation; Reconstructed Energy (MeV); Distance from Reactor (m)");
    //NullOscillation->Write();
    
    //Deltam2Hist->Write();
    
    detName.Form("Oscillation%i",i);
    TH2D* Oscillation = (TH2D*)NullOscillation->Clone(detName);
    // Subtract Oscillation
    std::cout <<fSinSq2Theta[s22tBin] << std::endl;
    Oscillation->Add(Deltam2Hist, -1.0*fSinSq2Theta[s22tBin]);
    //Oscillation->Write();
    
    
    detName.Form("LoverEOsc_%i",i);
    TH1D* LoverEOsc = new TH1D(detName.Data(),";L/E (m/MeV); Entries", nBins, minRange, maxRange);
    detName.Form("LoverENull_%i",i);
    TH1D* LoverENull = new TH1D(detName.Data(),";L/E (m/MeV); Entries", nBins, minRange, maxRange);
    
    double maxPbin = Deltam2Hist->FindLastBinAbove(10,2);
    double x=0;
    
    LoverENullCumulative->Sumw2();
    
    for (int j = 0; j < nBinsE; j++) {
      for (int k = 0; k < maxPbin; k++) {
        // Calculate True Antineutrino Energy of bin
        double EVal = Oscillation->GetXaxis()->GetBinCenter(j+1)+0.8;
        // Calculate the position of the bin
        double PVal = Oscillation->GetYaxis()->GetBinCenter(k+1);
        LoverEOsc->Fill(PVal/EVal, Oscillation->GetBinContent(j+1,k+1));
        LoverENull->Fill(PVal/EVal, NullOscillation->GetBinContent(j+1,k+1));
        
        LoverEOscCumulative->Fill(PVal/EVal, Oscillation->GetBinContent(j+1,k+1));
        LoverENullCumulative->Fill(PVal/EVal, NullOscillation->GetBinContent(j+1,k+1));
        LNullCumulative->Fill(PVal, NullOscillation->GetBinContent(j+1,k+1));
      }
    }
    
    LoverEOsc->Sumw2();
    LoverENull->Sumw2();
    
    LoverEOscCumulative->Sumw2();
    LoverENullCumulative->Sumw2();
    
    LoverEOscCumulative->Write();
    LoverENullCumulative->Write();
    LNullCumulative->Write();
  }
  outputFile->Close();
}
