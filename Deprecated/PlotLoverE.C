#include <iostream>

#include <TROOT.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>

#include "AnalysisInputs.H"

void PlotLoverE(double LiveTime = 0.41, int massBin = 34, int s22tBin = 28) {
  
  // LOI binning
  // Note that massBin  = 34 and s2tBin = 28 corresponds to delta m^2 = 1.78 eV^2 and sin^2 2theta = 0.09
  
  std::cout << "Running assuming " << LiveTime*100.0 << "% Live Time of the Reactor." << std::endl;
  
  // Ensure a File is loaded
  if (!gFile) {
    std::cout << "File must be attached" << std::endl;
    return;
  }
  
  TFile* origFile = gFile;
  TString inputName(gFile->GetName());
  std::cout << "Accessing " << inputName << std::endl;
  
  TH2F* NullOscillation = (TH2F*)gFile->Get("LvsE");
  TH2F* Deltam2Hist = (TH2F*)gFile->Get(Form("LvsE_%i",massBin));
  
  // Create Output File Name
  TString outputName(inputName);
  outputName.ReplaceAll("Oscillation", "LoverE");
  
  TFile* outputFile = new TFile(outputName, "RECREATE");

  // Setup array of Delta m^2 values to test.  The first value is 0.0 which refers to the unoscillated situation.
  ConstructDeltam2();
  
  // Find the correct massBin
  if (massBin > fNDeltam2) {
    std::cout << "Input mass bin must be less than " << fNDeltam2 << std::endl;
    std::cout << "Returning . . ." << std::endl;
    return;
  }
  std::cout << "Running with delta m^2 = " << fDeltam2[massBin] << std::endl;
  
  // Setup array of Sin2 2Theta values to test.
  ConstructSinSq2Theta();
  
  // Find the correct Sin Sq 2Theta bin
  if (s22tBin > fNSinSq2Theta) {
    std::cout << "Input sin^2 (2theta) bin must be less than " << fNSinSq2Theta << std::endl;
    std::cout << "Returning . . ." << std::endl;
    return;
  }
  std::cout << "Running with sin^2 (2theta) = " << fSinSq2Theta[s22tBin] << std::endl;
  
  
  // Get the number of bins
  int nBinsE = NullOscillation->GetXaxis()->GetNbins();
  int nBinsP = NullOscillation->GetYaxis()->GetNbins();
  
  // Fill the null oscillation histogram
  NullOscillation->SetName("NullOscillation");
  NullOscillation->SetTitle("Null Oscillation; Reconstructed Energy (MeV); Distance from Reactor (m)");
  // Scale Histogram by Reactor livetime
  NullOscillation->Scale(LiveTime);
  NullOscillation->Write();

  // Create Oscillation Histogram
  Deltam2Hist->SetName("Deltam2Hist");
  // Scale to Reactor Live Time
  Deltam2Hist->Scale(LiveTime);
  
  Deltam2Hist->Write();
  
  TH2D* Oscillation = (TH2D*)NullOscillation->Clone("Oscillation");
  // Subtract Oscillation
  Oscillation->Add(Deltam2Hist, -1.0*fSinSq2Theta[s22tBin]);
  
  Oscillation->Write();
  
  double minRange = 0;
  double maxRange = 9;
  int nBins = 72;
  TH1D* LoverEOsc = new TH1D("LoverEOsc",";L/E (m/MeV); Entries", nBins, minRange, maxRange);
  TH1D* LoverENull = new TH1D("LoverENull",";L/E (m/MeV); Entries", nBins, minRange, maxRange);

  double maxPbin = Deltam2Hist->FindLastBinAbove(10,2);
  
  for (int i = 0; i < nBinsE; i++) {
    for (int j = 0; j < maxPbin; j++) {
      // Calculate True Antineutrino Energy of bin
      double EVal = Oscillation->GetXaxis()->GetBinCenter(i+1)+0.8;
      // Calculate the position of the bin
      double PVal = Oscillation->GetYaxis()->GetBinCenter(j+1);
      LoverEOsc->Fill(PVal/EVal, Oscillation->GetBinContent(i+1,j+1));
      LoverENull->Fill(PVal/EVal, NullOscillation->GetBinContent(i+1,j+1));
    }
  }
  
  LoverEOsc->Sumw2();
  LoverENull->Sumw2();
  
  TH1D* LoverERatio = new TH1D("LoverERatio",";L/E (m/MeV); Osc/Null", nBins, minRange, maxRange);

  for (int i = 0; i < LoverEOsc->GetNbinsX(); i++) {
    if (LoverENull->GetBinContent(i+1) < 100.0) continue;
    if (LoverEOsc->GetBinContent(i+1) < 100.0) continue;
    LoverERatio->SetBinContent(i+1,LoverEOsc->GetBinContent(i+1)/LoverENull->GetBinContent(i+1));
    LoverERatio->SetBinError(i+1, TMath::Sqrt(LoverEOsc->GetBinContent(i+1))/LoverENull->GetBinContent(i+1));
  }
  
  LoverERatio->GetYaxis()->SetRangeUser(0.85, 1.02);
  LoverERatio->SetFillColor(kRed);
  LoverERatio->Draw("e2");
  std::cout << "Successfully reached end!" << std::endl;
  
  
  LoverENull->Write();
  LoverEOsc->Write();
  LoverERatio->Write();
  
  
}
