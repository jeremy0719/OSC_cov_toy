#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <map>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TList.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "PROSPECTStyle.hh"

void usage(int nInputs) {
  printf("Are you serious!, you think this program could run with %i inputs\n",nInputs);
  printf("Usage: ./CalculateSBLDeficit \n");
  printf("Optional Usage: ./CalculateSBLDeficit sin2theta_14 deltam2_14\n");
}

int main(int argc, char** argv) {
  if(argc != 1 && argc != 3) {
      usage(argc);
      return -1;
  }

  // Best fit for sterile neutrino (Kopp)
  double dSin22Theta14 = 0.09;
  double dDeltam214 = 1.78;
  
  // Allow optional arguments to specify a set of sterile parameters.
  if (argc == 3) {
    dSin22Theta14 = std::atof(argv[1]);
    dDeltam214 = std::atof(argv[2]);
  }
  
  // Create an output root file
  // Base output names on input parameters
  TFile* fOutput = new TFile(TString::Format("SBLOscillations_%.2f_%.2f.root", dSin22Theta14, dDeltam214), "RECREATE");
  
  TCanvas* c1 = new TCanvas();
  
  
  // 2 Neutrino Mixing:
  // Only assume sterile mixing
  // Set the energy at 4.0 MeV
  TF1* fOscillation = new TF1("fOscillation", "1-[0]*TMath::Power(TMath::Sin(1.27*[1]*x/[2]),2)", -5.0, +25.0);
  fOscillation->SetParameter(0, dSin22Theta14);
  fOscillation->SetParameter(1, dDeltam214);
  fOscillation->SetParameter(2, 4.0);
  
  // Average out baselines over 1.0 meter reactor
  TH1F* hReactorBaselines = new TH1F("hReactorBaselines", "Reactor Baselines; Range(m); Uniform events", 50, -.5, .5);
  for (int i = 0; i < hReactorBaselines->GetNbinsX(); i++) {
    hReactorBaselines->SetBinContent(i+1, 1.0);
  }
  
  // Do baselines from 0 to 20.
  int nSteps = 100;
  double dRatio[nSteps];
  double dDistance[nSteps];
  for (int i = 0; i < nSteps; i++) {
    dDistance[i] = 20.0/nSteps * i;
    
    // Loop through histogram and oscillate
    TH1F* hOscillatedBaseline = (TH1F*) hReactorBaselines->Clone("hOscillatedBaseline");
    
    for (int j = 0; j < hReactorBaselines->GetNbinsX(); j++) {
      double tempBaseline = dDistance[i] + hReactorBaselines->GetBinCenter(j+1);
      
      // Evaluate oscillation at baseline for 4 MeV.
      double dOscillationProbability = fOscillation->Eval(tempBaseline);
      
      // Reweight the baseline bins by the probability of oscilation at 4 MeV
      hOscillatedBaseline->SetBinContent(j+1, dOscillationProbability);
    }
    
    // Take the ratio of the reweighted osc to the default.
    dRatio[i] = hOscillatedBaseline->Integral()/hReactorBaselines->Integral();
    if ( i == 10 ) hOscillatedBaseline->Write();

    // Delete the histograms to remove memory issues
    hOscillatedBaseline->Delete();
  }
  
  // Create a graph
  TGraph* gBaselineRatio = new TGraph(nSteps, dDistance, dRatio);
  gBaselineRatio->SetTitle("Sterile Neutrino Mixing; Baseline (m); N_{osc}/N_{null}");
  gBaselineRatio->SetName("gBaselineRatio");
  gBaselineRatio->Draw("AC");
  c1->Print("SBLRatio.pdf");
  
  // Save output to a root file
  gBaselineRatio->Write();
  hReactorBaselines->Write();
  
  fOutput->Write();
  fOutput->Close();
  
} // End of main
