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
  printf("Usage: ./PlotFluxDeficit \n");
  printf("Optional Usage: ./PlotFluxDeficit sin2theta_14 deltam2_14\n");
}

int main(int argc, char** argv) {
  if(argc != 1 && argc != 3) {
      usage(argc);
      return -1;
  }
  
  // Two neutrino oscillation at short baselines
  double dSin22Theta13 = 0.084;
  double dDeltam213 = 2.5*TMath::Power(10,-3);

  // Best fit for sterile neutrino (Kopp)
  double dSin22Theta14 = 0.09;
  double dDeltam214 = 1.78;
  
  // Allow optional arguments to specify a set of sterile parameters.
  if (argc == 3) {
    dSin22Theta14 = std::atof(argv[1]);
    dDeltam214 = std::atof(argv[2]);
  }
  
  // Create an output root file
    TFile* fOutput = new TFile("./aux/SterileFluxDeficit.root", "RECREATE");
  
  TCanvas* c1 = new TCanvas();
  
  // Energy Threshold: Range starts at 1.8 MeV since it is the threshold for inverse beta decay (mass difference between the proton and the neutron and positron
  double EMin = 1.8;
  double EMax = 10.0;
  
  // Equation 22 in "On the determination of antineutrino spectra from nuclear reactors" by P. Huber
  // PRD Vol. 39 No. 11
  TF1* fNeutrinoFlux = new TF1("fNeutrinoFlux","TMath::Exp([0] + [1]*TMath::Power(x,1) + [2]*TMath::Power(x,2) + [3]*TMath::Power(x,3) + [4]*TMath::Power(x,4) + [5]*TMath::Power(x,5) + [6]*TMath::Power(x,6))", EMin, EMax);
  
  // Values for U-235 taken from Table 3 of above.
  fNeutrinoFlux->SetParameters(4.367, -4.577, 2.100, -0.5294, 0.06186, -0.002777);

  fNeutrinoFlux->Draw();
  c1->Print("NeutrinoFlux.pdf");
  
  // From "Angular distribution of neutron inverse beta decay" by P Vogel and J Beacom.
  // Missing a factor of 10^-42 because I don't care about absolute normalization, just relative normalization.
  TF1* fCrossSection = new TF1("fCrossSection"," [0] * (x - [1]) * TMath::Sqrt(TMath::Power(x - [1],2) - TMath::Power([2],2))", EMin, EMax);
  fCrossSection->SetParameters(0.0952, 1.293, 0.511);
  
  fCrossSection->Draw();
  c1->Print("CrossSection.pdf");
  
  // We just care about the comparison between two shapes osc and null osc so we don't really care about the overall normalization
  
  // Make spectrum shape into histogram
  TH1F* hNull = new TH1F("hNull", "Null Oscillation; Energy (MeV); Rate (arb)", 1000, EMin, EMax);
  
  for (int i = 0; i < hNull->GetNbinsX(); i++) {
    double energy = hNull->GetBinCenter(i+1);
    hNull->SetBinContent(i+1, fNeutrinoFlux->Eval(energy) * fCrossSection->Eval(energy));
  }

  hNull->Draw();
  c1->Print("Spectrum.pdf");
  
  // 3 Neutrino Mixing:
  // Assumes delta m^2_{12} is too small to contribute
  // And that delta m^2_{13} is approximately delta m^2_{23.}
  TF1* fNoSterile = new TF1("fNoSterile", "1-[0]*TMath::Power(TMath::Sin(1.27*[1]*[2]/x),2)", EMin, EMax);
  fNoSterile->SetParameter(0, dSin22Theta13);
  fNoSterile->SetParameter(1, dDeltam213);
  
  // 3+1 Neutrino
  // Assumes delta m^2_{12} is too small to contribute
  // And that delta m^2_{13} is approximately delta m^2_{23}.
  // And that delta m^2_{14} is approximately delta m^2_{24} and delta m^2_{34}.
  TF1* fSterile = new TF1("fSterile", "1-[0]*TMath::Power(TMath::Sin(1.27*[1]*[2]/x),2)-[3]*TMath::Power(TMath::Sin(1.27*[4]*[2]/x),2)", EMin, EMax);
  fSterile->SetParameter(0, dSin22Theta14);
  fSterile->SetParameter(1, dDeltam214);
  fSterile->SetParameter(3, dSin22Theta13);
  fSterile->SetParameter(4, dDeltam213);
  
  // Do baselines from 1 to 2000 by log steps.
  int nSteps = 100;
  double dRatio3[nSteps];
  double dRatio4[nSteps];
  double dDistance[nSteps];
  for (int i = 0; i < nSteps; i++) {
    dDistance[i] = TMath::Power(10, 3.3/nSteps * i);
    
    // Set Baseline
    fNoSterile->SetParameter(2, dDistance[i]);
    fSterile->SetParameter(2, dDistance[i]);
    // Loop through histogram and oscillate
    TH1F* h3Osc = (TH1F*) hNull->Clone("h3Osc");
    TH1F* h4Osc = (TH1F*) hNull->Clone("h4Osc");
    
    for (int i = 0; i < hNull->GetNbinsX(); i++) {
      double prob3 = fNoSterile->Eval(hNull->GetBinCenter(i+1));
      h3Osc->SetBinContent(i+1, hNull->GetBinContent(i+1) * prob3);
      double prob4 = fSterile->Eval(hNull->GetBinCenter(i+1));
      h4Osc->SetBinContent(i+1, hNull->GetBinContent(i+1) * prob4);
    }
    
    dRatio3[i] = h3Osc->Integral()/hNull->Integral();
    dRatio4[i] = h4Osc->Integral()/hNull->Integral();
    
    h3Osc->Delete();
    h4Osc->Delete();
  }
  
  TGraph* gRatio3 = new TGraph(nSteps, dDistance, dRatio3);
  gRatio3->SetTitle("3 Neutrino Mixing; Baseline (m); N_{osc}/N_{null}");
  gRatio3->SetName("gRatio3");
  c1->SetLogx();
  gRatio3->Draw("AC");
  c1->Print("Ratio3.pdf");
  
  TGraph* gRatio4 = new TGraph(nSteps, dDistance, dRatio4);
  gRatio4->SetTitle("3+1 Neutrino Mixing; Baseline (m); N_{osc}/N_{null}");
  gRatio4->SetName("gRatio4");
  gRatio4->Draw("AC");
  c1->Print("Ratio4.pdf");
  
  
  // Save output to a root file
  fCrossSection->Write();
  fNeutrinoFlux->Write();
  hNull->Write();
  gRatio3->Write();
  gRatio4->Write();
  
  fOutput->Write();
  fOutput->Close();
  
} // End of main
