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
#include "TStyle.h"

void usage(int nInputs) {
  printf("Are you serious!, you think this program could run with %i inputs\n",nInputs);
  printf("Usage: ./CalculateSBLDeficit \n");
}

int main(int argc, char** argv) {
  if(argc != 1) {
    (void) argv;
    usage(argc);
    return -1;
  }
  
  SetupProspectStyle();
  
  // Minor style adjustments
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleOffset(1.0,"y");
  
  // Make Canvas
  TCanvas *c1 = new TCanvas("c1","c1",0,0,700,400);
  
  // Load graphs from files
  TFile* fOscSmall = new TFile("SBLOscillations_0.09_0.50.root", "READ");
  if (!fOscSmall->IsOpen() || fOscSmall->IsZombie()) exit(1);
  TGraph* gOscSmall = (TGraph*) fOscSmall->Get("gBaselineRatio");
  gOscSmall->SetName("gOscSmall");
  gOscSmall->SetTitle(";Baseline (m); Oscillated/Unoscillated");
  gOscSmall->SetLineWidth(2);
  gOscSmall->SetLineColor(kBlue);
  
  TFile* fOscMedium = new TFile("SBLOscillations_0.09_1.80.root", "READ");
  if (!fOscMedium->IsOpen() || fOscMedium->IsZombie()) exit(1);
  TGraph* gOscMedium = (TGraph*) fOscMedium->Get("gBaselineRatio");
  gOscMedium->SetName("gOscMedium");
  gOscMedium->SetTitle(";Baseline (m); Oscillated/Unoscillated");
  gOscMedium->SetLineWidth(2);
  gOscMedium->SetLineColor(kRed);
  
  TFile* fOscLarge = new TFile("SBLOscillations_0.09_5.00.root", "READ");
  if (!fOscLarge->IsOpen() || fOscLarge->IsZombie()) exit(1);
  TGraph* gOscLarge = (TGraph*) fOscLarge->Get("gBaselineRatio");
  gOscLarge->SetName("gOscLarge");
  gOscLarge->SetTitle(";Baseline (m); Oscillated/Unoscillated");
  gOscLarge->SetLineWidth(2);
  gOscLarge->SetLineColor(kBlack);
  
  
  // Construct the frame
  TH1F *hFrame = c1->DrawFrame(0,0.9,20,1.05);
  hFrame->GetXaxis()->SetTitle("Baseline (m)");
  hFrame->GetYaxis()->SetTitle("Oscillated/Unoscillated");
  
  // Draw the graphs
  gOscSmall->Draw("CSame");
  gOscMedium->Draw("CSame");
  gOscLarge->Draw("CSame");
  
  // Detector position
  double aDetectorX[] = {9.5, 18};
  double aDetectorEX[] = {2.5, 2};
  // Large 0.9275 - 0.9825
  double aLargeY[] = {(0.9825+0.9275)/2, (0.9825+0.9275)/2};
  double aLargeEY[] = {(0.9825-0.9275)/2, (0.9825-0.9275)/2};
  // Medium 0.9125 - 0.9975
  double aMediumY[] = {(0.9125 + 0.9975)/2, (0.9125 + 0.9975)/2};
  double aMediumEY[] = {(0.9125 - 0.9975)/2, (0.9125 - 0.9975)/2};
  
  // Small 0.91 - 0.94 and 0.97 - 1.0
  double aSmallY[] = {(0.91 + 0.94)/2, (0.97 + 1.0)/2};
  double aSmallEY[] = {(0.91 - 0.94)/2, (0.97 - 1.0)/2};
  
  // Add a legend
  TLegend* leg = new TLegend(0.75,0.7,0.95,0.9);
  leg->SetTextFont(132);
  //leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(gOscSmall,"0.5 eV^{2}", "l");
  leg->AddEntry(gOscMedium,"1.8 eV^{2}", "l");
  leg->AddEntry(gOscLarge,"5.0 eV^{2}", "l");
  
  leg->Draw("same");
  
  // Save the output
  TString outputTitle("BaselineOscNoDetector");
  c1->SaveAs(outputTitle + ".png");
  c1->SaveAs(outputTitle + ".eps");
  c1->SaveAs(outputTitle + ".pdf");
  c1->SaveAs(outputTitle + ".C");
  
  TGraphErrors* gDetectorLarge = new TGraphErrors(1, aDetectorX, aLargeY, aDetectorEX, aLargeEY);
  gDetectorLarge->SetTitle("gDetectorLargeI");
  gDetectorLarge->SetFillStyle(3004);
  gDetectorLarge->SetFillColor(kBlack);
  gDetectorLarge->Draw("Same2");

  TGraphErrors* gDetectorMedium = new TGraphErrors(1, aDetectorX, aMediumY, aDetectorEX, aMediumEY);
  gDetectorMedium->SetTitle("gDetectorMediumI");
  gDetectorMedium->SetFillStyle(3005);
  gDetectorMedium->SetFillColor(kRed);
  gDetectorMedium->Draw("Same2");
  
  TGraphErrors* gDetectorSmall = new TGraphErrors(1, aDetectorX, aSmallY, aDetectorEX, aSmallEY);
  gDetectorSmall->SetTitle("gDetectorSmallI");
  gDetectorSmall->SetFillStyle(3006);
  gDetectorSmall->SetFillColor(kBlue);
  gDetectorSmall->Draw("Same2");
  
  // Save the output
  outputTitle = "BaselineOscPhaseI";
  c1->SaveAs(outputTitle + ".png");
  c1->SaveAs(outputTitle + ".eps");
  c1->SaveAs(outputTitle + ".pdf");
  c1->SaveAs(outputTitle + ".C");
  
  // Redraw
  hFrame->Draw();
  gOscSmall->Draw("CSame");
  gOscMedium->Draw("CSame");
  gOscLarge->Draw("CSame");
  leg->Draw("same");
  
  gDetectorLarge = new TGraphErrors(2, aDetectorX, aLargeY, aDetectorEX, aLargeEY);
  gDetectorLarge->SetTitle("gDetectorLargeII");
  gDetectorLarge->SetFillStyle(3004);
  gDetectorLarge->SetFillColor(kBlack);
  gDetectorLarge->Draw("Same2");
  
  gDetectorMedium = new TGraphErrors(2, aDetectorX, aMediumY, aDetectorEX, aMediumEY);
  gDetectorMedium->SetTitle("gDetectorMediumII");
  gDetectorMedium->SetFillStyle(3005);
  gDetectorMedium->SetFillColor(kRed);
  gDetectorMedium->Draw("Same2");
  
  gDetectorSmall = new TGraphErrors(2, aDetectorX, aSmallY, aDetectorEX, aSmallEY);
  gDetectorSmall->SetTitle("gDetectorSmallII");
  gDetectorSmall->SetFillStyle(3006);
  gDetectorSmall->SetFillColor(kBlue);
  gDetectorSmall->Draw("Same2");
  
  // Save the output
  outputTitle = "BaselineOscPhaseII";
  c1->SaveAs(outputTitle + ".png");
  c1->SaveAs(outputTitle + ".eps");
  c1->SaveAs(outputTitle + ".pdf");
  c1->SaveAs(outputTitle + ".C");

  
}
