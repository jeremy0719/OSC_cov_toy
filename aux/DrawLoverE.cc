#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TLegend.h"

#include "PROSPECTStyle.hh"

void usage(int nInputs)
{
  printf("Are you serious!, you think this program could run with %i inputs\n",nInputs);
  printf("Need an input files\n");
  printf("Usage: ./DrawLoverE inputfile \n");
  exit(1);
}

int main(int argc,char** argv) {
  
  //Load style
  SetupProspectStyle();
  
  if (argc!=2) usage(argc);
  
  TString inputFileName1 = argv[1];
  
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,400);
  //c1->SetGrid();
  
  //Get previously made energy-position graphs
  TFile *f = new TFile(inputFileName1.Data());
  
  c1->cd();
  
  TH2D *phaseI = (TH2D*)f->Get("LoverERatioCumulative");

  
  phaseI->SetLineColor(kRed+2);
  phaseI->SetMarkerStyle(kFullCircle);
  phaseI->SetMarkerColor(kRed+2);
  
  phaseI->GetXaxis()->SetRangeUser(0.0, 4.8);
  
  phaseI->SetLineWidth(2);
  
  phaseI->GetYaxis()->SetRangeUser(0.85,1.02);
  phaseI->GetYaxis()->SetTitle("Osc./Unosc.");
  
  
  phaseI->Draw("E0");
  
  TLine* nullLine = new TLine(0, 1.0, 4.8, 1.0);
  nullLine->SetLineWidth(2);
  nullLine->SetLineStyle(2);
  nullLine->Draw();
  
  TLegend* leg = new TLegend(0.46,0.15,0.92,0.38,"");
  leg->SetTextFont(132);
  leg->SetTextSize(0.07);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);
  
  //leg->SetHeader("#bf{Mass Splitting: 0.00 eV^{2}; Osc. Amplitude: 0.004}");
  leg->SetHeader("#bf{Mass Splitting: 1.78 eV^{2}; Osc. Amplitude: 0.09}");
  //leg->SetHeader("#splitline{#bf{Mass Splitting: 0.10 eV^{2}}}{#bf{Osc. Amplitude: 0.09}}");
  leg->AddEntry(phaseI,"PROSPECT (3 yr)","ELP");
  leg->Draw();
  
  TString histName(inputFileName1);
  histName.ReplaceAll(".root",".pdf");
  c1->SaveAs(histName);
  histName.ReplaceAll(".pdf",".C");
  c1->SaveAs(histName);
  
}
