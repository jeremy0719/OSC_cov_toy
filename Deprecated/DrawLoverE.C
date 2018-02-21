#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"

#include "prospect_style.H"

void DrawLoverE(){
  
  //Load style
  SetupProspectStyle();
  
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,400);
  //c1->SetGrid();
  
  //Get previously made energy-position graphs
  TFile *f = new TFile("LoverE_HFIR_0_1.00.root");
  TFile *f2 = new TFile("LoverE_HFIR_0_1.00.root");
  
  c1->cd();
  
  TH2F *phaseI = (TH2F*)f->Get("LoverERatio");
  TH2F *phaseII = (TH2F*)f2->Get("LoverERatio");
  
  phaseI->SetLineColor(kBlack);
  phaseII->SetFillColorAlpha(kRed, 0.5);
  
  phaseII->GetXaxis()->SetLabelSize(0.08);
  phaseII->GetYaxis()->SetLabelSize(0.08);
  
  phaseII->GetXaxis()->SetTitleSize(0.08);
  phaseII->GetYaxis()->SetTitleSize(0.08);
  
  phaseII->SetTitleOffset(0.6,"Y");
  
  phaseI->GetXaxis()->SetRangeUser(0.0, 4.8);
  
  phaseI->SetLineWidth(2);
  phaseII->SetLineWidth(2);

  phaseII->GetYaxis()->SetTitle("Osc./Unosc.");
  
  phaseII->Draw("E2");
  phaseI->Draw("same");

  TLine* nullLine = new TLine(0, 1.0, 9, 1.0);
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
  leg->AddEntry(phaseI,"Front Position","l");
  leg->AddEntry(phaseII,"Front Position","f");
  leg->Draw();


  c1->SaveAs("LoverE.eps");
  c1->SaveAs("LoverE.png");
  c1->SaveAs("LoverE.pdf");
  c1->SaveAs("LoverE.C");

}
