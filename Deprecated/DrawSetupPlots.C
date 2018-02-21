// C++ Includes
#include <iostream>

// ROOT Includes
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>

// Local Header Files Includes
//#include "CovarianceSetup/GiljeStyle.H"
#include "prospect_style.H"

////////////////////////////////
// Karin Gilje, February 12, 2016
//
// Create plots for technical note
//
/////////////////////////////////


void DrawSetupPlots() {
  
  // Load Style
  //SetGlobalStyle();
  SetupProspectStyle();
  SetGradientPalette();
  
  // Ensure a File is loaded
  if (!gFile) {
    std::cout << "File must be attached" << std::endl;
    return;
  }
  
  // Load necessary histograms
  TH2F* BaselineToSegment = (TH2F*)gFile->Get("BaselineToSegment");
  TH2F* PositionToSegment = (TH2F*)gFile->Get("PositionToSegment");
  TH2F* TrueEToRecoEFine = (TH2F*)gFile->Get("TrueEToRecoEFine");
  TH2F* TrueEToRecoE = (TH2F*)gFile->Get("TrueEToRecoE");
  TH2F* LvsENull = (TH2F*)gFile->Get("LvsE");
  TH2F* SegvsENull = (TH2F*)gFile->Get("SegvsE");
  
  // Load background histograms to pass to the next stage
  TH2F* BkgL = (TH2F*)gFile->Get("BkgL");
  TH2F* BkgSeg = (TH2F*)gFile->Get("BkgSeg");
  
  BaselineToSegment->Draw("COLZ");
  gPad->Print("SetupPlots/BaselineToSegment.pdf",".pdf");
  
  PositionToSegment->Draw("COLZ");
  gPad->Print("SetupPlots/PositionToSegment.pdf",".pdf");

  TrueEToRecoEFine->Draw("COLZ");
  gPad->Print("SetupPlots/TrueEToRecoEFine.pdf",".pdf");

  TrueEToRecoE->Draw("COLZ");
  gPad->Print("SetupPlots/TrueEToRecoE.pdf",".pdf");

  LvsENull->Draw("COLZ");
  gPad->Print("SetupPlots/LvsENull.pdf",".pdf");

  SegvsENull->Draw("COLZ");
  gPad->Print("SetupPlots/SegvsENull.pdf",".pdf");

  BkgL->Draw("COLZ");
  gPad->Print("SetupPlots/BkgL.pdf",".pdf");

  BkgSeg->Draw("COLZ");
  gPad->Print("SetupPlots/BkgSeg.pdf",".pdf");

  // Reset Right Margin for non-COLZ plots
  gPad->SetRightMargin(0.05);
  gPad->Modified();
  gPad->Update();
  
  TH1F* Bkg1D = (TH1F*)gFile->Get("Bkg1D");

  Bkg1D->SetLineColor(kBlue+2);
  Bkg1D->Draw();
  gPad->Print("SetupPlots/Bkg1D.pdf",".pdf");

}