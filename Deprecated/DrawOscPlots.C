// C++ Includes
#include <iostream>

// ROOT Includes
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>

// Local Header Files Includes
#include "prospect_style.H"

////////////////////////////////
// Karin Gilje, February 12, 2016
//
// Create plots for technical note
//
/////////////////////////////////


void DrawOscPlots() {
  
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
  TH2F* LvsENull = (TH2F*)gFile->Get("LvsE");
  // 33 is delta m2 of 1.78
  TH2F* LvsE_33 = (TH2F*)gFile->Get("LvsE_33");

  LvsENull->Draw("COLZ");
  gPad->Print("OscPlots/LvsENull.pdf",".pdf");

  TH2F* LvsEOsc = (TH2F*)LvsENull->Clone("LvsEOsc");
  // Assuming sin2 2theta_14 = 0.1
  LvsEOsc->Add(LvsE_33, -0.1);
  
  LvsEOsc->Draw("COLZ");
  gPad->Print("OscPlots/LvsEOsc.pdf",".pdf");
  
  LvsEOsc->Divide(LvsENull);
  
  LvsEOsc->GetZaxis()->SetRangeUser(0.9, 1.0);
  
  LvsEOsc->Draw("COLZ");
  gPad->Print("OscPlots/LvsERatio.pdf",".pdf");

}