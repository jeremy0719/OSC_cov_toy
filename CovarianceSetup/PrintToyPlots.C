#include <iostream>

#include "GiljeStyle.H"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TString.h"

TCanvas* gCanvas;


void Print2DPlots(TH2F* inputHist, TString inputName,
                  TString inputTitle, TH1F* EnergyBins,
                  bool segment) {
  // Get the location of the first filled bin
  double firstBinX = inputHist->FindFirstBinAbove(0, 1);
  double firstBinY = inputHist->FindFirstBinAbove(0, 2);
  
  // Get the number of Energy bins
  double nEBins = EnergyBins->GetXaxis()->GetNbins();
  // Setup the energy binning
  double xBins[(int) (nEBins+1)];
  for (int i = 0; i < nEBins+1; i++) {
    xBins[i] = EnergyBins->GetBinLowEdge(i+1);
  }
  
  // Setup the long and short versions of the string describing the type.
  // This will be used in the naming conventions
  TString typelong;
  TString typeshort;
  if (segment) {
    typeshort = "Seg";
    typelong = "Segment";
  }
  else {
    typelong = "Position";
    typeshort = "Pos";
  }
  
  
  /////////////////////////////////////////////////////
  //  Make Matrix Histogram for First Position Bin
  /////////////////////////////////////////////////////
  
  TH2F* OnePos = new TH2F("OnePos","; Reconstructed Energy (MeV); Reconstructed Energy (MeV)", nEBins, xBins, nEBins, xBins);
  
  OnePos->SetTitle(TString::Format("%s Matrix for First Filled %s Bin",
                                   inputTitle.Data(), typelong.Data()));
  
  for (int i = 0; i < nEBins; i++) {
    for (int j = 0; j < nEBins; j++) {
      double content = inputHist->GetBinContent(firstBinX+i, firstBinY+j);
      OnePos->SetBinContent(i+1, j+1, content);
    }
  }
  
  OnePos->SetStats(kFALSE);
  if (inputName.Contains("Corr")) OnePos->GetZaxis()->SetRangeUser(-1.0, 1.0);
  
  OnePos->Draw("COLZ");

  gCanvas->Print(TString::Format("plots/%sOne%sBin.png",
                                 inputName.Data(), typeshort.Data()), "png");
  gCanvas->Print(TString::Format("plots/%sOne%sBin.pdf",
                                 inputName.Data(), typeshort.Data()), "pdf");
  gCanvas->Print(TString::Format("plots/%sOne%sBin.C",
                                 inputName.Data(), typeshort.Data()), "C");
  
  OnePos->Delete();
  
  /////////////////////////////////////////////////////
  //  Make Histogram for two Position Bins
  /////////////////////////////////////////////////////
  
  TH2F* TwoPos = new TH2F("TwoPos","; Energy Position Bin; Energy Position Bin",
                                 nEBins*2, -0.5, nEBins*2-0.5,
                                 nEBins*2, -0.5, nEBins*2-0.5);
  
  TwoPos->SetTitle(TString::Format("%s Matrix for First Two Filled %s Bins; Energy %s Bin; Energy %s Bin",
                           inputTitle.Data(), typelong.Data(), typelong.Data(), typelong.Data()));
                          
  for (int i = firstBinX; i < nEBins*2+firstBinX; i++) {
    for (int j = firstBinY; j < nEBins*2+firstBinY; j++) {
      double content = inputHist->GetBinContent(i, j);
      TwoPos->Fill(i-firstBinX, j-firstBinX, content);
    }
  }
  
  TwoPos->SetStats(kFALSE);
  if (inputName.Contains("Corr")) TwoPos->GetZaxis()->SetRangeUser(-1.0, 1.0);
  
  TwoPos->Draw("COLZ");
  
  gCanvas->Print(TString::Format("plots/%sTwo%sBin.png",
                                 inputName.Data(), typeshort.Data()), "png");
  gCanvas->Print(TString::Format("plots/%sTwo%sBin.pdf",
                                 inputName.Data(), typeshort.Data()), "pdf");
  gCanvas->Print(TString::Format("plots/%sTwo%sBin.C",
                                 inputName.Data(), typeshort.Data()), "C");
  
  TwoPos->Delete();
  
  /////////////////////////////////////////////////////
  //  Make Correlation Matrix Histogram
  /////////////////////////////////////////////////////
  
  inputHist->UseCurrentStyle();
  inputHist->SetStats(kFALSE);
  if (inputName.Contains("Corr")) inputHist->GetZaxis()->SetRangeUser(-1.0, 1.0);
  
  inputHist->SetTitle(TString::Format("%s Matrix; Energy %s Bin; Energy %s Bin",
                           inputTitle.Data(), typelong.Data(), typelong.Data()));
  inputHist->Draw("COLZ");
  
  gCanvas->Print(TString::Format("plots/%sMatrix.png",
                                 inputName.Data()), "png");
  gCanvas->Print(TString::Format("plots/%sMatrix.pdf",
                                 inputName.Data()), "pdf");
  gCanvas->Print(TString::Format("plots/%sMatrix.C",
                                 inputName.Data()), "C");
  
  return;
  
}

void PrintToyPlots() {
  SetGlobalStyle();

  // Create Canvas
  gCanvas = new TCanvas("gCanvas", "gCanvas", 600, 600);
  
  // Check that a file is attached.
  if (!gFile->IsOpen()) {
    std::cout << "No File Attached." << std::endl;
    std::cout << "Returning ..." << std::endl;
    return;
  }

  // Get histogram with energy binning.
  TH1F* EnergyBins = (TH1F*) gFile->Get("PositionBin1");
  
  // Get nominal histogram
  TH1F* nominal = (TH1F*) gFile->Get("SegmentBin14");
  bool segment = true;
  int binNum = 14;
  TString binType = "Segment";
  
  // If the histogram is position based, switch to position based histogram.
  if (!nominal) {
    segment = false;
    nominal = (TH1F*) gFile->Get("PositionBin4");
    binNum = 4;
    binType = "Position";
  }
  
  // Get Toy spectra Histograms
  std::vector<TH1F*> ToyHists;
  ToyHists.clear();

  for (int i = 0; i < 20; i++) {
    TString name = Form("Toy%dBin%d",i, binNum);
    TH1F* temp = (TH1F*) gFile->Get(name);
    if (!temp) {
      std::cout << "Histogram of name: " << name << " not found." << std::endl;
      std::cout << "Quit Loop..." << std::endl;
      break;
    }
    ToyHists.push_back(temp);
  }
  
  // Get sigma histogram
  TH1F* SigmaValues = (TH1F*) gFile->Get("SigmaValues");

  // Get correlation histogram
  TH2F* CorrMatrix = (TH2F*) gFile->Get("CorrMatrix");
  
  // Get correlation histogram
  TH2F* CovMatrix = (TH2F*) gFile->Get("CovMatrix");

  // Get correlation histogram
  TH2F* ReducedCov = (TH2F*) gFile->Get("ReducedCovMatrix");
  
  // Get the distribution of the throws
  TH1F* ThrowValues = (TH1F*) gFile->Get(TString::Format("Throw%d", binNum));
  if (!ThrowValues) ThrowValues = (TH1F*) gFile->Get("Throw1");
  
  // Get the segment based covariance matrix.
  // This is a HUGE matrix, don't look at it in TBrowser!
  TH2F* CovSegMatrix = (TH2F*) gFile->Get("CovSegMatrix");

  /////////////////////////////////////////////////////
  //  Make Spectra Comparison Histograms
  /////////////////////////////////////////////////////
  
  // Find the maximum of the histograms
  double maximum = 0.0;
  for (unsigned int i = 0; i < ToyHists.size(); i++) {
    if (maximum < ToyHists[i]->GetMaximum())
      maximum = ToyHists[i]->GetMaximum();
  }

  nominal->UseCurrentStyle();
  
  nominal->SetTitle(TString::Format("Toy spectra for %s Bin %d", binType.Data(), binNum));
  nominal->GetYaxis()->SetRangeUser(0.0, maximum*(1.05));
  nominal->SetStats(kFALSE);
  nominal->SetLineColor(kMagenta+2);
  nominal->SetLineWidth(2);

  nominal->Draw();

  for (unsigned int i = 0; i < ToyHists.size(); i++) {
    ToyHists[i]->SetStats(kFALSE);
    ToyHists[i]->Draw("same");
  }

  nominal->Draw("same");

  gCanvas->Print("plots/ToySpectra.png", "png");
  gCanvas->Print("plots/ToySpectra.pdf", "pdf");
  gCanvas->Print("plots/ToySpectra.C", "C");
  
  /////////////////////////////////////////////////////
  //  Make Fractional Ratio Histograms
  /////////////////////////////////////////////////////
  
  
  TH1F* nominalRatio = (TH1F*)nominal->Clone("nominalRatio");
  
  nominalRatio->Add(nominal, -1.0);
  nominalRatio->Divide(nominal);
  
  std::vector<TH1F*> ToyRatios;
  ToyRatios.clear();
  
  maximum = 0.0;
  double minimum = 0.0;
  for (unsigned int i = 0; i < ToyHists.size(); i++) {
    TString name = Form("Toy%dBin%dRatio",i, binNum);
    TH1F* temp = (TH1F*)ToyHists[i]->Clone(name);
    temp->Add(nominal, -1.0);
    temp->Divide(nominal);
    
    ToyRatios.push_back(temp);
    
    if (maximum < temp->GetMaximum())
      maximum = temp->GetMaximum();
    if (minimum > temp->GetMinimum())
      minimum = temp->GetMinimum();
  }
  
  nominalRatio->UseCurrentStyle();
  
  nominalRatio->SetTitle(TString::Format("Toy spectra Ratio for %s Bin %d", binType.Data(), binNum));
  nominalRatio->GetYaxis()->SetRangeUser(minimum*(1.05), maximum*(1.05));
  nominalRatio->SetStats(kFALSE);
  nominalRatio->SetLineColor(kMagenta+2);
  nominalRatio->SetLineWidth(2);
  
  nominalRatio->Draw();
  
  for (unsigned int i = 0; i < ToyRatios.size(); i++) {
    ToyRatios[i]->SetStats(kFALSE);
    ToyRatios[i]->Draw("same");
  }
  
  nominalRatio->Draw("same");
  
  gCanvas->Print("plots/ToySpectraRatio.png", "png");
  gCanvas->Print("plots/ToySpectraRatio.pdf", "pdf");
  gCanvas->Print("plots/ToySpectraRatio.C", "C");
  
  /////////////////////////////////////////////////////
  //  Make Throw Values histogram
  /////////////////////////////////////////////////////
  
  ThrowValues->UseCurrentStyle();
  ThrowValues->SetStats(false);
  ThrowValues->SetLineColor(kMagenta+2);
  ThrowValues->SetLineWidth(2);
  
  ThrowValues->SetTitle(Form("Throw for %s Bin %d", binType.Data(), binNum));
  
  ThrowValues->Draw();
  
  gCanvas->Print("plots/ThrowValues.png", "png");
  gCanvas->Print("plots/ThrowValues.pdf", "pdf");
  gCanvas->Print("plots/ThrowValues.C", "C");
  
  /////////////////////////////////////////////////////
  //  Make Fractional Sigma Histogram
  /////////////////////////////////////////////////////
  
  TH1F* fracSigma = (TH1F*)nominal->Clone("fracSigma");
  
  fracSigma->Scale(0.0);
  
  double nEBins = fracSigma->GetNbinsX();

  // Get the location of the first filled bin
  double firstBinX = CovMatrix->FindFirstBinAbove(0, 1);
  double firstBinY = CovMatrix->FindFirstBinAbove(0, 2);

  for (int i = 0; i < nEBins+1; i++) {
    double content = SigmaValues->GetBinContent(firstBinX+i);
    fracSigma->SetBinContent(i+1, content);
  }

  fracSigma->SetTitle("Fractional Sigma");
  fracSigma->Draw();
  
  gCanvas->Print("plots/FractionalSigma.png", "png");
  gCanvas->Print("plots/FractionalSigma.pdf", "pdf");
  gCanvas->Print("plots/FractionalSigma.C", "C");
  
  /////////////////////////////////////////////////////
  //  Make 2D Matrices
  /////////////////////////////////////////////////////
  
  // Print Correlation Matrix
  Print2DPlots(CorrMatrix, "Corr", "Correlation", nominal, false);
  
  // Change palette for other histograms
  SetBlackBodyPalette();
  
  // Print Covariance Matrices
  Print2DPlots(CovMatrix, "Cov", "Covariance", nominal, false);
  Print2DPlots(ReducedCov, "RCov", "Reduced Covariance", nominal, false);

  if (segment) Print2DPlots(CovSegMatrix, "CovSeg", "Segment-based Covariance", nominal, true);
  
  return;
}
