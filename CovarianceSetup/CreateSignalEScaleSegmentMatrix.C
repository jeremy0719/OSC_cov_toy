#include <iostream>
#include <algorithm>
#include <fstream>
#include <map>

#include "CovarianceTools.H"

#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"

//////////////////////////////////////////////////////////////
// Karin Gilje, January 13, 2015
//
// This function will create an element of the systematic
// covariance matrix, V_{sys}.
// This deals with the segment-by-segment energy scale
// variation of the signal.
//////////////////////////////////////////////////////////////

void ShiftHistogram(TH1F* inputHist, TH1F* toyHist, double scale) {
  toyHist->Scale(0);
  
  // Create a temporary histogram.
  TH1F* tempHist = (TH1F*)inputHist->Clone("tempHist");
  
  double nBins = inputHist->GetNbinsX();
  
  for (int i = 0; i < nBins; i++) {
    // Clear temporary histogram
    tempHist->Scale(0);
    
    double content = inputHist->GetBinContent(i+1);
    double center = inputHist->GetBinCenter(i+1);
    double width = inputHist->GetBinWidth(i+1);
    double lowEdge = inputHist->GetBinLowEdge(i+1);
    double highEdge = lowEdge + width;
    
    // Find the bin that corresponds to the low edge scaled
    double newLowEdge = lowEdge*scale;
    double newLoBin = tempHist->FindBin(newLowEdge);
    
    // Find the bin that corresponds to the high edge scaled
    double newHighEdge = highEdge*scale;
    double newHiBin = tempHist->FindBin(newHighEdge);
    
    // Find the size of the new width for the content
    double newWidth = width * scale;
    
    // Loop through narrow bin range
    for (int j = newLoBin; j < newHiBin+1; j++) {
      // Find Bin Edges
      double low = inputHist->GetBinLowEdge(j);
      double high = low+width;
      
      double tempWidth = 0;
      
      if (low < newLowEdge) tempWidth = high-newLowEdge;
      else if (high > newHighEdge) tempWidth = newHighEdge - low;
      else tempWidth = width;
      
      // Shift the histogram by a scaling factor
      tempHist->SetBinContent(j, content * tempWidth / newWidth);
    }
    
    toyHist->Add(tempHist);
  }
}


void CreateSignalEScaleSegmentMatrix(int nToys = 100) {
  
  // Check if file is attached.
  if (!gFile->IsOpen()) {
    std::cout << "No File Attached." << std::endl;
    std::cout << "Returning ..." << std::endl;
    return;
  }
  
  std::cout << "Retrieve the Energy Spectrum Shape" << std::endl;
  
  TFile* origFile = gFile;
  TString inputName(gFile->GetName());
  std::cout << "Accessing " << inputName << std::endl;
  
  // Load necessary histograms
  TH2F* SegvsENull = (TH2F*)gFile->Get("SegvsE");
  TH2F* LvsENull = (TH2F*)gFile->Get("LvsE");
  TH2F* LToSeg = (TH2F*)gFile->Get("PositionToSegmentFid");
  
  // Create Output File
  TString outputName(inputName);
  outputName.ReplaceAll("Setup", "Cov_Sys_EScaleSegment");
  
  TFile* outputFile = new TFile(outputName, "RECREATE");
  
  // Write the original Histogram to the file
  //SegvsENull->Write();
  //LvsENull->Write();
  //LToSeg->Write();

  // Create vector of Energy Spectra for varying position
  double nEneBins = SegvsENull->GetXaxis()->GetNbins();
  double nSegBins = SegvsENull->GetYaxis()->GetNbins();
  double nPosBins = LvsENull->GetYaxis()->GetNbins();
  
  // Fill the position maps.
  CreatePositionMaps(LToSeg);

  // Vector of Histograms for each segment bin
  std::vector<TH1F*> EnergySpectrumSeg = Create1DProjections(SegvsENull, false);
  
  // Vector of Histograms for each position bin
  std::vector<TH1F*> EnergySpectrumL = Create1DProjections(LvsENull, true);

  std::cout << "Spectra Created." << std::endl;
  
  // Calculate the number of position energy bins for the covariance matrix
  double covMatSegBins = nEneBins * nSegBins;
  double covMatBins = nEneBins * nPosBins;
  
  // Create Covariance Matrix
  TH2F* CovSegMatrix = new TH2F("CovSegMatrix", "Covariance Matrix by Segment; Neutrino Energy (MeV); Neutrino Energy (MeV)",
                             covMatSegBins, -0.5, covMatSegBins-0.5,
                             covMatSegBins, -0.5, covMatSegBins-0.5);
  
  TH2F* CovMatrix = new TH2F("CovMatrix", "Covariance Matrix; Neutrino Energy (MeV); Neutrino Energy (MeV)",
                                covMatBins, -0.5, covMatBins-0.5,
                                covMatBins, -0.5, covMatBins-0.5);
  
  // Set the sigma to use below
  double sigma = 0.1;
  
  // Create a histogram to record the energy scale throw values.
  std::vector<TH1F*> throwHists = CreateThrowHistograms(nSegBins, sigma);
  
  // Create Temporary histograms...
  // These will be modified for each toy.
  std::vector<TH1F*> tempToyHists = CreateToyHistograms(EnergySpectrumSeg);

  // Set the random seed to computer time
  gRandom->SetSeed(0);
 
  std::cout << "Begin Toy Creation... " << std::endl;
  
  // Create Toys
  for (int i = 0; i < nToys; i++) {
    
    if ((i%50) == 0) {
      std::cout << "Working on Toy " << i << std::endl;
    }
    
    // Assume Energy Throw is correlated among position bins. So one value per toy per energy bin.
    std::vector<double> throws;
    throws.clear();
    for (int j = 0; j < nSegBins; j++) {
      double pull = sigma * gRandom->Gaus(0.0, 1.0);
      throwHists[j]->Fill(pull);
      throws.push_back(pull);
    }
    
    // Loop through segments and scale
    for (unsigned int j = 0; j < EnergySpectrumSeg.size(); j++) {
      ShiftHistogram(EnergySpectrumSeg[j], tempToyHists[j], (1+throws[j]));
    }
    
    // Write the first 20 toys to file, for explaination plots.
    // Due to the large number of segments, we will only save one segment.
    if (i < 20) {
      TString name = Form("Toy%dBin%d", i, 14);
      tempToyHists[13]->SetName(name);
      tempToyHists[13]->Write();
    }
    
    // At this point, there are two vectors of histograms,
    // one describes the "nominal" state, the other describes
    // an bin by bin adjusted state.  Next, we use those
    // histograms to calculate a covariance matrix.
    
    // Calculate Covariance per bin
    for (int k = 0; k < covMatSegBins; k++) {
      // The kth bin relates to what energy and position bin?
      int segIndexK = TMath::Floor((double) k / nEneBins);
      int eneIndexK = k - segIndexK * nEneBins;
      
      // Skip the empty bins (happens because of fiducialization)
      double statsK = EnergySpectrumSeg[segIndexK]->GetBinContent(eneIndexK+1);
      if (statsK == 0) continue;
      
      double deviationK = EnergySpectrumSeg[segIndexK]->GetBinContent(eneIndexK+1) - tempToyHists[segIndexK]->GetBinContent(eneIndexK+1);

      for (int l = 0; l < covMatSegBins; l++) {
        // The lth bin relates to what energy and position bin?
        int segIndexL = TMath::Floor((double) l / nEneBins);
        int eneIndexL = l - segIndexL * nEneBins;

        // Skip the empty bins (happens because of fiducialization)
        double statsL = EnergySpectrumSeg[segIndexL]->GetBinContent(eneIndexL+1);
        if (statsL == 0) continue;
        
        double deviationL = EnergySpectrumSeg[segIndexL]->GetBinContent(eneIndexL+1) - tempToyHists[segIndexL]->GetBinContent(eneIndexL+1);
        
        // Calculate Covariance
        double element = deviationK * deviationL;
        CovSegMatrix->Fill(k, l, element);
        
      }
    }

  } // End of loop through toys
  
  // Scale the Segment based covariance matrix by the number of Toys
  CovSegMatrix->Scale(1.0/nToys);
  
  // Scale Covariance matrix by the number of Toys
  CovSegMatrix->Write();
  
  // Convert from Segment based to Position Based
  ConvertCovarianceMatrix(CovSegMatrix, CovMatrix, nEneBins);
  
  // Scale Covariance matrix by the number of Toys
  //CovMatrix->Scale(1.0/nToys);

  // Save Histograms
  for (unsigned int i = 0; i < throwHists.size(); i++) throwHists[i]->Write();
  CovMatrix->Write();
  
  CalculateReducedCovarianceMatrix(CovMatrix, EnergySpectrumL);
}
