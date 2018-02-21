#include <iostream>
#include <algorithm>
#include <fstream>

//#include "GiljeStyle.H"
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
// Karin Gilje, February 4, 2016
//
// This function will create an element of the background
// covariance matrix, V_{bkg}.
// This deals with the overall normalization of the background.
// This error comes from our understanding of the variation of
// the background.
//////////////////////////////////////////////////////////////


void CreateBackgroundSpecNormMatrix(int nToys = 100) {
  
  // Check if file is attached.
  if (!gFile->IsOpen()) {
    std::cout << "No File Attached." << std::endl;
    std::cout << "Returning ..." << std::endl;
    return;
  }
  
  std::cout << "Retrieve the Background Energy Spectrum Shape" << std::endl;
  
  TFile* origFile = gFile;
  TString inputName(gFile->GetName());
  std::cout << "Accessing " << inputName << std::endl;
  
  // Load necessary histograms
  TH2F* LvsEBkg = (TH2F*)gFile->Get("BkgL");

  // Create Output File
  TString outputName(inputName);
  outputName.ReplaceAll("Setup", "Cov_Bkg_SpecNorm");
  
  TFile* outputFile = new TFile(outputName, "RECREATE");
  
  // Write the original Histogram to the file
  LvsEBkg->Write();
  
  // Create vector of Energy Spectra for varying position
  double nEneBins = LvsEBkg->GetXaxis()->GetNbins();
  double nPosBins = LvsEBkg->GetYaxis()->GetNbins();
  
  // Vector of Histograms for each position bin
  std::vector<TH1F*> EnergySpectrum = Create1DProjections(LvsEBkg, true);

  std::cout << "Spectra Created." << std::endl;
  
  // Calculate the number of position energy bins for the covariance matrix
  double covMatBins = nEneBins * nPosBins;
  
  // Create Covariance Matrix
  TH2F* CovMatrix = new TH2F("CovMatrix", "Covariance Matrix; Neutrino Energy (MeV); Neutrino Energy (MeV)",
                             covMatBins, -0.5, covMatBins-0.5,
                             covMatBins, -0.5, covMatBins-0.5);
  
  // Set the sigma value to use below
  double sigma = 0.02;
  
  // Create a histogram to record the energy scale throw values.
  std::vector<TH1F*> throwHists = CreateThrowHistograms(1, sigma);
  
  // Create Temporary histograms...
  // These will be modified for each toy.
  std::vector<TH1F*> tempToyHists = CreateToyHistograms(EnergySpectrum);

  // Set the random seed to computer time
  gRandom->SetSeed(0);
 
  std::cout << "Begin Toy Creation... " << std::endl;
  
  // Create Toys
  for (int i = 0; i < nToys; i++) {
    
    if ((i%50) == 0) {
      std::cout << "Working on Toy " << i << std::endl;
    }
    
    // Using a total overall normalization factor (so one element)
    double pull = sigma * gRandom->Gaus(0.0, 1.0);
    throwHists[0]->Fill(pull);

    // Loop through positions
    for (unsigned int j = 0; j < EnergySpectrum.size(); j++) {
      // Loop through energy bins
      for (int k = 0; k < nEneBins; k++) {
        double content = EnergySpectrum[j]->GetBinContent(k+1);
        tempToyHists[j]->SetBinContent(k+1, content * (1.0+pull));
      }
    }

    // Write the first 20 toys to file, for explaination plots.
    // Only save one position bin for comparison
    if (i < 20) {
      TString name = Form("Toy%dBin%d", i, 4);
      tempToyHists[3]->SetName(name);
      tempToyHists[3]->Write();
    }

    // At this point, there are two vectors of histograms,
    // one describes the "nominal" state, the other describes
    // an energy scale adjusted state.  Next, we use those
    // histograms to calculate a covariance matrix.
    
    // Calculate Covariance per bin
    for (int k = 0; k < covMatBins; k++) {
      // The kth bin relates to what energy and position bin?
      int posIndexK = TMath::Floor((double) k / nEneBins);
      int eneIndexK = k - posIndexK * nEneBins;

      // Skip the empty bins (happens because of fiducialization)
      double statsK = EnergySpectrum[posIndexK]->GetBinContent(eneIndexK+1);
      if (statsK == 0) continue;

      double deviationK = EnergySpectrum[posIndexK]->GetBinContent(eneIndexK+1) - tempToyHists[posIndexK]->GetBinContent(eneIndexK+1);

      for (int l = 0; l < covMatBins; l++) {
        // The lth bin relates to what energy and position bin?
        int posIndexL = TMath::Floor((double) l / nEneBins);
        int eneIndexL = l - posIndexL * nEneBins;

        // Skip the empty bins (happens because of fiducialization)
        double statsL = EnergySpectrum[posIndexL]->GetBinContent(eneIndexL+1);
        if (statsL == 0) continue;
        
        double deviationL = EnergySpectrum[posIndexL]->GetBinContent(eneIndexL+1) - tempToyHists[posIndexL]->GetBinContent(eneIndexL+1);
        
        // Calculate Covariance
        double element = deviationK * deviationL;
        CovMatrix->Fill(k, l, element);
        
      }
    }

  } // End of loop through toys
  
  // Scale Covariance matrix by the number of Toys
  CovMatrix->Scale(1.0/nToys);

  // Save Histograms
  for (unsigned int i = 0; i < throwHists.size(); i++) throwHists[i]->Write();
  CovMatrix->Write();

  CalculateReducedCovarianceMatrix(CovMatrix, EnergySpectrum);
}
