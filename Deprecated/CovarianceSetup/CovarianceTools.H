#include <vector>
#include <map>


#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"


//////////////////////////////////////////////////////////////
// Karin Gilje, December 11, 2015
//
// This file contains some commonly used functions among the
// covariance matrix calculations.
//////////////////////////////////////////////////////////////

// Map of the segment to the position.
std::map<int, double> SegToLMap;

// This function creates the maps of position to segment.
void CreatePositionMaps(TH2F* inputHist) {

  SegToLMap.clear();
  
  double nSegBins = inputHist->GetXaxis()->GetNbins();
  double nPosBins = inputHist->GetYaxis()->GetNbins();

  for (int i = 0; i < nSegBins; i++) {
    // Skip segments that have no contribution
    double integral = inputHist->Integral(i+1, i+1, 0, -1);
    if (integral == 0) continue;
    
    for (int j = 0; j < nPosBins; j++) {
      // Make sure bin has information
      double content = inputHist->GetBinContent(i+1, j+1);
      if (content == 0) continue;

      SegToLMap[i] = j;
    }
  }
  return;
}

// This function takes a 2D histogram and returns a vector of
// the 1D projections. The flag is true when we want position bins
// and false when we want segment bins.
std::vector<TH1F*> Create1DProjections(TH2F* inputHist, bool flag) {
  int nBins = inputHist->GetNbinsY();
  
  std::vector<TH1F*> outputVector;
  outputVector.clear();
  
  for (int i = 0; i < nBins; i++) {
    double binWidth = inputHist->GetYaxis()->GetBinWidth(i+1);
    double binValue = inputHist->GetYaxis()->GetBinCenter(i+1);
    
    TString title = Form("Position %f; Energy (MeV)", binValue);
    TString name = Form("PositionBin%d", i+1);
    
    if (!flag) {
      title = Form("Segment %i; Energy (MeV)", (int)binValue);
      name = Form("SegmentBin%d", i+1);
    }
    
    // Pull the spectra associated with the position
    TH1F* temp = (TH1F*)inputHist->ProjectionX(name, i+1,  i+1);
    temp->SetTitle(title);
    
    outputVector.push_back(temp);
    
    // Save Histogram to File
    temp->Write();
  }
  
  return outputVector;
}

// This function takes an input of a number of different throws in each pull.
std::vector<TH1F*> CreateThrowHistograms(int nThrowTypes, double sigma) {
  std::vector<TH1F*> outputVector;
  outputVector.clear();
  
  for (int i = 0; i < nThrowTypes; i++) {
    TString title = Form("Throw %d; Throw (%%)", i+1);
    TString name = Form("Throw%d", i+1);
    TH1F* tempThrows = new TH1F(name, title, 50, -5.0*sigma, 5.0*sigma);
    outputVector.push_back(tempThrows);
  }
  
  return outputVector;
}

// This function takes a vector of histograms and clones their structure.
std::vector<TH1F*> CreateToyHistograms(std::vector<TH1F*> inputHistograms) {
  std::vector<TH1F*> outputVector;
  outputVector.clear();

  for (unsigned int i = 0; i < inputHistograms.size(); i++) {
    TString name = Form("tempBin%d", i+1);
    TH1F* tempHist = (TH1F*)inputHistograms[i]->Clone(name);
    tempHist->Scale(0.0);
    outputVector.push_back(tempHist);
  }

  return outputVector;
}


// This function takes the segment based covariance and
// turns it into a position based covariance matrix.
void ConvertCovarianceMatrix(TH2F* inputCovariance, TH2F* outputCovariance, double nEneBins) {
  std::cout << "Coverting Segment based covariance to Position based..." << std::endl;
  
  double covMatSegBins = inputCovariance->GetXaxis()->GetNbins();
  
  std::map<int, double>::iterator mapIt;
  
  // Convert Covariance matrix from energy-segment bins to energy-position bins.
  for (int i = 0; i < covMatSegBins; i++) {
    // Find the segment and energy bin
    int segIndexI = TMath::Floor((double) i / nEneBins);
    int eneIndexI = i - segIndexI * nEneBins;
    
    mapIt = SegToLMap.find(segIndexI);
    if (mapIt == SegToLMap.end()) continue;
    
    // Get the L position of the segment contribution
    double positionI = SegToLMap[segIndexI];
    
    // Calculate the new bin value
    double newI = positionI * nEneBins + eneIndexI;

    for (int j = 0; j < covMatSegBins; j++) {
      // Find the segment and energy bin
      int segIndexJ = TMath::Floor((double) j / nEneBins);
      int eneIndexJ = j - segIndexJ * nEneBins;
      
      mapIt = SegToLMap.find(segIndexJ);
      if (mapIt == SegToLMap.end()) continue;
      
      // Get the L position of the segment contribution
      double positionJ = SegToLMap[segIndexJ];
      
      // Calculate the new bin value
      double newJ = positionJ * nEneBins + eneIndexJ;
      
      double content = inputCovariance->GetBinContent(i+1, j+1);
      
      // Calculate the new covariance...
      outputCovariance->Fill(newI, newJ, content);
    }
  }
  return;
}

// This function takes the covariance matrix and produces the correlation matrix, a matrix of fractional sigma values and a statistics free covariance matrix.
void CalculateReducedCovarianceMatrix(TH2F* inputCovariance, std::vector<TH1F*> inputSpectra) {
  int covMatBins = inputCovariance->GetNbinsX();
  int nEneBins = inputSpectra[0]->GetNbinsX();
  
  // Make Correlation Matrix (visually easier to interpret)
  TH2F* CorrMatrix = new TH2F("CorrMatrix", "Correlation Matrix; Energy Position Bin; Energy Position Bin",
                              covMatBins, -0.5, covMatBins-0.5,
                              covMatBins, -0.5, covMatBins-0.5);
  
  // Create a 1D Histogram of the 1 sigma variation of each bin
  TH1F* SigmaValues = new TH1F("SigmaValues", "1 #sigma Variation; Energy Position Bin",
                               covMatBins, -0.5, covMatBins-0.5);
  
  // Create a reduced covariance matrix with the information from the correlation matrix and the fractional sigmas.
  TH2F* ReducedCovMatrix = new TH2F("ReducedCovMatrix", "Reduced Covariance Matrix; Energy Position Bin; Energy Position Bin",
                                    covMatBins, -0.5, covMatBins-0.5,
                                    covMatBins, -0.5, covMatBins-0.5);
  
  for (int i = 0; i < covMatBins; i++) {
    int posIndexI = TMath::Floor((double) i / nEneBins);
    int eneIndexI = i - posIndexI * nEneBins;
    
    double sigmaI = TMath::Sqrt(inputCovariance->GetBinContent(i+1, i+1));
    double statsI = inputSpectra[posIndexI]->GetBinContent(eneIndexI+1);
    
    if (statsI == 0.0) continue;
    
    // Fill sigma values with fractional sigma
    if (statsI != 0) SigmaValues->Fill(i, sigmaI/statsI);
    
    for (int j = 0; j < covMatBins; j++) {
      int posIndexJ = TMath::Floor((double) j / nEneBins);
      int eneIndexJ = j - posIndexJ * nEneBins;
      
      double sigmaJ = TMath::Sqrt(inputCovariance->GetBinContent(j+1, j+1));
      double statsJ = inputSpectra[posIndexJ]->GetBinContent(eneIndexJ+1);
      
      double correlation = 0.0;
      double reducedCov = 0.0;
      
      if (statsJ == 0.0) continue;
      
      if (sigmaI != 0.0 && sigmaJ != 0.0) {
        // Calculate the correlation element
        correlation = inputCovariance->GetBinContent(i+1, j+1) / (sigmaI*sigmaJ);
        // Calculate the reduced covariance element
        reducedCov = inputCovariance->GetBinContent(i+1, j+1) / (statsI*statsJ);
      }
      
      // Fill the Correlation matrix
      CorrMatrix->Fill(i, j, correlation);
      ReducedCovMatrix->Fill(i, j, reducedCov);
    } // End of loop through j bins
  } // End of loop through i bins

  SigmaValues->Write();
  CorrMatrix->Write();
  ReducedCovMatrix->Write();

}


