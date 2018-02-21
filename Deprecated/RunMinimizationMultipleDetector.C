// C++ Includes
#include <iostream>
#include <algorithm>
#include <fstream>

// ROOT Includes
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixF.h>
#include <TVectorF.h>
#include "TObjArray.h"
#include "TObjString.h"

// Local Header Files Includes
#include "AnalysisInputs.H"

////////////////////////////////
// Karin Gilje, March, 2016
//
// Minimize samples: Uses a covariance matrix method.
// This combines the inputs from multiple detectors.
//
/////////////////////////////////

// This function takes an input list of TFiles, opens each one and pushes them to a vector.
std::vector<TString> UnpackNames(TString inputList) {
  std::vector<TString> outputVector;
  outputVector.clear();
  
  ifstream inputFile(inputList.Data());
  if (!inputFile.is_open()) {
    std::cout << "File " << inputList << " not found!" << std::endl;
    exit(1);
  }
  
  TString fileName;
  while (!inputFile.eof()) {
    inputFile >> fileName;
    // Need in case of empty line, ensures file name exists.
    if (fileName.IsWhitespace()) continue;
    outputVector.push_back(fileName);
  }
  
  inputFile.close();
  return outputVector;
}

// This function takes an input list of TFiles, opens each one and pushes them to a vector.
std::vector<TFile*> UnpackFiles(std::vector<TString> inputNames) {
  std::vector<TFile*> outputVector;
  outputVector.clear();
  
  for (unsigned int i = 0; i < inputNames.size(); i++) {
    TFile* temp = new TFile(inputNames[i].Data(), "READ");
    outputVector.push_back(temp);
  }
  
  return outputVector;
}

// This function accesses the detector ids embedded in the filenames.
// The input names must have the format: Oscillation_HFIR_DetID_*
std::vector<int> UnpackDetectorIds(std::vector<TString> inputNames) {
  std::vector<int> outputVector;
  outputVector.clear();
  
  for (unsigned int i = 0; i < inputNames.size(); i++) {
    // Remove Directory Information.
    TObjArray* subStrings = inputNames[i].Tokenize("/");
    TString fileName = (((TObjString*)subStrings->At(subStrings->GetEntries() -1))->GetString());
    
    // Look at file name to extract the detector id.
    TObjArray* detectorString = fileName.Tokenize("_");
    
    if (detectorString->GetEntries() < 3) {
      std::cout << "Incorrect input file name: " << inputNames[i] << std::endl;
      std::cout << "Exiting ..." << std::endl;
      return outputVector;
    }
    
    int DetectorType = (((TObjString*)detectorString->At(2))->GetString()).Atoi();
    outputVector.push_back(DetectorType);
  }

  return outputVector;
}

void RunMinimizationMultipleDetector(TString inputList) {
  
  // Check that the input file name ends with .list
  // It should be a list of .root files that you want to analyze together
  if (!inputList.EndsWith(".list")) {
    std::cout << "Incorrect argument structure." << std::endl;
    std::cout << " Argument: " << inputList << " must end in '.list' " << std::endl;
    exit(1);
  }
  
  // Define the number of cycles to use
  double OneCycle = 25./365.;
  const int nDetectors = 4;
  // Explicitly state the number of cycles used for each detector.
  double aRxOnFraction[nDetectors] = {17.0 * OneCycle, 8.0 * OneCycle, 3.0 * OneCycle, 0.0 * OneCycle};
  double aRxOffFraction[nDetectors] = {17.0 * OneCycle, 8.0 * OneCycle, 3.0 * OneCycle, 0.0 * OneCycle};
  
  for (int i = 0; i < nDetectors; i++) {
    std::cout << "Detector " << i << ":" << std::endl;
    std::cout << " Running assuming " << aRxOnFraction[i]*100.0 << "% Reactor On Time." << std::endl;
    std::cout << " Running assuming " << aRxOffFraction[i]*100.0 << "% Cosmic Reactor Off Time." << std::endl;
  }
  
  // Access the TFiles that contain information from the different detectors.
  std::vector<TString> inputNames = UnpackNames(inputList);
  std::vector<TFile*> inputFiles = UnpackFiles(inputNames);
  std::vector<int> detIds = UnpackDetectorIds(inputNames);
  
  // Ensure that all these vectors are the same size.
  if (inputNames.size() != inputFiles.size() || inputNames.size() != detIds.size()) {
    std::cout << "Something went horribly wrong extracting the files from the input list!" << std::endl;
    return;
  }
  
  // As a check, report the files and detectors that are being used in the analysis.
  std::cout << "The number of detectors in the analysis is: " << inputNames.size() << std::endl;
  std::cout << "The the detector ids and their inputs files that are used in the analysis are:" << std::endl;
  for (unsigned int i = 0; i < inputNames.size(); i++) {
    std::cout << "   " << detIds[i] << ": " << inputNames[i] << std::endl;
  }
  
  // Build the arrays of Dm2 and sin22theta
  ConstructDeltam2();
  ConstructSinSq2Theta();
  
  // Load the necessary histograms
  std::vector<TH2F*> LvsENullOsc;
  LvsENullOsc.clear();
  std::vector<TH2F*> LvsEBackgroundRxOn;
  LvsEBackgroundRxOn.clear();
  std::vector<TH2F*> LvsEBackgroundRxOff;
  LvsEBackgroundRxOff.clear();
  std::vector<std::vector<TH2F*> > LvsEDeltam2;
  LvsEDeltam2.clear();
  
  for (unsigned int i = 0; i < inputFiles.size(); i++) {
    inputFiles[i]->cd();
    
    // Access Histograms
    TH2F* tempSigRxOn = (TH2F*)inputFiles[i]->Get("LvsE");
    TH2F* tempBkgRxOn = (TH2F*)inputFiles[i]->Get("LvsEBackground");
    TH2F* tempBkgRxOff = (TH2F*)inputFiles[i]->Get("LvsEBackground");
    
    std::vector<TH2F*> tempDm2;
    tempDm2.clear();
    
    for (int j = 0; j < fNDeltam2; j++) {
      TH2F* LvsETemp = (TH2F*)gFile->Get(TString::Format("LvsE_%d", j));
      // Give Histograms Unique Names
      LvsETemp->SetName(TString::Format("Det%d_LvsE_%d", i, j));
      // Rescale histograms with RxOnFraction
      LvsETemp->Scale(aRxOnFraction[detIds[i]]);
      tempDm2.push_back(LvsETemp);
    }
    
    // Give Histograms Unique Names
    tempSigRxOn->SetName(TString::Format("Det%d_LvsENullOsc_RxOn", i));
    tempBkgRxOn->SetName(TString::Format("Det%d_LvsEBackground_RxOn", i));
    tempBkgRxOff->SetName(TString::Format("Det%d_LvsEBackground_RxOff", i));
    
    // Rescale histograms based on predicted reactor live time
    tempSigRxOn->Scale(aRxOnFraction[detIds[i]]);
    tempBkgRxOn->Scale(aRxOnFraction[detIds[i]]);
    tempBkgRxOff->Scale(aRxOffFraction[detIds[i]]);
    
    // Add Histograms to Vectors
    LvsENullOsc.push_back(tempSigRxOn);
    LvsEBackgroundRxOn.push_back(tempBkgRxOn);
    LvsEBackgroundRxOff.push_back(tempBkgRxOff);
    LvsEDeltam2.push_back(tempDm2);
    
  }
  
  // Access the Covariance Matrices
  
  // Access the Uncorrelated Matrices
  // These are reduced covariance matrices (see documentation).
  std::vector<TH2F*> rCovUCSigMatrix;
  rCovUCSigMatrix.clear();
  std::vector<TH2F*> rCovUCBkgMatrix;
  rCovUCBkgMatrix.clear();
  for (unsigned int d = 0; d < detIds.size(); d++) {
    TFile* tempUCSigFile = new TFile(TString::Format("Cov_Sys_HFIR_%d_UC.root", detIds[d]), "READ");
    if (!tempUCSigFile) {
      std::cout << "The Uncorrelated Signal Matrix File has not been found!!!" << std::endl;
      std::cout << "Returning..." << std::endl;
      return;
    }
    TH2F* tempSigMatrix = (TH2F*)tempUCSigFile->Get("ReducedCovMatrix");
    rCovUCSigMatrix.push_back(tempSigMatrix);
    
    TFile* tempUCBkgFile = new TFile(TString::Format("Cov_Bkg_HFIR_%d_UC.root", detIds[d]), "READ");
    if (!tempUCBkgFile) {
      std::cout << "The Uncorrelated Background Matrix File has not been found!!!" << std::endl;
      std::cout << "Returning..." << std::endl;
      return;
    }
    TH2F* tempBkgMatrix = (TH2F*)tempUCBkgFile->Get("ReducedCovMatrix");
    rCovUCBkgMatrix.push_back(tempBkgMatrix);
  }
  
  // Access the Fully Correlated Matrices
  TFile* tempFCSigFile = new TFile("Cov_Sys_HFIR_FC.root", "READ");
  if (!tempFCSigFile) {
    std::cout << "The Fully Correlated Signal Matrix File has not been found!!!" << std::endl;
    std::cout << "Returning..." << std::endl;
    return;
  }
  TH2F* rCovFCSigMatrix = (TH2F*)tempFCSigFile->Get("ReducedCovMatrix");

  TFile* tempFCBkgFile = new TFile("Cov_Bkg_HFIR_FC.root", "READ");
  if (!tempFCBkgFile) {
    std::cout << "The Fully Correlated Background Matrix File has not been found!!!" << std::endl;
    std::cout << "Returning..." << std::endl;
    return;
  }
  TH2F* rCovFCBkgMatrix = (TH2F*)tempFCBkgFile->Get("ReducedCovMatrix");
  
  // There's probably a simpler way of doing this . . . but this works?
  // Calculate the number of bins per detector
  std::vector<int> detCovBins;
  // Calculate the maximum bin of each detector
  std::vector<int> maxCovBins;
  // Calculate the miniminum bin of each detector
  std::vector<int> minCovBins;
  // The total of covariance matrix bins.
  int covMatBins = 0;
  for (unsigned int i = 0; i < detIds.size(); i++) {
    int eneBins = LvsENullOsc[i]->GetNbinsX();
    int posBins = LvsENullOsc[i]->GetNbinsY();
    int detBins = eneBins * posBins;
    
    minCovBins.push_back(covMatBins);
    covMatBins += detBins;
    detCovBins.push_back(detBins);
    maxCovBins.push_back(covMatBins);
  }
  
  std::cout << "The total number of covariance bins is: " << covMatBins << std::endl;
  
  // Create the Covariance Matrices needed.
  // FC Matrices are Fully Correlated
  // UC Matrices are UnCorrelated
  // PC Matrices are Partially Correlated, maybe needed in Future?
  TMatrixF CovMatrixSigFC(covMatBins, covMatBins);
  TMatrixF CovMatrixSigUC(covMatBins, covMatBins);
  TMatrixF CovMatrixBkgFC(covMatBins, covMatBins);
  TMatrixF CovMatrixBkgUC(covMatBins, covMatBins);
  TMatrixF CovMatrixSigStat(covMatBins, covMatBins);
  TMatrixF CovMatrixBkgStat(covMatBins, covMatBins);
  TMatrixF CovMatrixTot(covMatBins, covMatBins);
  
  // Create Output File
  TString outputName(inputList);
  outputName.ReplaceAll(".list", ".root");
  
  TFile* outputFile = new TFile(outputName, "RECREATE");
  
  // Calculate the number of energy bins (the same for every position!)
  int nEneBins = LvsENullOsc[0]->GetNbinsX();
  
  // Final Covariance matrices used in the analysis.
  TH2F* hSigCovMatrixFC = new TH2F("SigCovMatrixFC", "Fully Correlated Signal Covariance Matrix", covMatBins, -0.5, covMatBins-0.5,covMatBins, -0.5, covMatBins-0.5);
  TH2F* hSigCovMatrixUC = new TH2F("SigCovMatrixUC", "Uncorrelated Signal Covariance Matrix", covMatBins, -0.5, covMatBins-0.5,covMatBins, -0.5, covMatBins-0.5);
  TH2F* hBkgCovMatrixFC = new TH2F("BkgCovMatrixFC", "Fully Correlated Background Covariance Matrix", covMatBins, -0.5, covMatBins-0.5,covMatBins, -0.5, covMatBins-0.5);
  TH2F* hBkgCovMatrixUC = new TH2F("BkgCovMatrixUC", "Uncorrelated Background Covariance Matrix", covMatBins, -0.5, covMatBins-0.5,covMatBins, -0.5, covMatBins-0.5);
  TH2F* hSigStatCovMatrix = new TH2F("SigStatCovMatrix", "Full Statistical Covariance Matrix", covMatBins, -0.5, covMatBins-0.5,covMatBins, -0.5, covMatBins-0.5);
  TH2F* hBkgStatCovMatrix = new TH2F("BkgStatCovMatrix", "Background Statistical Covariance Matrix", covMatBins, -0.5, covMatBins-0.5,covMatBins, -0.5, covMatBins-0.5);
  TH2F* hTotalCovMatrix = new TH2F("TotalCovMatrix", "Full Total Covariance Matrix", covMatBins, -0.5, covMatBins-0.5,covMatBins, -0.5, covMatBins-0.5);
  
  // Build the covariance matrix
  std::cout << "Building the Covariance Matrix." << std::endl;
  for (int i = 0; i < covMatBins; i++) {
    // Access the detector index
    int detIndexI = 0;
    for (unsigned int k = 0; k < detIds.size(); k++) {
      if (i < maxCovBins[k]) {
        detIndexI = k;
        break;
      }
    }
    
    // Access the position bin and energy bin indices.
    int posIndexI = TMath::Floor((double) (i - minCovBins[detIndexI]) / nEneBins);
    int eneIndexI = (i - minCovBins[detIndexI]) - posIndexI * nEneBins;
    
    // Get the total statistics in bin I Reactor On
    double predSigRxOnI = LvsENullOsc[detIndexI]->GetBinContent(eneIndexI+1, posIndexI+1);
    double predBkgRxOnI = LvsEBackgroundRxOn[detIndexI]->GetBinContent(eneIndexI+1, posIndexI+1);
    double predTotRxOnI = predSigRxOnI + predBkgRxOnI;
    
    // Get the total statistics in bin I Reactor Off
    double predBkgRxOffI = LvsEBackgroundRxOff[detIndexI]->GetBinContent(eneIndexI+1, posIndexI+1);
    double predTotRxOffI = predBkgRxOffI;
    
    for (int j = 0; j < covMatBins; j++) {
      // Access the detector index
      int detIndexJ = 0;
      for (unsigned int k = 0; k < detIds.size(); k++) {
        if (j < maxCovBins[k]) {
          detIndexJ = k;
          break;
        }
      }
      
      // Access the position bin and energy bin indices.
      int posIndexJ = TMath::Floor((double) (j - minCovBins[detIndexJ]) / nEneBins);
      int eneIndexJ = (j - minCovBins[detIndexJ]) - posIndexJ * nEneBins;
      
      // Get the total statistics in bin J Reactor On
      double predSigRxOnJ = LvsENullOsc[detIndexJ]->GetBinContent(eneIndexJ+1, posIndexJ+1);
      double predBkgRxOnJ = LvsEBackgroundRxOn[detIndexJ]->GetBinContent(eneIndexJ+1, posIndexJ+1);
      double predTotRxOnJ = predSigRxOnJ + predBkgRxOnJ;
      
      // Get the total statistics in bin J Reactor Off
      double predBkgRxOffJ = LvsEBackgroundRxOff[detIndexJ]->GetBinContent(eneIndexJ+1, posIndexJ+1);
      double predTotRxOffJ = predBkgRxOffJ;
      
      // Build up stats matrix
      // The error on S+B is Sqrt(S+B).
      CovMatrixSigStat(i, j) = TMath::Sqrt(predTotRxOnI*predTotRxOnJ) * (i==j);
      hSigStatCovMatrix->SetBinContent(i+1, j+1, CovMatrixSigStat(i, j));
      
      // Build up signal systematic matrix
      // Need to convert the reduced covariance matrix back into a full covariance matrix
      // Uncorrelated Errors
      if (detIndexI == detIndexJ) {
        CovMatrixSigUC(i, j) = rCovUCSigMatrix[detIndexI]->GetBinContent(i+1, j+1) * predSigRxOnI * predSigRxOnJ;
        hSigCovMatrixUC->SetBinContent(i+1, j+1, CovMatrixSigUC(i, j));
      }
      
      // Correlated Errors
      CovMatrixSigFC(i, j) = rCovFCSigMatrix->GetBinContent(eneIndexI+1, eneIndexJ+1) * predSigRxOnI * predSigRxOnJ;
      hSigCovMatrixFC->SetBinContent(i+1, j+1, CovMatrixSigFC(i, j));
      
      // Build up background systematic matrix
      // Need to convert the reduced covariance matrix back into a full covariance matrix
      // Uncorrelated Errors
      if (detIndexI == detIndexJ) {
        CovMatrixBkgUC(i, j) = rCovUCBkgMatrix[detIndexI]->GetBinContent(i+1, j+1) * predBkgRxOnI * predBkgRxOnJ;
        hBkgCovMatrixUC->SetBinContent(i+1, j+1, CovMatrixBkgUC(i, j));
        
        // Since the background is "measured" when the reactor is off,
        // the statistical error on the background should be determined
        // by the statistics of the reactor off period.
        CovMatrixBkgStat(i, j) = TMath::Sqrt(aRxOnFraction[detIds[detIndexI]]/aRxOffFraction[detIds[detIndexI]] *
                                             predBkgRxOffI *
                                             aRxOnFraction[detIds[detIndexI]]/aRxOffFraction[detIds[detIndexI]] *
                                             predBkgRxOffJ) * (i==j);
        hBkgStatCovMatrix->SetBinContent(i+1, j+1, CovMatrixBkgStat(i, j));
      }
      
      // Correlated Errors
      CovMatrixBkgFC(i, j) = rCovFCBkgMatrix->GetBinContent(eneIndexI+1, eneIndexJ+1) * predBkgRxOnI * predBkgRxOnJ;
      hBkgCovMatrixFC->SetBinContent(i+1, j+1, CovMatrixBkgFC(i, j));
      
    } // End of loop through j bins.
  } // End of loop through i bins.
  
  std::cout << "Inversion Starting..." << std::endl;
  
  // Add Matrices
  // Options: CovMatrixSigStat, CovMatrixBkgStat, CovMatrixSigFC, CovMatrixSigUC, CovMatrixBkgFC, CovMatrixBkgUC
  CovMatrixTot = CovMatrixSigStat + CovMatrixBkgStat + CovMatrixSigFC + CovMatrixSigUC + CovMatrixBkgFC;
  
  // Loop through bins of total covariance matrix and fill the associated histogram.
  // Do this before inversion.
  for (int i = 0; i < covMatBins; i++) {
    for (int j = 0; j < covMatBins; j++) {
      hTotalCovMatrix->SetBinContent(i+1, j+1, CovMatrixTot(i,j));
    }
  }
  
  // Invert matrix
  CovMatrixTot.Invert();
  
  std::cout << "Covariance Matrix Constructed." << std::endl;
  
  // Create Chi Square Histogram
  TH2D* ChiSquareMap = new TH2D("ChiSquareMap","#chi^{2} Values; sin^{2} 2#theta_{14}; #Delta m^{2}_{14}", fNSinSq2Theta, fSinSq2ThetaBins, fNDeltam2-1, fDeltam2Bins);
  
  // Loop through the delta m values.
  for (int m = 0; m < fNDeltam2; m++) {
    std::cout << "Evaluating Delta m2 = " << fDeltam2[m] << std::endl;
    
    // Loop through sin2 2theta values.
    for (int s = 0; s < fNSinSq2Theta; s++) {
      
      // Create the "X" vector
      TVectorF XVector(covMatBins);
      
      // Loop through the detectors.
      for (unsigned int d = 0; d < detIds.size(); d++) {

        // Build temporary histogram
        TH2F* OscHistTemp = (TH2F*)LvsENullOsc[d]->Clone("OscHistTemp");
        OscHistTemp->Add(LvsEDeltam2[d][m], -1.0 * fSinSq2Theta[s]);

        for (int i = minCovBins[d]; i < maxCovBins[d]; i++) {
          // Access the position bin and energy bin indices.
          int posIndexI = TMath::Floor((double) (i - minCovBins[d]) / nEneBins);
          int eneIndexI = (i - minCovBins[d]) - posIndexI * nEneBins;
          
          // Get the total statistics in bin I
          double predSigRxOnI = LvsENullOsc[d]->GetBinContent(eneIndexI+1, posIndexI+1);
          double predBkgRxOnI = LvsEBackgroundRxOn[d]->GetBinContent(eneIndexI+1, posIndexI+1);
          double predTotRxOnI = predSigRxOnI + predBkgRxOnI;
          
          double predBkgRxOffI = LvsEBackgroundRxOff[d]->GetBinContent(eneIndexI+1, posIndexI+1);
          double predTotRxOffI = predBkgRxOffI;
          
          double obsSigRxOnI = OscHistTemp->GetBinContent(eneIndexI+1, posIndexI+1);
          double obsTotRxOnI = obsSigRxOnI + predBkgRxOnI;
          
          // Fill X vector
          XVector(i) = obsSigRxOnI-predSigRxOnI;
          
        } // End of loop through i bins.
      } // End of loop through detectors.
      
      // Find Transpose
      TVectorF XTVector = XVector;
      
      // Perform X^T * C^{-1} * X
      XTVector *= CovMatrixTot;
      double CovarianceChiSquare = XTVector * XVector;
      
      // Fill Chi Square Histogram
      ChiSquareMap->Fill(fSinSq2Theta[s], fDeltam2[m], CovarianceChiSquare);
      
    } // End of loop through sin22theta values.
  } // End of loop through delta m2 values.
  
  // Write covariance matrices used
  hSigStatCovMatrix->Write();
  hBkgStatCovMatrix->Write();
  hSigCovMatrixFC->Write();
  hSigCovMatrixUC->Write();
  hBkgCovMatrixFC->Write();
  hBkgCovMatrixUC->Write();
  
  hTotalCovMatrix->Write();
  
  // Write Chi Square Histogram
  ChiSquareMap->Write();
}
