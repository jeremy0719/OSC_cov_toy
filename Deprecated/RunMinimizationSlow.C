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

// Local Header Files Includes
#include "AnalysisInputs.H"

////////////////////////////////
// Karin Gilje, October 14, 2015
//
// Minimize samples: Uses a covariance matrix method.
//
/////////////////////////////////

void RunMinimizationSlow(double RxOnFraction = fRxOnFraction) {
  
  std::cout << "Running assuming " << RxOnFraction*100.0 << "% Live Time of the Reactor." << std::endl;
  
  // Ensure a File is loaded
  if (!gFile) {
    std::cout << "File must be attached" << std::endl;
    return;
  }
  
  TFile* origFile = gFile;
  TString inputName(gFile->GetName());
  std::cout << "Accessing " << inputName << std::endl;
  
  // Load necessary histograms
  TH2F* LvsENull = (TH2F*)gFile->Get("LvsE");
  TH2F* LvsEBackground = (TH2F*)gFile->Get("LvsEBackground");
  
  // Rescale histograms based on predicted reactor live time
  LvsENull->Scale(RxOnFraction);
  LvsEBackground->Scale(RxOnFraction);
  
  // Create a vector of delta m2 histograms.
  std::vector<TH2F*> LvsEDeltam2;
  LvsEDeltam2.clear();
  
  ConstructDeltam2();
  ConstructSinSq2Theta();
  
  for (int i = 0; i < fNDeltam2; i++) {
    TH2F* LvsETemp = (TH2F*)gFile->Get(TString::Format("LvsE_%d", i));
    
    // Rescale histograms with RxOnFraction
    LvsETemp->Scale(RxOnFraction);
    LvsEDeltam2.push_back(LvsETemp);
  }
  
  // Access the covariance matrix
  TFile* CovMatFile = new TFile("Cov_Sys_HFIR_0_1.00.root", "READ");
  
  // Pull covariance matrix from file
  // This is the "reduced" covariance matrix.
  // We need to divide by the number of events for each bin.
  // There are three parts, signal, background and statistical
  // The statistical matrix is created in this file.
  // The other two should be calculated beforehand.
  TH2F* ReducedSigMatrix = (TH2F*)CovMatFile->Get("ReducedCovMatrix");
  TH2F* ReducedBkgMatrix = (TH2F*)CovMatFile->Get("ReducedCovMatrix");

  // Find the total number of bins.
  double nEneBins = LvsENull->GetXaxis()->GetNbins();
  double nPosBins = LvsENull->GetYaxis()->GetNbins();
 
  // Find the fiducialized range.
  double firstPosBin = LvsENull->FindFirstBinAbove(1, 2);
  double lastPosBin = LvsENull->FindLastBinAbove(1, 2);
  double nFidBins = lastPosBin - firstPosBin + 1;
  
  // The number of bins in the covariance matrix.
  int covMatBins = nFidBins * nEneBins;

  // Create the Covariance Matrices needed.
  TMatrixF CovMatrixSig(covMatBins, covMatBins);
  TMatrixF CovMatrixBkg(covMatBins, covMatBins);
  TMatrixF CovMatrixStat(covMatBins, covMatBins);
  TMatrixF CovMatrixTot(covMatBins, covMatBins);
  
  // Create Output File
  TString outputName(inputName);
  outputName.ReplaceAll("Oscillation", "Minimization");
  
  TFile* outputFile = new TFile(outputName, "RECREATE");
  
  // Final Covariance matrices used in the analysis.
  TH2F* hSigCovMatrix = (TH2F*)ReducedSigMatrix->Clone("SigCovMatrix");
  hSigCovMatrix->Scale(0.0);
  TH2F* hBkgCovMatrix = (TH2F*)ReducedSigMatrix->Clone("BkgCovMatrix");
  hBkgCovMatrix->Scale(0.0);
  TH2F* hStatCovMatrix = (TH2F*)ReducedSigMatrix->Clone("StatCovMatrix");
  hStatCovMatrix->Scale(0.0);
  TH2F* hTotalCovMatrix = (TH2F*)ReducedSigMatrix->Clone("TotalCovMatrix");
  hTotalCovMatrix->Scale(0.0);
 
  // Create Chi Square Histogram
  TH2D* ChiSquareMap = new TH2D("ChiSquareMap","#chi^{2} Values; sin^{2} 2#theta_{14}; #Delta m^{2}_{14}", fNSinSq2Theta, fSinSq2ThetaBins, fNDeltam2-1, fDeltam2Bins);
  
  // Loop through the delta m values.
  for (int m = 0; m < fNDeltam2; m++) {
    std::cout << "Evaluating Delta m2 = " << fDeltam2[m] << std::endl;
    
    // Loop through sin2 2theta values.
    for (int s = 0; s < fNSinSq2Theta; s++) {
      
      // Build temporary histogram
      TH2F* OscHistTemp = (TH2F*)LvsENull->Clone("OscHistTemp");
      OscHistTemp->Add(LvsEDeltam2[m], -1.0 * fSinSq2Theta[s]);
  
      // Create the "X" vector
      TVectorF XVector(covMatBins);
      
      for (int i = 0; i < covMatBins; i++) {
        // Adjust the position bin to be the first position bin with entries
        int posIndexI = TMath::Floor((double) i / nEneBins);
        int eneIndexI = i - posIndexI * nEneBins;
        
        // Get the total statistics in bin I
        double predSigI = LvsENull->GetBinContent(eneIndexI+1, posIndexI+firstPosBin);
        double predBkgI = LvsEBackground->GetBinContent(eneIndexI+1, posIndexI+firstPosBin);
        double predTotI = predSigI + predBkgI;
        
        double obsSigI = OscHistTemp->GetBinContent(eneIndexI+1, posIndexI+firstPosBin);
        double obsTotI = obsSigI + predBkgI;
        
        // Fill X vector
        XVector(i) = obsSigI-predSigI;
        
        for (int j = 0; j < covMatBins; j++) {
          int posIndexJ = TMath::Floor((double) j / nEneBins);
          int eneIndexJ = j - posIndexJ * nEneBins;
          
          // Get the total statistics in bin J
          double predSigJ = LvsENull->GetBinContent(eneIndexJ+1, posIndexJ+firstPosBin);
          double predBkgJ = LvsEBackground->GetBinContent(eneIndexJ+1, posIndexJ+firstPosBin);
          double predTotJ = predSigJ + predBkgJ;
          
          double obsSigJ = OscHistTemp->GetBinContent(eneIndexJ+1, posIndexJ+firstPosBin);
          double obsTotJ = obsSigJ + predBkgJ;
          
          // Get Equivalent Covariance Matrix Bins
          // Due to skipping some position bins,
          // the matrix i, j will not match the histogram elements.
          double covI = (nEneBins) * (posIndexI+firstPosBin-1) + eneIndexI + 1;
          double covJ = (nEneBins) * (posIndexJ+firstPosBin-1) + eneIndexJ + 1;
          
          // Build up stats matrix
          // The error on S+B is Sqrt(S+B), if we assume that
          // the error on the background is known, the error
          // on just the signal is S/Sqrt(S+B).
          if (i == j) {
            // This error considers all error comes from the number of signal events
            //CovMatrixStat(i, j) = obsSigI * obsSigJ / TMath::Sqrt(obsTotI*obsTotJ);
            // This error is exactly the same as the chi2 approach
            CovMatrixStat(i, j) = TMath::Sqrt(predTotI*predTotJ);
            hStatCovMatrix->SetBinContent(covI, covJ, CovMatrixStat(i, j));
            
          }
                    
          // Build up signal systematic matrix
          // Need to convert the reduced covariance matrix back into a full covariance matrix
          //CovMatrixSig(i, j) = ReducedSigMatrix->GetBinContent(covI, covJ) * obsSigI * obsSigJ;
          CovMatrixSig(i, j) = TMath::Power(1.0, 2) * predSigI * predSigJ;
          
          hSigCovMatrix->SetBinContent(covI, covJ, CovMatrixSig(i, j));
          
        } // End of loop through j bins.
      } // End of loop through i bins.
      
      // Find Transpose
      TVectorF XTVector = XVector;
      
      CovMatrixTot = CovMatrixStat + CovMatrixSig;
      
      // Invert matrix
      CovMatrixTot.Invert();

      // Perform X^T * C^{-1} * X
      XTVector *= CovMatrixTot;
      double CovarianceChiSquare = XTVector * XVector;
      
      // Fill Chi Square Histogram
      ChiSquareMap->Fill(fSinSq2Theta[s], fDeltam2[m], CovarianceChiSquare);

    } // End of loop through sin22theta values.
  } // End of loop through delta m2 values.
  
  // Write covariance matrices used
  hStatCovMatrix->Write();
  hSigCovMatrix->Write();
  
  hTotalCovMatrix->Add(hStatCovMatrix);
  hTotalCovMatrix->Add(hSigCovMatrix);
  //hTotalCovMatrix->Add(hBkgCovMatrix);
  hTotalCovMatrix->Write();
  
  // Write Chi Square Histogram
  ChiSquareMap->Write();
  
}
