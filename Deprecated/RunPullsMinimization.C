#include <iostream>

#include <TROOT.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>

#include "chiSquare_simplified.H"

void RunPullsMinimization(double RxOnFraction = 0.41) {
  
  std::cout << "Running assuming " << RxOnFraction*100.0 << "% Live Time of the Reactor." << std::endl;
  
  // Ensure a File is loaded
  if (!gFile) {
    std::cout << "File must be attached" << std::endl;
    return;
  }
  
  TFile* origFile = gFile;
  TString inputName(gFile->GetName());
  std::cout << "Accessing " << inputName << std::endl;
  

  // Fill the null oscillation histogram
  TH2D* NullOscillation = (TH2D*)gFile->Get("LvsE");

  // Already has the right binning!
  TH2F* Background = (TH2F*)gFile->Get("LvsEBackground");

  // Load necessary histograms
  std::vector<TH2F*> LvsEDeltam2;
  LvsEDeltam2.clear();
  
  for (int i = 0; i < 57; i++) {
    TH2F* LvsETemp = (TH2F*)gFile->Get(TString::Format("LvsE_%d", i));
    
    // Rescale histograms with RxOnFraction
    LvsETemp->Scale(RxOnFraction);
    LvsEDeltam2.push_back(LvsETemp);
  }
  
  // Create Output File Name
  TString outputName(inputName);
  outputName.ReplaceAll("Oscillation", "MinimizedChi");
  
  TFile* outputFile = new TFile(outputName, "RECREATE");
  
  // Get the chi-Square Fitting Object
  ChiSquareDayaBay chiSquare;

  // Setup array of Delta m^2 values to test.  The first value is 0.0 which refers to the unoscillated situation.
  int nDeltam2 = 57;
  double deltam2[nDeltam2]; // eV^2
  deltam2[0] = 0.0; // Null Oscillation
  for (int i = 1; i < nDeltam2; i++) {
    deltam2[i] = TMath::Power(10, -1.4+0.05*(i-1));
  }
  
  // Setup array of Sin2 2Theta values to test.
  int nSinSq2Theta = 48;
  double sinSq2Theta[nSinSq2Theta];
  for(int i = 0; i < nSinSq2Theta; i++) {
    sinSq2Theta[i] = TMath::Power(10, -2.4 + 0.05*i);
  }

  // Set Bounds of the Minimization.
  chiSquare.SetBounds(NullOscillation);

  // Add Null Oscillation to chiSquare object
  NullOscillation->Scale(RxOnFraction);
  chiSquare.SetSignalTestOscillation(NullOscillation);
  chiSquare.SetSignalNoOscillation(NullOscillation);
  
  // Add Background to chiSquare object
  Background->Scale(RxOnFraction);
  chiSquare.SetBackground(Background);
  
  chiSquare.Nevents(NullOscillation);
  chiSquare.Nevents(Background);
  
  // Create a ChiSquare map
  // First make bins
  double m2Bins[nDeltam2];
  for (int i = 0; i < nDeltam2; i++) {
    m2Bins[i] = TMath::Power(10, -1.425 + 0.05*i);
  }
  double sinBins[nSinSq2Theta+1];
  for (int i = 0; i < nSinSq2Theta+1; i++) {
    sinBins[i] = TMath::Power(10, -2.425 + i*0.05);
  }
  
  TH2D* ChiSquareMap = new TH2D("ChiSquareMap","#chi^{2} Values; sin^{2} 2#theta; #Delta m^{2}", nSinSq2Theta, 0.004, 1.0, nDeltam2-1, 0.04, 25.0);
  (ChiSquareMap->GetYaxis())->Set(nDeltam2-1, m2Bins);
  (ChiSquareMap->GetXaxis())->Set(nSinSq2Theta, sinBins);
  
  chiSquare.SetSigmaEbin();
  chiSquare.SetSigmaSource(1.0);
  
  for (int i = 1; i < nDeltam2; i++) {
    std::cout << "Working on Delta m2 " << deltam2[i] << " at position " << i << std::endl;
    chiSquare.SetSignalMaxOscillation(LvsEDeltam2[i]);

    for (int j = 0; j < nSinSq2Theta; j++) {
      chiSquare.SetSinSq2Theta(sinSq2Theta[j]);
      chiSquare.FindMinimum();
      double chiSq = chiSquare.Minimum();
      //double chiSq = chiSquare.EvalStat();
      
      ChiSquareMap->Fill(sinSq2Theta[j], deltam2[i], chiSq);

    } // End loop through Sin^2 2\theta Values
  } // End loop through Delta m^2 values
  
  ChiSquareMap->Write();
  
  std::cout << "Successfully reached end!" << std::endl;
  
}
