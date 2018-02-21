// C++ Includes
#include <iostream>
#include <algorithm>
#include <fstream>

// ROOT Includes
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

// Local Header Files Includes
#include "AnalysisInputs.H"
#include "DefinePositionTree.H"

////////////////////////////////
// Karin Gilje, October 9, 2015
//
// Create Oscillated sample:
//
/////////////////////////////////

void SimulateOscillation() {

  // Ensure a File is loaded
  if (!gFile) {
    std::cout << "File must be attached" << std::endl;
    return;
  }
  
  TFile* origFile = gFile;
  TString inputName(gFile->GetName());
  std::cout << "Accessing " << inputName << std::endl;
  
  // Load necessary histograms
  TH2F* BaselineToSegment = (TH2F*)gFile->Get("BaselineToSegment");
  TH2F* TrueEToRecoE = (TH2F*)gFile->Get("TrueEToRecoE");
  TH2F* LvsENull = (TH2F*)gFile->Get("LvsE");
  
  // Load background histograms to pass to the next stage
  TH2F* LvsEBackground = (TH2F*)gFile->Get("BkgL");
  LvsEBackground->SetName("LvsEBackground");
  
  // Load the Position Mapping
  ReadPositionTree();
  std::cout << "Position Tree: " << PositionTree->GetEntries() << std::endl;;

  // Create Output File
  TString outputName(inputName);
  outputName.ReplaceAll("Setup", "Oscillation");
  
  TFile* outputFile = new TFile(outputName, "RECREATE");
  
  // Copy the null case to new file
  LvsENull->Write();
  
  // Copy the background to new file
  LvsEBackground->Write();
  
  ConstructDeltam2();
  
  // Create the delta m2 list of histograms
  for (int i = 0; i < fNDeltam2; i++) {
    
    double deltam2 = fDeltam2[i];
    
    // Copy the L vs E structure
    TH2F* LvsETemp = (TH2F*)LvsENull->Clone(TString::Format("LvsE_%d", i));
    LvsETemp->Scale(0.0);
    
    // Loop through segments and get measured baselines
    for (int j = 0; j < BaselineToSegment->GetNbinsX(); j++) {
      // Make sure Delta m2 is nonzero
      if (deltam2 == 0) continue;
      
      PositionTree->GetEntry(j);
      double measuredBaseline = gBaseline;

      // Loop through true baselines
      for (int k = 0; k < BaselineToSegment->GetNbinsY(); k++) {
        double trueBaseline = BaselineToSegment->GetYaxis()->GetBinCenter(k+1);
        double baselineWeight = BaselineToSegment->GetBinContent(j+1, k+1);
        if (baselineWeight == 0) continue;

        // Loop through reconstructed energy
        for (int l = 0; l < TrueEToRecoE->GetNbinsX(); l++) {
          double measuredEnergy = TrueEToRecoE->GetXaxis()->GetBinCenter(l+1);

          // Loop through neutrino energies
          for (int m = 0; m < TrueEToRecoE->GetNbinsY(); m++) {
            double trueEnergy =  TrueEToRecoE->GetYaxis()->GetBinCenter(m+1);
            double energyWeight = TrueEToRecoE->GetBinContent(l+1, m+1);
            if (energyWeight == 0) continue;
            
            // Calculate the Delta m2 contribution
            double dm2Term = TMath::Power(TMath::Sin(1.27 * deltam2 * trueBaseline / trueEnergy), 2);

            // Fill the histogram
            LvsETemp->Fill(measuredEnergy, measuredBaseline, baselineWeight * energyWeight * dm2Term);
            
          } // End of loop through neutrino energies
        } // End of loop through measured energies
      } // End of loop through true baselines
    } // End of loop through measured baseline
    
    // Write the histograms
    LvsETemp->Write();
    
  } // End of loop through delta m2 values
  
}
