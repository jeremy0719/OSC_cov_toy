// C++ Includes
#include <iostream>
#include <algorithm>
#include <fstream>

// ROOT Includes
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>

// Local Header Files Includes
#include "AnalysisInputs.H"
#include "DefineEventTree.H"
#include "DefinePositionTree.H"

////////////////////////////////
// Karin Gilje, October 5, 2015
//
// Main Sensitivity Code:
//
// Detector: 0 - Front Near Detector Only
//           1 - 1.5m Back Near Detector Only
//           2 - 3m Back Near Detector Only
//           3 - Far Detector Only
//
// Exposure:  Time of Exposure in years.
/////////////////////////////////

void CreateExperiment(int DetectorType, double Exposure) {
  
  // Pass the detector exposure to memory.
  fExposure = Exposure;

  // Set the seed to pull from the computer time
  gRandom->SetSeed(0);

  // Create the Output File
  TString outputName = TString::Format("Setup_HFIR_%d_%.2f.root", DetectorType, Exposure);
  TFile* outputFile = new TFile(outputName.Data(), "RECREATE");
  
  std::cout << "Using HFIR Reactor" << std::endl;

  // Build the Tree of Events
  CreateEventTree();

  // Read in the Energy Spectrum Information
  ReadEnergySpectrum();
  
  // Set the global detector values
  SetDetectorVariables(DetectorType);
  
  // Fill the Tree of Events
  if (DetectorType < 4) {
      FillEventTree();
  }
  else {
    std::cout << "No viable option chosen." << std::endl;
    std::cout << "Options are: " << std::endl;
    std::cout << "\t 0- default near position" << std::endl;
    std::cout << "\t 1- 1.5 m back from default near position" << std::endl;
    std::cout << "\t 2- 3 m back from default near position" << std::endl;
    std::cout << "\t 3- default far position" << std::endl;
    std::cout << "Exiting ..." << std::endl;
    std::exit(1);
  }
  
  // Write the Event Tree
  WriteEventTree();
  
  // At this point all the events have been generated
  std::cout << "Events Generated" << std::endl;
  
  // Map the true baseline to individual segments
  std::cout << "Mapping the true baseline position to the center of the segment" << std::endl;
  int TotalSegments = fSegX*fSegZ;
  double minBaseline = TMath::Sqrt(TMath::Power(fDistX,2) + TMath::Power(fDistY,2) + TMath::Power(fDistZ,2))-fSegWidth;
  double maxBaseline = TMath::Sqrt(TMath::Power(fDistX + fSegX*fSegWidth,2)
                                   + TMath::Power(fDistY - fSegLength/2.0,2)
                                   + TMath::Power(fDistZ + fSegZ*fSegWidth,2))+fSegWidth;

  // Make the width of the bin the size of the segment
  double nBinsL = TMath::Ceil((maxBaseline - minBaseline)/fSegWidth);
  maxBaseline = minBaseline + nBinsL*fSegWidth;
  double posBinning[(int)nBinsL+1];
  
  // Make array of bins
  for (int i = 0; i < nBinsL+1; i++) {
    posBinning[i] = minBaseline + i * fSegWidth;
  }
  
  // Maps the true baseline to the segment
  TH2F* BaselineToSegment = new TH2F("BaselineToSegment",
                                     "Map of Segments to true Baselines; Segment Index; True Baseline (m)",
                                     TotalSegments, -0.5, TotalSegments-0.5, nBinsL*5, minBaseline, maxBaseline);
  
  // Maps the measured baseline bins to the segment
  TH2F* PositionToSegment = new TH2F("PositionToSegment",
                                     "Map of Segments to Positions; Segment Index; Position (m)",
                                     TotalSegments , -0.5, TotalSegments-0.5, nBinsL, minBaseline, maxBaseline);
  
  // Next step, use the generated events to figure out which ones are in the "fiducial" volume.
  TH3F* PseudoDetector = new TH3F("PseudoDetector", "Fake Detector; X (m); Y(m); Z(m)",
                                  fSegX, fDistX, fDistX + fSegX*fSegWidth,
                                  106, fDistY - fSegLength/2.0, fDistY + fSegLength/2.0,
                                  fSegZ, fDistZ, fDistZ + fSegZ*fSegWidth);
  
  // This tree will map a segment integer and map it to a measured baseline from the reactor.
  CreatePositionTree();
  FillPositionTree();
  WritePositionTree();
  
  // Fill the position to segment map
  // First get the total position weight
  for (int i = 0; i < fNSeg; i++) {
    PositionTree->GetEntry(i);

    // Use fiducialized volume
    if (gSegmentX == 0 || gSegmentX == fSegX-1) continue;
    if (gSegmentZ == 0 || gSegmentZ == fSegZ-1) continue;
  
    PositionToSegment->Fill(gSegment, gBaseline);
  }
  
  // Find the fiducial binning range.
  int minFidBin = PositionToSegment->FindFirstBinAbove(0., 2);
  int maxFidBin = PositionToSegment->FindLastBinAbove(0., 2) + 1;

  double minFiducial = PositionToSegment->GetYaxis()->GetBinLowEdge(minFidBin);
  int nBinsFid = maxFidBin - minFidBin;
  double fidBinning[nBinsFid+1];
  
  // Make array of Bins
  for (int i = 0; i < nBinsFid+1; i++) {
    fidBinning[i] = minFiducial + i * fSegWidth;
  }
  
  
  // Maps the measured baseline bins to the fiducial segment
  // This cuts out the position bins that refer to non-fiducial segments
  TH2F* PositionToSegmentFid = new TH2F("PositionToSegmentFid",
                                     "Map of Segments to Positions; Segment Index; Position (m)",
                                     TotalSegments , -0.5, TotalSegments-0.5, nBinsFid, fidBinning[0], fidBinning[nBinsFid]);
  
  // Fill the position to segment map
  // First get the total position weight
  for (int i = 0; i < fNSeg; i++) {
    PositionTree->GetEntry(i);
    
    // Use fiducialized volume
    if (gSegmentX == 0 || gSegmentX == fSegX-1) continue;
    if (gSegmentZ == 0 || gSegmentZ == fSegZ-1) continue;
    
    PositionToSegmentFid->Fill(gSegment, gBaseline);
  }
  
  // Save Histogram to File
  PositionToSegment->Write();
  PositionToSegmentFid->Write();
  
  // Fill the PseudoDetector Histogram with the positions of the generated events.
  for (int ev = 0; ev < fNEvents; ev++) {
    EventTree->GetEvent(ev);
    
    PseudoDetector->Fill(tDetectorPos->X(), tDetectorPos->Y(), tDetectorPos->Z(), tWeight);
    
    //Calculate Segment
    int binX = 0;
    int binY = 0;
    int binZ = 0;
    int globalBin = PseudoDetector->FindBin(tDetectorPos->X(), tDetectorPos->Y(), tDetectorPos->Z());
    PseudoDetector->GetBinXYZ(globalBin, binX, binY, binZ);
    int segment = (binX-1) + (binZ-1) * fSegX;
    
    // Perform Fiducial Cuts
    if (binX == 1 || binX == fSegX) continue;
    if (binZ == 1 || binZ == fSegZ) continue;
    if (TMath::Abs(tDetectorPos->Y()) > fSegLength/2.0-fSegWidth) continue;

    BaselineToSegment->Fill(segment, tTrueBaseline, tWeight);
  }

  PseudoDetector->Write();
  BaselineToSegment->Write();

  std::cout << "Mapping the true neutrino energy to the reconstructed energy" << std::endl;
  
  // Create a true neutrino energy to reconstructed energy map
  // Make the variable reconstructed binning.
  // This is based on Daya Bay binning with 0.2 MeV bins between
  // 1.0 and 8.0 MeV and one bin for events above 8.0 MeV.
  // This gives a total of 36 Energy Bins.
  int nRecoEBins = 37;
  double recoEBinning[nRecoEBins+1];
  
  for (int i = 0; i < nRecoEBins; i++) {
    recoEBinning[i] = 0.8 + i * 0.2;
  }
  recoEBinning[nRecoEBins] = 10.0;
  
  // Construct the true energy binning
  int nTrueEBins = 100;
  double trueEBinning[nTrueEBins+1];
  
  for (int i = 0; i < nTrueEBins+1; i++) {
    trueEBinning[i] = 0.0 + i * 0.1;
  }
  
  // Build the energy map histogram
  TH2F* TrueEToRecoEFine = new TH2F("TrueEToRecoEFine", "Map of Reconstructed Energy to True #nu Energy; Reconstructed Energy (MeV); True #nu Energy (MeV)",
                                    nTrueEBins, trueEBinning, nTrueEBins, trueEBinning);
  TH2F* TrueEToRecoE = new TH2F("TrueEToRecoE", "Map of Reconstructed Energy to True #nu Energy; Reconstructed Energy (MeV); True #nu Energy (MeV)",
                                nRecoEBins, recoEBinning, nTrueEBins, trueEBinning);
  
  // Create a temporary energy mapping to hold the smeared values for one energy
  TH2F* tempEnergyMap = new TH2F("tempEnergyMap", "Map of Reconstructed Energy to True #nu Energy; Reconstructed Energy (MeV); True #nu Energy (MeV)",
                                    nTrueEBins, trueEBinning, nTrueEBins, trueEBinning);
  // Reconstructed Energy smearing Function
  TF1* EnergyResolution = new TF1("EnergyResolution", "1/(TMath::Sqrt(2*TMath::Pi())*[1]) * TMath::Exp(-0.5*TMath::Power((x-[0])/[1], 2))", 0.0, 10.0);

  // Fill the fine binned histogram, we are doing small steps to more accurately portray the smearing.
  for (unsigned int i = 0; i < fAntiNuEnergy.size(); i++) {
    
    // Clear the temporary energy map
    tempEnergyMap->Scale(0.0);
    
    double trueEnergy = fAntiNuEnergy[i];

    // Skip events below event threshold
    if (trueEnergy < 1.8) continue;
    
    // Shift Energy by the difference between neutrino energy and prompt energy (0.8 MeV).
    // This comes from the 1.8 MeV threshold - 2*0.511 MeV of the annihillation of the positron.
    // This is the "reconstructed" energy
    double recoEnergy = fAntiNuEnergy[i] - 0.8;
    double eneWeight = fAntiNuFlux[i]*fIBDCrossSection[i];
    
    // Set smearing values
    EnergyResolution->SetParameter(0, recoEnergy);
    EnergyResolution->SetParameter(1, fEneRes / TMath::Sqrt(recoEnergy));
    
    // Calculate the integral of this bit of the map
    double tempIntegral = 0.0;
    
    // Loop through bins and smear
    for (int j = 0; j < nTrueEBins; j++) {
      double tempEnergy = TrueEToRecoEFine->GetXaxis()->GetBinCenter(j+1)-0.8;
      double reWeight = EnergyResolution->Eval(tempEnergy);
      if (reWeight == 0) continue;
      
      tempIntegral += eneWeight * reWeight;
      
      tempEnergyMap->Fill(tempEnergy, trueEnergy, eneWeight * reWeight);
    } // End of smearing
    
    // Ensure there was content in the temporary histogram
    if (tempIntegral != 0 && eneWeight != 0.0) TrueEToRecoEFine->Add(tempEnergyMap, eneWeight / tempIntegral);

  } // End of loop through all energies
  
  TrueEToRecoEFine->Write();
  
  // Rebin to variable binning
  for (int i = 0; i < nTrueEBins; i++) {
    for (int j = 0; j < nTrueEBins; j++) {
      double recoEnergy = TrueEToRecoEFine->GetXaxis()->GetBinCenter(i+1);
      double trueEnergy = TrueEToRecoEFine->GetYaxis()->GetBinCenter(j+1);
      double weight = TrueEToRecoEFine->GetBinContent(i+1, j+1);
      TrueEToRecoE->Fill(recoEnergy, trueEnergy, weight);
    }
  }
  
  TrueEToRecoE->Write();
  
  std::cout << "Making the null oscillation L-E histogram examples." << std::endl;
  
  // Create a segment vs E plot used in covariance matrix calculations.
  int nBinsSeg = fNSeg;
  double segBinning[nBinsSeg+1];
  for (int i = 0; i < nBinsSeg+1; i++) {
    segBinning[i] = -0.5 + i;
  }
  
  TH2F* SegvsE = new TH2F("SegvsE","Null Oscillation; Prompt Energy (MeV); Segment Index",
                          nRecoEBins, recoEBinning, nBinsSeg, segBinning);
  
  // Create the measured L vs E plot
  TH2F* LvsE = new TH2F("LvsE", "Null Oscillation; Prompt Energy (MeV); Measured Baseline(m)",
                        nRecoEBins, recoEBinning, nBinsFid, fidBinning);
  
  // First get the total position weight
  for (int i = 0; i < fNSeg; i++) {
    double posWeight = BaselineToSegment->Integral(i+1, i+1, 0, -1);
    PositionTree->GetEntry(i);
    double position = gBaseline;
    if (posWeight == 0) continue;
    // Now get the total Energy weight
    for (int j = 0; j < nRecoEBins; j++) {
      double energy = TrueEToRecoE->GetXaxis()->GetBinCenter(j+1);
      double eneWeight = TrueEToRecoE->Integral(j+1, j+1, 0, -1);

      double totWeight = posWeight*eneWeight;
      LvsE->Fill(energy, position, totWeight);
      SegvsE->Fill(energy, gSegment, totWeight);
    }
  }
  
  LvsE->Write();
  SegvsE->Write();
  
  std::cout << "Constructing the Shape of the background from txt input" << std::endl;

  // The background comes from Michael Mendenhall and includes detector response.
  
  // Pull Information from file (comes from Michael Mendenhall)
  // The background is in mHz/MeV and normalized over the fiducialized near detector.
  std::ifstream ifs ("SimulatedBackground.txt");
  
  double elow;
  double ehi;
  double bkg;
  double bkgerr;
  double sig;
  double sigerr;
  
  std::vector< std::vector<double> > BackgroundValues; // Input is mHz/MeV
  
  if (ifs.is_open()) {
    std::string line;
    while (getline(ifs, line)) {
      switch (line[0]) { // check first character, skip all lines that start with #
        case '#' :
          break;
        default :
          ifs >> elow >> ehi >> bkg >> bkgerr >> sig >> sigerr;
          std::vector<double> temp;
          temp.clear();
          temp.push_back(elow);
          temp.push_back(ehi);
          temp.push_back(bkg);
          BackgroundValues.push_back(temp);
          break;
      } // end of switch
    } // end of while loop
  }
  else {
    std::cout << "Background File not available!" << std::endl;
    std::cout << "Returning... " << std::endl;
    return;
  }

  // This describes the cosmic background PER SEGMENT! in the fiducial volume.
  TH1F* Bkg1D = new TH1F("Bkg1D", "Background; Energy (MeV); Entries (mHz/MeV)",
                                nRecoEBins, recoEBinning);
  for (unsigned int i = 0; i < BackgroundValues.size(); i++) {
    double pullE = (BackgroundValues[i][0]+BackgroundValues[i][1])/2.0;
    double eWidth = BackgroundValues[i][1]-BackgroundValues[i][0];
    double bkgWeight = BackgroundValues[i][2] *
    fYearToSeconds/1000.0 * // extra factor converts mHz into 1/yr
    Exposure * // Include the amount of time exposed.
    eWidth / // Include the bin width in MeV
    80.0; // Include the average over the number of fiducial segments.
    // The background was calculated with fiducialized 12x10 segments or 10x8 fiducial segments.
    Bkg1D->Fill(pullE, bkgWeight);
  }

  Bkg1D->Write();
  
  // Create the background histogram for fiducialized volume.
  TH2F* BkgSeg = new TH2F("BkgSeg",
                          "Simulated Background; Energy (MeV); Segment Number",
                          nRecoEBins, recoEBinning, nBinsSeg, segBinning);
  TH2F* BkgL = new TH2F("BkgL",
                        "Simulated Background; Energy (MeV); Baseline (m)",
                        nRecoEBins, recoEBinning, nBinsFid, fidBinning);
  
  // Fill the background histogram for each segment.
  for (int i = 0; i < fNSeg; i++) {
    // Get the measured position
    PositionTree->GetEntry(i);
    double pullL = gBaseline;
    double pullS = gSegment;
    
    // Use fiducialized volume
    if (gSegmentX == 0 || gSegmentX == fSegX-1) continue;
    if (gSegmentZ == 0 || gSegmentZ == fSegZ-1) continue;
    
    for (int i = 0; i < Bkg1D->GetNbinsX(); i++) {
      double pullE = Bkg1D->GetBinCenter(i+1);
      double bkgWeight = Bkg1D->GetBinContent(i+1);
      BkgL->Fill(pullE, pullL, bkgWeight);
      BkgSeg->Fill(pullE, pullS, bkgWeight);
    } // End loop through Bkg bins
  }
  
  BkgSeg->Write();
  BkgL->Write();
  
  // Check the signal to background
  
  // Work between 0.8 MeV and 7.2 MeV
  double lowEBin = 1;
  double highEBin = 32;
  double lowPBin = LvsE->FindFirstBinAbove(10,2);
  double highPBin = LvsE->FindLastBinAbove(10,2);
  
  // Calculate the normalization of each histogram.
  double nEventsNull = LvsE->Integral(lowEBin, highEBin, lowPBin, highPBin);
  // In mHz/MeV with 0.2 MeV bins.
  double nEventsBkg = BkgL->Integral(lowEBin, highEBin, lowPBin, highPBin);
  
  // Ensure far detector has same S:B as near.
  if (DetectorType == 3) {
    double SignalToBackground = 3.25;
    BkgL->Scale(nEventsNull / SignalToBackground / nEventsBkg);
    BkgSeg->Scale(nEventsNull / SignalToBackground / nEventsBkg);
    nEventsBkg = BkgL->Integral(lowEBin, highEBin, lowPBin, highPBin);
  }
  
  std::cout << "Writing Background Histogram." << std::endl;
  std::cout << " Min E: " << BkgL->GetXaxis()->GetBinLowEdge(lowEBin) << std::endl;
  std::cout << " Max E: " << BkgL->GetXaxis()->GetBinLowEdge(highEBin) + BkgL->GetXaxis()->GetBinWidth(highEBin) << std::endl;
  std::cout << " Number of Signal: " << nEventsNull << std::endl;
  std::cout << " Number of Background: " << nEventsBkg << std::endl;
  std::cout << " double check S/B : " << nEventsNull / nEventsBkg << std::endl;
  
}
