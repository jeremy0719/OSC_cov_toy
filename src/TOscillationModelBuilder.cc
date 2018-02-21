#include "TOscillationModelBuilder.hh"
#include "TCanvas.h"

TOscillationModelBuilder::TOscillationModelBuilder(TExperiment & experiment):fExperiment(experiment)
{
  // Load necessary histograms and position tree
  for(auto it:(fExperiment.fDetectors)){
    int detNumber = it->fDetectorCode;
        hSimulatedNueEnergy[detNumber]=(TH1D*)(fExperiment.hSimulatedNueEnergy.at(detNumber)->Clone());
    //    hLvsENull[detNumber] = (TH2D*)(fExperiment.hLvsENull.at(detNumber)->Clone());
  }
}

TOscillationModelBuilder::~TOscillationModelBuilder()
{
  //Since the arrays are dynamic, it leads to memory leaks if not deleted
  delete [] fDeltam2Bins;
  delete [] fSinSq2ThetaBins;
}

int TOscillationModelBuilder::EstimateBinNumbers(double deltam2Value, double sinSq2ThetaValue, int& deltam2Bin,int& sin22theta2Bin)
{
  deltam2Bin = (std::log10(deltam2Value)+1.425)*fNDeltam2/(std::log10(28.2)+1.425);
  sin22theta2Bin = (std::log10(sinSq2ThetaValue)+2.425)*fNSinSq2Theta/(std::log10(2.425)+2.4);
  return 1;
}

void TOscillationModelBuilder::ConstructDeltam2Bins()
{
  fDeltam2Bins = new double[fNDeltam2+1];
  for (int i = 1; i <= fNDeltam2; i++) {
    fDeltam2Bins[i-1] = fDeltam2.at(i-1)*TMath::Power(10, -0.025);
  }
  fDeltam2Bins[fNDeltam2] = fDeltam2.at(fNDeltam2-1)*TMath::Power(10, 0.025);
}

void TOscillationModelBuilder::ConstructDeltam2()
{
  fDeltam2.clear();
  fDeltam2.push_back(0.0);
  double deltam2StepSize = (std::log10(28.2) +1.4)/fNDeltam2;
  for (int i = 1; i <= fNDeltam2; i++) {
    fDeltam2.push_back(TMath::Power(10, -1.4+ deltam2StepSize*(i-1)));
  }
  ConstructDeltam2Bins();
}

void TOscillationModelBuilder::ConstructSinSq2ThetaBins()
{
  fSinSq2ThetaBins = new double[fNSinSq2Theta+1];
  for (int i = 1; i <= fNSinSq2Theta; i++) {
    fSinSq2ThetaBins[i-1] = fSinSq2Theta.at(i-1)*TMath::Power(10, -0.025);
  }
  fSinSq2ThetaBins[fNSinSq2Theta] = fSinSq2Theta.at(fNSinSq2Theta-1)*TMath::Power(10, 0.025);
}

void TOscillationModelBuilder::ConstructSinSq2Theta()
{
  fSinSq2Theta.clear();
  fSinSq2Theta.push_back(0.0);
  double sinSq2ThetaStepSize = 2.45/fNSinSq2Theta;
  for (int i = 1; i <= fNSinSq2Theta; i++) {
    fSinSq2Theta.push_back(TMath::Power(10, -2.4+ sinSq2ThetaStepSize*(i)));
  }
  ConstructSinSq2ThetaBins();
}

void TOscillationModelBuilder::ConstructBinsFromFile(TString input)
{
  printf("Using file %s for deltam2 and sin22theta binning \n",input.Data());
  if(input.EndsWith(".txt",TString::kIgnoreCase))
  {
    // read the input file
    std::ifstream inputFile(input.Data());
    //Check if the file is open
    if (!inputFile.is_open()) {
      printf("Error opening the file\n");
      exit(1);
    }
    
    bool fillDeltam2=false;
    
    std::string line;
    while(std::getline(inputFile, line)){
      TString inputLine(line);
      inputLine.ReplaceAll(" ","");//remove any spaces
      if (inputLine.IsWhitespace()) continue;//skip if the line is empty
      if (inputLine.BeginsWith("#")) continue;// skip if the line starts with #
      //Since the first set of bins have to be deltam2, set fillDeltam2 to true so that deltam2 bins will be filled
      if (inputLine.Contains("deltam2",TString::kIgnoreCase)){
        fDeltam2.clear();
        fillDeltam2=true;
      }
      // When encountered with SinSq2Theta fillDeltam2 is set to false to fill sin22theta bins
      else if (inputLine.Contains("SinSq2Theta",TString::kIgnoreCase)){
        fSinSq2Theta.clear();
        fillDeltam2=false;
      }
      else{
        if (fillDeltam2==true){
          fDeltam2.push_back(inputLine.Atof());
        }
        else{
          fSinSq2Theta.push_back(inputLine.Atof());
        }
      }
    }
    // Change number of bins based on the inputs
    fNDeltam2=fDeltam2.size();
    fNSinSq2Theta=fSinSq2Theta.size();
  }
  
  if(input.EndsWith(".root",TString::kIgnoreCase))
  {
    // read the input file
    TFile *inputFile = TFile::Open(input.Data());
    //Check if the file is open
    if(!inputFile->IsOpen() || inputFile->IsZombie()){
      printf("Error opening the file\n");
      exit(1);
    }
    // Obtain bin numbers from TH2D in the root file
    TH2D* binningFile = (TH2D*)inputFile->Get("OscBins");
    fNSinSq2Theta=binningFile->GetNbinsX();
    fNDeltam2=binningFile->GetNbinsY();
    
    for (int i = 1; i <= fNDeltam2; i++) {
      double binCenter=binningFile->GetYaxis()->GetBinCenter(i);
      fDeltam2.push_back(binCenter);
    }
    for (int i = 1; i <= fNSinSq2Theta; i++) {
      double binCenter=binningFile->GetXaxis()->GetBinCenter(i);
      fSinSq2Theta.push_back(binCenter);
    }
  }
  // Make sure to call these methods to fill Deltam2 and Sin22theta bins
  ConstructDeltam2Bins();
  ConstructSinSq2ThetaBins();
}

void TOscillationModelBuilder::SetUpHistograms()
{
  for(auto it:(fExperiment.fDetectors)){
    int detNumber = it->fDetectorCode;
    for (int i = 0; i < fNDeltam2; i++) {
      int deltam2Number = 1000*i + detNumber;
      hLvsEDeltam2[deltam2Number] = ((TH2D*)fExperiment.hLvsENull.at(detNumber)->Clone());
      hLvsEDeltam2Relative[deltam2Number] = ((TH2D*)fExperiment.hLvsENull.at(detNumber)->Clone());
      hLvsEDeltam2.at(deltam2Number)->Scale(0.0);
      hLvsEDeltam2Relative.at(deltam2Number)->Scale(0.0);
    }
  }
}

void TOscillationModelBuilder::SetUpModelOscillation(int ndeltam2, int nSinSq2Theta)
{
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  // Set bins for deltam2 ans sin22theta
  TString value;
  PROSPECTMacroInterface.RetrieveValue("InputOscBins",value);
  
  if(value.CompareTo("YES",TString::kIgnoreCase)==0){
    PROSPECTMacroInterface.RetrieveValue("OscBinsInputFile",value);
    ConstructBinsFromFile(value);
  }
  
  else{
    fNDeltam2=ndeltam2;
    fNSinSq2Theta=nSinSq2Theta;
    // Set bins for deltam2 ans sin22theta
    PROSPECTMacroInterface.RetrieveValue("NDeltam2",fNDeltam2);
    PROSPECTMacroInterface.RetrieveValue("NSinSq2Theta",fNSinSq2Theta);
    
    printf("Setting up oscillation parameters with %i bins in δm2 and %i bins in sin22θ \n",fNDeltam2,fNSinSq2Theta);
    ConstructDeltam2();
    ConstructSinSq2Theta();
  }
  SetUpHistograms();
}

void TOscillationModelBuilder::SimulateModelOscillation(const TExperiment& experiment,TDetector& detector,TH2D &hOsc, double deltam2)
{
  int detNumber = detector.fDetectorCode;
  hOsc.Scale(0.0);
  
  
  int segmentIndex = 0;   // The index of the segment
  int segmentIndexX = 0;  // The x index of the segment
  int segmentIndexZ = 0;  // The z index of the segment
  double measuredBaseline = 0.0;  // The "measured" baseline to the segment
  detector.ReadPositionTree(segmentIndex, segmentIndexX, segmentIndexZ, measuredBaseline);
  
  
  TH2D* hSegvsEOsc = (TH2D*)fExperiment.hSegvsETrue.at(detNumber)->Clone("hSegvsETrue");
  TH2D* hSegvsEOutOsc = (TH2D*)fExperiment.hSegvsENull.at(detNumber)->Clone("hSegvsEOsc");
  hSegvsEOsc->Clear();
  hSegvsEOutOsc->Clear();
  hSegvsEOsc->Scale(0);
  if (deltam2 == 0) return;
  
  detector.ResetPositionTreeValues();
  for (int i = 0; i < detector.fPositionTree->GetEntries(); i++) {
    // Get position weight of the tree
    double posWeight =  fExperiment.hBaselineToSegment.at(detNumber)->Integral(i+1, i+1, 0, -1);
    detector.fPositionTree->GetEntry(i);

    if (posWeight == 0) continue;
    // Now get the total reconstructed energy weight
    for (int j = 0; j < hSimulatedNueEnergy.at(detNumber)->GetNbinsX(); j++) {
      double trueEnergy= hSimulatedNueEnergy.at(detNumber)->GetBinCenter(j+1);
      double eneWeight= hSimulatedNueEnergy.at(detNumber)->GetBinContent(j+1);
      
      if (posWeight == 0) continue;
      double totWeight = posWeight*eneWeight;

      // Calculate the Delta m2 contribution
      double dm2Term = fOscillationSimulator.GetDeltam2Term(deltam2,trueEnergy,measuredBaseline);
      hSegvsEOsc->Fill(trueEnergy, segmentIndex, totWeight*dm2Term);
    }
  }
  
  // Apply detector response to the segvsE hist
  TDetectorResponseInterface fDetRespInt(detNumber);
  fDetRespInt.ApplyDetectorResponse(*hSegvsEOsc,hSegvsEOutOsc);
  
  detector.ResetPositionTreeValues();
  for (int i = 0; i < hSegvsEOutOsc->GetNbinsX(); i++) {
    double energy = hSegvsEOutOsc->GetXaxis()->GetBinCenter(i+1);
    // Loop through segments and get measured baselines
    for (int j = 0; j < detector.fPositionTree->GetEntries(); j++) {
      detector.fPositionTree->GetEntry(j);
      // Fill LvsE hists from SegvsE hist
      hOsc.Fill(energy, measuredBaseline, hSegvsEOutOsc->GetBinContent(i+1,j+1));
    }
  }
}

/*
 void TOscillationModelBuilder::SimulateModelOscillation(const TExperiment& experiment,TDetector& detector,TH2D* hRef, double deltam2)
 {
 int detNumber = detector.fDetectorCode;
 hRef->Clear();
 
 int segmentIndex = 0;   // The index of the segment
 int segmentIndexX = 0;  // The x index of the segment
 int segmentIndexZ = 0;  // The z index of the segment
 double measuredBaseline = 0.0;  // The "measured" baseline to the segment
 detector.ReadPositionTree(segmentIndex, segmentIndexX, segmentIndexZ, measuredBaseline);
 
 // Loop through segments and get measured baselines
 for (int j = 0; j < detector.fPositionTree->GetEntries(); j++) {
 // Make sure Delta m2 is nonzero
 if (deltam2 == 0) continue;
 detector.fPositionTree->GetEntry(j);
 // Loop through true baselines
 for (int k = 0; k < hBaselineToSegment.at(detNumber)->GetNbinsY(); k++) {
 double trueBaseline = hBaselineToSegment.at(detNumber)->GetYaxis()->GetBinCenter(k+1);
 double baselineWeight = hBaselineToSegment.at(detNumber)->GetBinContent(j+1, k+1);
 if (baselineWeight == 0) continue;
 
 // Loop through reconstructed energy
 for (int l = 0; l < hTrueEToRecoE.at(detNumber)->GetNbinsY(); l++) {
 double measuredEnergy = hTrueEToRecoE.at(detNumber)->GetYaxis()->GetBinCenter(l+1);
 // Loop through neutrino energies
 
 for (int m = 0; m < hTrueEToRecoE.at(detNumber)->GetNbinsX(); m++) {
 double trueEnergy =  hTrueEToRecoE.at(detNumber)->GetXaxis()->GetBinCenter(m+1);
 double energyWeight = hTrueEToRecoE.at(detNumber)->GetBinContent(m+1,l+1);
 
 if (energyWeight == 0) continue;
 // Calculate the Delta m2 contribution
 double dm2Term = fOscillationSimulator.GetDeltam2Term(deltam2,trueEnergy,trueBaseline);
 // Fill the histogram
 hRef->Fill(measuredEnergy, measuredBaseline,baselineWeight * energyWeight * dm2Term);
 } // End of loop through neutrino energies
 } // End of loop through measured energies
 } // End of loop through true baselines
 } // End of loop through measured baseline
 TDetectorResponse fDetResp;
 fDetResp.ApplyResponse(hRef,0.045);
 }*/

void TOscillationModelBuilder::SimulateModelOscillation()
{
  for(auto it:(fExperiment.fDetectors)){
    int detNumber = it->fDetectorCode;
    for (int i = 0; i < fNDeltam2; i++) {
      int deltam2Number = 1000*i + detNumber;
      // Copy the L vs E structure
      double deltam2 = fDeltam2[i];
      SimulateModelOscillation(fExperiment,*it,*hLvsEDeltam2.at(deltam2Number),deltam2);
      }
  }
}

void TOscillationModelBuilder::WriteHistograms(TFile& outputFile)
{
  printf("Writing oscillation histograms \n");
  if(!outputFile.IsOpen() || outputFile.IsZombie()){
    printf("File not open");
    return;
  }
  outputFile.cd();
  TString histName;
  for(auto it:(fExperiment.fDetectors)){
    int detNumber = it->fDetectorCode;
//    fExperiment.hBaselineToSegment.at(detNumber)->Write();
    //    fExperiment.hTrueEToRecoEFine.at(detNumber)->Write();
    //    fExperiment.hTrueEToRecoE.at(detNumber)->Write();
    fExperiment.hLvsENull.at(detNumber)->Write();
    for (int i = 0; i < fNDeltam2; i++) {
      int deltam2Number = 1000*i + detNumber;
      histName.Form("LvsEDeltam2_%i_%i",detNumber,i);
      hLvsEDeltam2.at(deltam2Number)->SetName(histName.Data());
      hLvsEDeltam2.at(deltam2Number)->Write();
    }
  }
  Printf("TOscillationModelBuilder files written to the file %s\n",outputFile.GetName());
}
