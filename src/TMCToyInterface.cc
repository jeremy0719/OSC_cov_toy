#include "TMCToyInterface.hh"

// Default Class Constructor
TMCToyInterface::TMCToyInterface(TExperiment &experiment):fExperiment(experiment),fCovMatInt(experiment,true) {
  
  // Make vectors from the experiment class
  v_ReactorOnSignal = ExtractVector(fExperiment.hLvsERef);
  v_ReactorOnBackground = ExtractVector(fExperiment.hBkgRxOnLvsE);
  v_ReactorOffBackground = ExtractVector(fExperiment.hBkgRxOffLvsE);
  
  std::cout << "Number of Vector Elements: " << v_ReactorOnSignal->GetNrows() << std::endl;
  
  // Make Matrix from covariance histograms
  // Grab histograms from file
  int eneBins;
  for (auto it:fExperiment.fDetectors) {
    int detId = it->fDetectorCode;
    
    TString histName;
    histName.Form("LvsE_%i",detId);
    TH2D* hLvsENull = (TH2D*)(fExperiment.hLvsENull.at(detId)->Clone(histName));
    // Calculate the number of energy bins (the same for every position!)
    eneBins = hLvsENull->GetNbinsX();
    double posBins = hLvsENull->GetNbinsY();
    fCovMatInt.fDetCovMatBins[detId] = eneBins * posBins;
    fCovMatInt.fDetCovMatBins[detId] = eneBins * posBins;
    fCovMatInt.fCumulativeCovMatBins += fCovMatInt.fDetCovMatBins.at(detId);
//    SetupCovarianceMatrices(detId);
  }
  fCovMatInt.SetupCovarianceMatrices(eneBins);
  fCovMatInt.BuildCovarianceMatrices();
  
  std::cout << "Matrices retrieved." << std::endl;
  
  // Build Matrices into histogram format
//  BuildCovarianceMatrices();
  
  std::cout << "Matrices built and inflated." << std::endl;
  
  // Ensure matrix and vectors have the same number of elements
  // Vectors should be 1xN and Matrices should be NxN
  int NumberOfMatrixRows = (fCovMatInt.fReactorOnSigCovMatrix).GetNrows();
  int NumberOfVectorEntries = v_ReactorOnSignal->GetNrows();
  
  if (NumberOfMatrixRows != NumberOfVectorEntries) {
    std::cout << "Something went horribly wrong!" << std::endl;
    std::cout << "Number of Matrix Rows: " << NumberOfMatrixRows << std::endl;
    std::cout << "Number of Vector Rows: " << NumberOfVectorEntries << std::endl;
    std::cout << "Returning..." << std::endl;
    exit(1);
  }
  
  // Create Toy throwing instances
  toyReactorOnSignal = new TThrowMCToy(*v_ReactorOnSignal, (fCovMatInt.fReactorOnSigCovMatrix));
  toyReactorOnBackground = new TThrowMCToy(*v_ReactorOnBackground, (fCovMatInt.fReactorOnBkgCovMatrix));
  toyReactorOffBackground = new TThrowMCToy(*v_ReactorOffBackground, (fCovMatInt.fReactorOffBkgCovMatrix));
  
  return;
}

// Class Deconstructor
TMCToyInterface::~TMCToyInterface() {
  
  // Delete vectors
  if (v_ReactorOnSignal != NULL) delete v_ReactorOnSignal;
  if (v_ReactorOnBackground != NULL) delete v_ReactorOnBackground;
  if (v_ReactorOffBackground != NULL) delete v_ReactorOffBackground;
  
//  // Delete matrices
//  if (fCovMatInt.fReactorOnSigCovMatrix != NULL) delete fCovMatInt.fReactorOnSigCovMatrix;
//  if (fCovMatInt.fReactorOnBkgCovMatrix != NULL) delete fCovMatInt.fReactorOnBkgCovMatrix;
//  if (fCovMatInt.fReactorOffBkgCovMatrix != NULL) delete fCovMatInt.fReactorOffBkgCovMatrix;
  
  // Delete MC Toy instances
  if (toyReactorOnSignal != NULL) delete toyReactorOnSignal;
  if (toyReactorOnBackground != NULL) delete toyReactorOnBackground;
  if (toyReactorOffBackground != NULL) delete toyReactorOffBackground;
  
  return;
}

// Set the random seed
void TMCToyInterface::SetSeed(int seed = 1867) {
  toyReactorOnSignal->SetSeed(seed);
  toyReactorOnBackground->SetSeed(seed);
  toyReactorOffBackground->SetSeed(seed);
  return;
}

// Turn map of TH2Ds to a vector
TVectorD* TMCToyInterface::ExtractVector(std::map<int, TH2D*> inputMap) {
  
  int nBinsTot = 0;
  
  // Loop through the map to get the size of the output vector
  std::map<int, TH2D*>::iterator it;
  for (it = inputMap.begin(); it != inputMap.end(); it++) {
    
    // Calculate the number of bins.
    int nBinsX = (it->second)->GetNbinsX();
    int nBinsY = (it->second)->GetNbinsY();
    
    nBinsTot += nBinsX * nBinsY;
    
  } // End of loop through detectors
  // Create Output vector
  TVectorD* outputVector = new TVectorD(nBinsTot);
  
  // Index of position in output vector
  int iElement = 0;
  
  // Loop through the map again, filling the vector as you go.
  for (it = inputMap.begin(); it != inputMap.end(); it++) {
    
    // Calculate the number of bins.
    int nBinsX = (it->second)->GetNbinsX();
    int nBinsY = (it->second)->GetNbinsY();
    
    for (int j = 0; j < nBinsY; j++) {
      for (int i = 0; i < nBinsX; i++) {
        (*outputVector)(iElement) = (it->second)->GetBinContent(i+1, j+1);
        iElement++;
      } // End of loop through Y bins
    } // End of loop through X bins
    
  } // End of loop through detectors
  
  return outputVector;
}

// Get one complete toy
void TMCToyInterface::GetToy() {
  toyReactorOnSignal->ThrowExperiment(v_RxOnSignalToy);
  toyReactorOnBackground->ThrowExperiment(v_RxOnBackgroundToy);
  toyReactorOffBackground->ThrowExperiment(v_RxOffBackgroundToy);
  
  return;
}

// Fill a new map of histograms with toy information
std::map<int, TH2D*> TMCToyInterface::FillOutputHistogram(std::map<int, TH2D*> inputMap, TVectorD inputVector, TString RxStatus) {
  
  // Initialize output map
  std::map<int, TH2D*> outputMap;
  
  // Initialize the index of the vector
  int iIndex = 0;
  
  // Loop through the map
  std::map<int, TH2D*>::iterator it;
  for (it = inputMap.begin(); it != inputMap.end(); it++) {
    // Access detector id needed for naming of histograms
    int detId = it->first;
    
    // Clone the histogram to save in a new map
    TH2D* hTemp = (TH2D*)(it->second)->Clone(TString::Format("LvsERx%s_%d", RxStatus.Data(), detId));
    
    // Clear histogram values
    hTemp->Scale(0.);
    
    // Calculate the number of bins.
    int nBinsX = hTemp->GetNbinsX();
    int nBinsY = hTemp->GetNbinsY();
    
    for (int j = 0; j < nBinsY; j++) {
      for (int i = 0; i < nBinsX; i++) {
        
        // Check that the vector is big enough
        if (iIndex > inputVector.GetNoElements()) {
          std::cout << "Something went horribly wrong!" << std::endl;
          std::cout << iIndex <<  " expected number of elements: " << inputVector.GetNoElements() << std::endl;
          std::cout << "Exiting ..." << std::endl;
          exit(1);
        }
        
        hTemp->SetBinContent(i+1, j+1, inputVector(iIndex));
        iIndex++;
      } // End of loop through Y bins
    } // End of loop through X bins
    
    // Add elements to map
    outputMap[detId] = hTemp;
    
  } // End of loop through detectors
  return outputMap;
}

// Turn the toy into data-like information
void TMCToyInterface::CreateFakeData() {
  
  // Make maps of detector id to TH2Ds
  // Reverse the extract vector process
  
  TVectorD v_ToyReactorOn(v_RxOnSignalToy.size());
  TVectorD v_ToyReactorOff(v_RxOnSignalToy.size());
  
  for (unsigned int i = 0; i < v_RxOnSignalToy.size(); i++) {
    v_ToyReactorOn(i) = v_RxOnSignalToy[i] + v_RxOnBackgroundToy[i];
    v_ToyReactorOff(i) = v_RxOffBackgroundToy[i];
  }
  
  hReactorOn = FillOutputHistogram(fExperiment.hLvsERef, v_ToyReactorOn, "On");
  hReactorOff = FillOutputHistogram(fExperiment.hLvsERef, v_ToyReactorOff, "Off");
  
  return;
}

// The result of this function should be maps of detector TH2Ds,
// one for signal events and one for background events.
void TMCToyInterface::ThrowToy() {
  GetToy();
  CreateFakeData();
  
  return;
}

void TMCToyInterface::WriteHistograms(TFile& outputFile) {
  
  outputFile.cd();
  
  std::map<int, TH2D*>::iterator it;
  
  for (it = hReactorOn.begin(); it != hReactorOn.end(); it++) {
    it->second->Write();
  }
  
  for (it = hReactorOff.begin(); it != hReactorOff.end(); it++) {
    it->second->Write();
  }
  
  fCovMatInt.WriteHistograms(outputFile);
  
  return;
}
