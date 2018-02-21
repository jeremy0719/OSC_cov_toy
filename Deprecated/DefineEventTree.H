#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"

//////////////////////////////////////////////////////////////
// Karin Gilje, October 1, 2015
//
// Definition of tree and tree variables to access the
// information about the important baselines in the oscillation
// calculation.
//////////////////////////////////////////////////////////////

////////////////////////////////////////
// Define variables to be saved in tree.
////////////////////////////////////////

TTree* EventTree;

// Raw reactor coordinates
TVector3* tReactorPos = new TVector3;

// Raw detector coordinates
TVector3* tDetectorPos = new TVector3;

// Distance from Reactor-Pixel to Detector-Pixel
// Also known as the true baseline for the oscillation calculation.
Float_t tTrueBaseline = 0.0;

// Weight of event from baseline
Float_t tWeight = 0.0;

// Function sets the branches of the tree
void CreateEventTree() {
  EventTree = new TTree("EventTree","Flat Tree of the Interaction Baselines");
  
  EventTree->Branch("ReactorPos", "TVector3", &tReactorPos);
  EventTree->Branch("DetectorPos", "TVector3", &tDetectorPos);
  EventTree->Branch("TrueBaseline", &tTrueBaseline, "TrueBaseline/F");
  EventTree->Branch("Weight", &tWeight, "Weight");
}

// This function resets all the flat tree variables to default ones
void ResetEventTreeValues() {
  tReactorPos->SetXYZ(0., 0.,-5.);
  tDetectorPos->SetXYZ(0., 0., -5.);
  tTrueBaseline = -1.0;
  tWeight = 0.0;
}

// This function fills the tree
void FillEventTree() {
  
  // Find the total detector volume.
  double detVolume = fSegX*fSegWidth * fSegLength * fSegZ*fSegWidth;

  // Calculate normalization: applies to all events regardless
  // of energy or position
  // This, combined with the position weight, the energy weight
  // and the exposure, gives a weight for the event.
  double normalization = fProtonDensity * fDetectorEfficiency * fReactorPower * fYearToSeconds * fExposure;

  // Create all events
  for (int i = 0; i < fNEvents; i++) {
    ResetEventTreeValues();
    
    // Create the reactor position
    double reactorX = gRandom->Uniform(-fReactorRadius, fReactorRadius);
    double reactorY = gRandom->Uniform(-fReactorRadius, fReactorRadius);
    double reactorZ = gRandom->Uniform(-fReactorHeight/2.0, fReactorHeight/2.0);
  
    // Ensure position is within a cylinder
    while (reactorX*reactorX + reactorY*reactorY > fReactorRadius*fReactorRadius) {
      reactorX = gRandom->Uniform(-fReactorRadius, fReactorRadius);
      reactorY = gRandom->Uniform(-fReactorRadius, fReactorRadius);
    }

    tReactorPos->SetXYZ(reactorX, reactorY, reactorZ);

    // Create the detector position
    double detectorX = gRandom->Uniform(fDistX, fDistX + fSegX*fSegWidth);
    double detectorY = gRandom->Uniform(fDistY - fSegLength/2.0, fDistY + fSegLength/2.0);
    double detectorZ = gRandom->Uniform(fDistZ, fDistZ + fSegZ*fSegWidth);

    tDetectorPos->SetXYZ(detectorX, detectorY, detectorZ);
    
    // Calculate the True Baseline
    TVector3 baseline(tDetectorPos->X() - tReactorPos->X(),
                      tDetectorPos->Y() - tReactorPos->Y(),
                      tDetectorPos->Z() - tReactorPos->Z());
    tTrueBaseline = baseline.Mag();
    
    // Calculate the position weight
    double posWeight = detVolume / fNEvents * 1.0 / (4 * TMath::Pi() * tTrueBaseline * tTrueBaseline);
    
    // This weight is ONLY missing the weighting effect from the energy and the Exposure time.
    // This weight is equivalent to one year of running.,.
    tWeight = normalization * posWeight;
    
    EventTree->Fill();
    
  } // End of Creating Event Position
} // End of Fill Tree function

// This function gives access to reading the event tree
void ReadEventTree() {
  EventTree = (TTree*)gFile->Get("EventTree");
  
  // Check that the tree exists
  if (!EventTree) {
    std::cout << "Incomplete file, tree doesn't exist..." << std::endl;
    exit(1);
  }
  
  // Set the branch addresses
  EventTree->SetBranchAddress("ReactorPos", &tReactorPos);
  EventTree->SetBranchAddress("DetectorPos", &tDetectorPos);
  EventTree->SetBranchAddress("TrueBaseline", &tTrueBaseline);
  EventTree->SetBranchAddress("Weight", &tWeight);
  
}

void WriteEventTree() {
  EventTree->Write();
}