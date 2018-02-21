#include "TTree.h"

//////////////////////////////////////////////////////////////
// Karin Gilje, October 9, 2015
//
// Definition of tree and tree variables to access the
// information about the position of the segments.
//////////////////////////////////////////////////////////////

////////////////////////////////////////
// Define variables to be saved in tree.
////////////////////////////////////////

TTree* PositionTree;

// The global index of the segment
Int_t gSegment = 0;

// The x index of the segment
Int_t gSegmentX = 0;

// The z index of the segment
Int_t gSegmentZ = 0;

// The "measured" baseline to the segment
Float_t gBaseline = 0.0;

// Function sets the branches of the tree
void CreatePositionTree() {
  PositionTree = new TTree("PositionTree","Flat Tree of the Segment Baselines");
  PositionTree->Branch("Segment", &gSegment, "Segment/I");
  PositionTree->Branch("SegmentX", &gSegmentX, "SegmentX/I");
  PositionTree->Branch("SegmentZ", &gSegmentZ, "SegmentZ/I");
  PositionTree->Branch("Baseline", &gBaseline, "Baseline/F");
}

// This function resets all the flat tree variables to default ones
void ResetPositionTreeValues() {
  gSegment = 0;
  gSegmentX = 0;
  gSegmentZ = 0;
  gBaseline = 0.0;
}

// This function fills the tree
void FillPositionTree() {
  
  // Fill the map of segment to positon
  for (int j = 0; j < fSegZ; j++) {
    for (int i = 0; i < fSegX; i++) {
      ResetPositionTreeValues();
      
      gSegmentX = i;
      gSegmentZ = j;
      
      gSegment = i + j * fSegX;
      gBaseline = TMath::Sqrt(TMath::Power(fDistX + i*fSegWidth,2)
                              + TMath::Power(fDistZ + j*fSegWidth,2));
      
      PositionTree->Fill();
    }
  } // End of loop through segments
} // End of Fill Tree function

// This function gives access to reading the position tree
void ReadPositionTree() {
  PositionTree = (TTree*)gFile->Get("PositionTree");
  
  // Check that the tree exists
  if (!PositionTree) {
    std::cout << "Incomplete file, tree doesn't exist..." << std::endl;
    exit(1);
  }
  
  // Set the branch addresses
  PositionTree->SetBranchAddress("Segment", &gSegment);
  PositionTree->SetBranchAddress("SegmentX", &gSegmentX);
  PositionTree->SetBranchAddress("SegmentZ", &gSegmentZ);
  PositionTree->SetBranchAddress("Baseline", &gBaseline);
}

void WritePositionTree() {
  PositionTree->Write();
}