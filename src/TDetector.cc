#include "TDetector.hh"

int TDetector::DetectorCount;

void TDetector::ReadMacroInputs()
{
  TMacroInterface &PROSPECTMacroInterface = TMacroInterface::Instance();
  TString key;
  // Set the number of segments from macro
  key.Form("XSegments%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fNSegX);
  key.Form("ZSegments%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fNSegZ);
  
  // Set segment width and length from macro
  key.Form("SegmentWidth%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fSegWidth);
  key.Form("SegmentLength%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fSegLength);
  
  // Set efficiency, position and energy resolutions from macro
  key.Form("Efficiency%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fDetectorEfficiency);
  key.Form("YPositionResolution%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fPosResY);
  key.Form("EnergyResolution%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fEnergyResolution);
  
  // Set outer segments that will be taken out to make the fiducial volume from macro
  key.Form("XOuterSegmentsLeft%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fNOuterXLeft);
  key.Form("XOuterSegmentsRight%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fNOuterXRight);
  key.Form("ZOuterSegmentsTop%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fNOuterZTop);
  key.Form("ZOuterSegmentsBottom%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fNOuterZBottom);
  
  // Set position resolutions of the segment along x and z from macro
  // Only applies if we have better position resoultion than the size of the segment
  key.Form("XPositionResolution%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fPosResX);
  key.Form("ZPositionResolution%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fPosResZ);
  
  // Set the number of baseline bins to be used for the detector
  key.Form("NBinsL%i",fDetLocation.fDetectorType);
  PROSPECTMacroInterface.RetrieveValue(key,fNBinsL);
}

void TDetector::ComputeDetectorParameters(int detType, int detPosition)
{
  fDetectorCode = (100*detType)+detPosition;
  fDetLocation = TDetectorLocation(detType,detPosition);
  if(fDetLocation.fDetectorType==1){
    fDetectorName = "Near";
    
    // Default is 14x11
    fNSegX = 14;
    fNSegZ = 11;
    ReadMacroInputs();
    fNSeg = fNSegX*fNSegZ;
    
    fNFiducialSegments = (fNSegX-fNOuterXLeft-fNOuterXRight)* (fNSegZ-fNOuterZTop-fNOuterZBottom);
  }
  else if(fDetLocation.fDetectorType==2){
    fDetectorName = "Far";
    
    fNSegX = 28;
    fNSegZ = 20;
    ReadMacroInputs();
    fNSeg = fNSegX*fNSegZ;
    
    fNFiducialSegments = (fNSegX-fNOuterXLeft-fNOuterXRight)* (fNSegZ-fNOuterZTop-fNOuterZBottom);
  }
  else{
    printf("Detector number wrong!\n");
    exit(1);
  }

  fMinBaseline = TMath::Sqrt(TMath::Power(fDetLocation.fDistX - fNSegX * fSegWidth/2.0,2) + TMath::Power(fDetLocation.fDistY - fSegLength/2.0,2) + TMath::Power(fDetLocation.fDistZ - fNSegZ*fSegWidth/2.0,2));
  fMaxBaseline = TMath::Sqrt(TMath::Power(fDetLocation.fDistX + fNSegX*fSegWidth/2.0,2)
                             + TMath::Power(fDetLocation.fDistY - fSegLength/2.0,2)
                             + TMath::Power(fDetLocation.fDistZ + fNSegZ*fSegWidth/2.0,2));
  if(fNBinsL==0) fNBinsL = TMath::Ceil((fMaxBaseline - fMinBaseline)/fSegWidth);
  
  // Set the number of baseline bins using the width of the segment
  printf("Using the %s detector at position %i\n",fDetectorName.Data(),fDetLocation.fDetectorPosition);
  TDetector::DetectorCount++;
}

TDetector::TDetector(int detType, int detPosition)
{
  ComputeDetectorParameters(detType,detPosition);
}

void TDetector::DefineDetector(int detType, int detPosition)
{
  if(fDetectorCode!=0){
    printf("Cannot redefine a detector!\n");
    exit(1);
  }
  ComputeDetectorParameters(detType,detPosition);
}

void TDetector::CheckIfDetectorDefined()
{
  if(fDetectorCode==0){
    printf("TDetector: Detector not defined!\n");
    exit(1);
  }
}

void TDetector::LocalToGlobal(double& globalX, double& globalY, double& globalZ, double localX, double localY, double localZ)
{
  CheckIfDetectorDefined();
  // Rotate detector around Z-axis (leaves z unaffected)
  double detectorXnew = localX * TMath::Cos(fDetLocation.fAngleOffset * TMath::DegToRad())
  - localY * TMath::Sin(fDetLocation.fAngleOffset * TMath::DegToRad());
  double detectorYnew = localY * TMath::Cos(fDetLocation.fAngleOffset * TMath::DegToRad())
  + localX * TMath::Sin(fDetLocation.fAngleOffset * TMath::DegToRad());
  
  // Add in shift from reactor
  globalX = detectorXnew + fDetLocation.fDistX;
  globalY = detectorYnew + fDetLocation.fDistY;
  globalZ = localZ + fDetLocation.fDistZ;
}

void TDetector::LocalToGlobal()
{
  CheckIfDetectorDefined();
  LocalToGlobal(eventGlobalX, eventGlobalY, eventGlobalZ,eventLocalX, eventLocalY, eventLocalZ);
}



void TDetector::CreatePositionTree(int& sIndex,int& sIndexX,int& sIndexZ,double& bline)
{
  CheckIfDetectorDefined();
  TString positionTreeName ;
  positionTreeName.Form("PositionTree%i",fDetectorCode);
  fPositionTree = new TTree(positionTreeName.Data(),"Flat Tree of the Segment Baselines");
  fPositionTree->Branch("Segment", &sIndex, "Segment/I");
  fPositionTree->Branch("SegmentX", &sIndexX, "SegmentX/I");
  fPositionTree->Branch("SegmentZ", &sIndexZ, "SegmentZ/I");
  fPositionTree->Branch("Baseline", &bline, "Baseline/D");
  printf("Position tree for detector %i created!\n",fDetectorCode);
  isPositionTreeCreated=true;
}

void TDetector::CreatePositionTree()
{
  CheckIfDetectorDefined();
  CreatePositionTree(segmentIndex, segmentIndexX, segmentIndexZ, baseline);
}

void TDetector::ResetPositionTreeValues()
{
  CheckIfDetectorDefined();
  segmentIndex = 0;
  segmentIndexX = 0;
  segmentIndexZ = 0;
  baseline = 0.0;
}

void TDetector::FillPositionTree()
{
  CheckIfDetectorDefined();
  if(!isPositionTreeCreated){
    printf("Create position trees before filling them!\n");
    exit(1);
  }
  // Fill the map of segment to positon
  for (int j = 0; j < fNSegZ; j++) {
    for (int i = 0; i < fNSegX; i++) {
      ResetPositionTreeValues();
      
      segmentIndexX = i;
      segmentIndexZ = j;
      
      segmentIndex = i + j * fNSegX;
      // Find the local position in detector coordinates of the center of the segment
      double localX = fSegWidth*(i - (fNSegX-1)/2.0);
      double localY = 0.;
      double localZ = fSegWidth*(j - (fNSegZ-1)/2.0);
      
      // Shift to global coordinates
      double globalX = 0.;
      double globalY = 0.;
      double globalZ = 0.;
      
      LocalToGlobal(globalX, globalY, globalZ, localX, localY, localZ);
      
      baseline = TMath::Sqrt(TMath::Power(globalX,2) + TMath::Power(globalY,2) + TMath::Power(globalZ,2));
      fPositionTree->Fill();
    }
  } // End of loop through segments
  printf("Position tree for the detector %i filled!\n",fDetectorCode);
} // End of Fill Tree function

void TDetector::ReadPositionTree(int& sIndex,int& sIndexX,int& sIndexZ,double& bline)
{
  CheckIfDetectorDefined();
  if(!isPositionTreeCreated){
    printf("Create position trees before reading the position tree!\n");
    exit(1);
  }
  // Set the branch addresses
  fPositionTree->SetBranchAddress("Segment", &sIndex);
  fPositionTree->SetBranchAddress("SegmentX", &sIndexX);
  fPositionTree->SetBranchAddress("SegmentZ", &sIndexZ);
  fPositionTree->SetBranchAddress("Baseline", &bline);
}

void TDetector::ReadPositionTree(TString inputFileName)
{
  CheckIfDetectorDefined();
  TFile *inputFile = TFile::Open(inputFileName.Data());
  if(!inputFile->IsOpen() || inputFile->IsZombie()){
    printf("The file is not open or is a zombie\n");
    exit(1);
  }
  
  TString positionTreeName ;
  positionTreeName.Form("PositionTree%i",fDetectorCode);
  fPositionTree = (TTree*)inputFile->Get(positionTreeName.Data());
  
  // Check that the tree exists
  if (!fPositionTree) {
    printf("Incomplete file, tree doesn't exist\n");
    exit(1);
  }
  ReadPositionTree(segmentIndex, segmentIndexX, segmentIndexZ, baseline);
}

void TDetector::WritePositionTree(TFile& outputFile)
{
  CheckIfDetectorDefined();
  if(!isPositionTreeCreated){
    printf("Create position trees before writing them!\n");
    exit(1);
  }
  outputFile.cd();
  // Important to set CloneTree, or else only the headers will be copied not the data
  fPositionTree->Write();
  printf("Position tree for %i written to file %s !\n",fDetectorCode,outputFile.GetName());
}

void TDetector::CreateDetectorEventTree()
{
  printf("Creating the new tree!\n");
  CheckIfDetectorDefined();
  TString detectorEventTreeName ;
  detectorEventTreeName.Form("DetectorEventTree%i",fDetectorCode);
  fDetectorEventTree = new TTree(detectorEventTreeName.Data(),"Flat Tree of the Detector Events");
  fDetectorEventTree->Branch("Segment", &eventSegmentIndex, "Segment/I");
  fDetectorEventTree->Branch("SegmentX", &eventSegmentIndexX, "SegmentX/I");
  fDetectorEventTree->Branch("SegmentZ", &eventSegmentIndexZ, "SegmentZ/I");
  fDetectorEventTree->Branch("Fiducial", &isEventFiducial, "Fiducial/O");
  fDetectorEventTree->Branch("LocalX", &eventLocalX, "LocalX/D");
  fDetectorEventTree->Branch("LocalY", &eventLocalY, "LocalY/D");
  fDetectorEventTree->Branch("LocalZ", &eventLocalZ, "LocalZ/D");
  fDetectorEventTree->Branch("GlobalX", &eventGlobalX, "GlobalX/D");
  fDetectorEventTree->Branch("GlobalY", &eventGlobalY, "GlobalY/D");
  fDetectorEventTree->Branch("GlobalZ", &eventGlobalZ, "GlobalZ/D");
  printf("Detector event tree for detector %i created!\n",fDetectorCode);
  isDetectorEventTreeCreated=true;
}

void TDetector::ResetDetectorEventTree()
{
  CheckIfDetectorDefined();
  eventSegmentIndex = 0;
  eventSegmentIndexX = 0;
  eventSegmentIndexZ = 0;
  isEventFiducial = true;
  eventLocalX = 0.;
  eventLocalY = 0.;
  eventLocalZ = 0.;
  eventGlobalX = 0.;
  eventGlobalY = 0.;
  eventGlobalZ = 0.;
}

void TDetector::FillDetectorEventTree()
{
  
  CheckIfDetectorDefined();
  if(!isDetectorEventTreeCreated){
    printf("Create Detector Event Trees before filling them!\n");
    exit(1);
  }
  
  TRandom3 random(0);
  
  for (int i = 0; i < fNEvents; i++){
    ResetDetectorEventTree();
    // Create the detector position (in reference to center, local coordinates)
    eventLocalX = random.Uniform(-fNSegX / 2.0 * fSegWidth, fNSegX / 2.0 * fSegWidth);
    eventLocalY = random.Uniform(-fSegLength/2.0, fSegLength/2.0);
    eventLocalZ = random.Uniform(-fNSegZ / 2.0 * fSegWidth, fNSegZ / 2.0 * fSegWidth);
    
    // Figure out Segment for event position
    eventSegmentIndexX = TMath::Floor((eventLocalX + fNSegX*fSegWidth / 2.0) / fSegWidth);
    eventSegmentIndexZ = TMath::Floor((eventLocalZ + fNSegZ*fSegWidth / 2.0) / fSegWidth);    
    eventSegmentIndex = eventSegmentIndexX + eventSegmentIndexZ * fNSegX;
    
    // Check if the event is in the fiducial volume
    if (eventSegmentIndexX == 0 || eventSegmentIndexX == fNSegX-1) isEventFiducial = false;
    if (eventSegmentIndexZ == 0 || eventSegmentIndexZ == fNSegZ-1) isEventFiducial = false;
    if (TMath::Abs(eventLocalY) > fSegLength/2.0 - fSegWidth) isEventFiducial = false;
    
    // Change local detector coordinates to global coordinates
    LocalToGlobal();
    
    // Fill the tree
    fDetectorEventTree->Fill();
  }
  
  printf("Detector event tree for the detector %i filled!\n",fDetectorCode);

}

void TDetector::ReadDetectorEventTree(int& sIndex, int& sIndexX, int& sIndexZ, bool& isFid)
{
  CheckIfDetectorDefined();
  if(!isDetectorEventTreeCreated){
    printf("Create detector event tree before reading!\n");
    exit(1);
  }
  // Set the branch addresses
  fDetectorEventTree->SetBranchAddress("Segment", &sIndex);
  fDetectorEventTree->SetBranchAddress("SegmentX", &sIndexX);
  fDetectorEventTree->SetBranchAddress("SegmentZ", &sIndexZ);
  fDetectorEventTree->SetBranchAddress("Fiducial", &isFid);
}

void TDetector::ReadDetectorEventTree(double& posX, double& posY, double& posZ)
{
  CheckIfDetectorDefined();
  if(!isDetectorEventTreeCreated){
    printf("Create detector event tree before reading!\n");
    exit(1);
  }
  // Set the branch addresses
  fDetectorEventTree->SetBranchAddress("GlobalX", &posX);
  fDetectorEventTree->SetBranchAddress("GlobalY", &posY);
  fDetectorEventTree->SetBranchAddress("GlobalZ", &posZ);
}

void TDetector::ReadDetectorEventTree(TString inputFileName)
{
  CheckIfDetectorDefined();
  TFile *inputFile = TFile::Open(inputFileName.Data());
  if(!inputFile->IsOpen() || inputFile->IsZombie()){
    printf("The file is not open or is a zombie\n");
    exit(1);
  }
  
  TString detectorEventTreeName ;
  detectorEventTreeName.Form("DetectorEventTree%i",fDetectorCode);
  fDetectorEventTree = (TTree*)inputFile->Get(detectorEventTreeName.Data());
  
  // Check that the tree exists
  if (!fDetectorEventTree) {
    printf("Incomplete file, tree doesn't exist\n");
    exit(1);
  }
  
  // Set the branch addresses
  fDetectorEventTree->SetBranchAddress("Segment", &eventSegmentIndex);
  fDetectorEventTree->SetBranchAddress("SegmentX", &eventSegmentIndexX);
  fDetectorEventTree->SetBranchAddress("SegmentZ", &eventSegmentIndexZ);
  fDetectorEventTree->SetBranchAddress("Fiducial", &isEventFiducial);
  fDetectorEventTree->SetBranchAddress("LocalX", &eventLocalX);
  fDetectorEventTree->SetBranchAddress("LocalY", &eventLocalY);
  fDetectorEventTree->SetBranchAddress("LocalZ", &eventLocalZ);
  fDetectorEventTree->SetBranchAddress("GlobalX", &eventGlobalX);
  fDetectorEventTree->SetBranchAddress("GlobalY", &eventGlobalY);
  fDetectorEventTree->SetBranchAddress("GlobalZ", &eventGlobalZ);
}

void TDetector::WriteDetectorEventTree(TFile& outputFile)
{
  CheckIfDetectorDefined();
  if(!isDetectorEventTreeCreated){
    printf("Create detector event trees before writing them!\n");
    exit(1);
  }
  outputFile.cd();
  fDetectorEventTree->Write();
  printf("Detector Event tree for %i written to file %s !\n",fDetectorCode,outputFile.GetName());
}

void TDetector::Print()
{
  CheckIfDetectorDefined();
  fDetLocation.Print();
  printf("PROSPECT Detector : Type is %i\n",fDetLocation.fDetectorType);
  printf("PROSPECT Detector : Name is %s\n",fDetectorName.Data());
  printf("PROSPECT Detector : Position identifier is %i\n",fDetLocation.fDetectorPosition);
  printf("PROSPECT Detector : Segment width is %f\n",fSegWidth);
  printf("PROSPECT Detector : Efficiency is %2.2f\n",fDetectorEfficiency*100.0);
  printf("PROSPECT Detector : Position resolution along X, Y and Z directions are %f, %f and %f respectively\n",fPosResX,fPosResY,fPosResZ);
  printf("PROSPECT Detector : Energy resolution is %f\n",fEnergyResolution);
  printf("PROSPECT Scintillator : Proton density of the liquid scinitllator is %f\n",fProtonDensity);
  printf("PROSPECT Detector : Number of segments in x and z directions are %i and %i respectively\n",fNSegX,fNSegZ);
  printf("PROSPECT Detector : Total number of segments is %i\n",fNSeg);
  printf("PROSPECT Detector : Number of fiducial segments is %i\n",fNFiducialSegments);
  printf("PROSPECT Detector : Distance to first segment in x, y and z are %f, %f and %f respectively\n",fDetLocation.fDistX,fDetLocation.fDistY,fDetLocation.fDistZ);
}


