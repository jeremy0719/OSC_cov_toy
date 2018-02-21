#include "TReactor.hh"

void TReactor::ReadMacroInputs()
{
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  // Set the reactor related values from the macro
  PROSPECTMacroInterface.RetrieveValue("ReactorRadius",fReactorRadius);
  PROSPECTMacroInterface.RetrieveValue("ReactorHeight",fReactorHeight);
  PROSPECTMacroInterface.RetrieveValue("ReactorPower",fReactorPower);
}

TReactor::TReactor()
{
  ReadMacroInputs();
}

TVector3 TReactor::GetRandReactorPosition()
{
  /// A random position of the reactor. For a cylindrical assumption of the shape of the reactor, this is the random point in the cylinder.
  TVector3 randReactorPos;
  TRandom3 random(0);
  // Create the reactor position
  double reactorX = random.Uniform(-fReactorRadius, fReactorRadius);
  double reactorY = random.Uniform(-fReactorRadius, fReactorRadius);
  double reactorZ = random.Uniform(-fReactorHeight/2.0, fReactorHeight/2.0);
  
  // Ensure position is within a cylinder
  while (reactorX*reactorX + reactorY*reactorY > fReactorRadius*fReactorRadius) {
    reactorX = random.Uniform(-fReactorRadius, fReactorRadius);
    reactorY = random.Uniform(-fReactorRadius, fReactorRadius);
  }
  randReactorPos.SetXYZ(reactorX, reactorY, reactorZ);
  return randReactorPos;
}

void TReactor::CreateReactorEventTree()
{
  TString reactorEventTreeName ;
  reactorEventTreeName.Form("ReactorEventTree");
  fReactorEventTree = new TTree(reactorEventTreeName.Data(),"Flat Tree of the Reactor Events");
  fReactorEventTree->Branch("GlobalX", &reactionGlobalX, "GlobalX/D");
  fReactorEventTree->Branch("GlobalY", &reactionGlobalY, "GlobalY/D");
  fReactorEventTree->Branch("GlobalZ", &reactionGlobalZ, "GlobalZ/D");
  printf("Reactor event tree created!\n");
  isReactorEventTreeCreated=true;
}

void TReactor::ResetReactorEventTree()
{
  reactionGlobalX = 0.;
  reactionGlobalY = 0.;
  reactionGlobalZ = 0.;
}

void TReactor::FillReactorEventTree()
{
  
  if(!isReactorEventTreeCreated){
    printf("Create Reactor Event Tree before filling them!\n");
    exit(1);
  }
  
  TRandom3 random(0);
  
  for (int i = 0; i < fNEvents; i++){
    ResetReactorEventTree();
    // Create the reactor position
    reactionGlobalX = random.Uniform(-fReactorRadius, fReactorRadius);
    reactionGlobalY = random.Uniform(-fReactorRadius, fReactorRadius);
    reactionGlobalZ = random.Uniform(-fReactorHeight/2.0, fReactorHeight/2.0);
    
    // Ensure position is within a cylinder
    while (reactionGlobalX*reactionGlobalX + reactionGlobalY*reactionGlobalY > fReactorRadius*fReactorRadius) {
      reactionGlobalX = random.Uniform(-fReactorRadius, fReactorRadius);
      reactionGlobalY = random.Uniform(-fReactorRadius, fReactorRadius);
    }
    
    // Fill the tree
    fReactorEventTree->Fill();
  }

  printf("Reactor event tree filled!\n");
  
}

void TReactor::ReadReactorEventTree(double& posX, double& posY, double& posZ)
{
  if(!isReactorEventTreeCreated){
    printf("Create reactor event tree before reading!\n");
    exit(1);
  }
  // Set the branch addresses
  fReactorEventTree->SetBranchAddress("GlobalX", &posX);
  fReactorEventTree->SetBranchAddress("GlobalY", &posY);
  fReactorEventTree->SetBranchAddress("GlobalZ", &posZ);
}

void TReactor::ReadReactorEventTree(TString inputFileName)
{
  TFile *inputFile = TFile::Open(inputFileName.Data());
  if(!inputFile->IsOpen() || inputFile->IsZombie()){
    printf("The file is not open or is a zombie\n");
    exit(1);
  }
  
  TString reactorEventTreeName ;
  reactorEventTreeName.Form("ReactorEventTree");
  fReactorEventTree = (TTree*)inputFile->Get(reactorEventTreeName.Data());
  
  // Check that the tree exists
  if (!fReactorEventTree) {
    printf("Incomplete file, tree doesn't exist\n");
    exit(1);
  }
  
  // Set the branch addresses
  fReactorEventTree->SetBranchAddress("GlobalX", &reactionGlobalX);
  fReactorEventTree->SetBranchAddress("GlobalY", &reactionGlobalY);
  fReactorEventTree->SetBranchAddress("GlobalZ", &reactionGlobalZ);
}

void TReactor::WriteReactorEventTree(TFile& outputFile)
{
  if(!isReactorEventTreeCreated){
    printf("Create reactor event trees before writing them!\n");
    exit(1);
  }
  outputFile.cd();
  // Important to set CloneTree, or else only the headers will be copied not the data
  fReactorEventTree->Write();
  printf("Reactor Event tree written to file %s !\n",outputFile.GetName());
}

void TReactor::Print()
{
  printf("Reactor: Radius (m) of core %f\n",fReactorRadius);
  printf("Reactor: Heigth (m) of core %f\n",fReactorHeight);
  printf("Reactor: Power of core (GeV) %f\n",fReactorPower);
}
