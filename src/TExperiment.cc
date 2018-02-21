#include "TExperiment.hh"

void TExperiment::ReadMacroInputs(int detectorCode)
{
  // We mainly setup the exposure time here
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  // Set events from the macros
  PROSPECTMacroInterface.RetrieveValue("Events",fNEvents);
  TString key;
  
  
  //  if(isReferenceData){
  //    // TODO:These numbers are supposed to come from the data table
  //    fRxOnExposure[detectorCode]=2134080;
  //    fRxOffExposure[detectorCode]=2134080;
  //  }
  //
  //  else{
  double eRxOnSeconds=0;
  double eRxOnDays=0;
  double eRxOnCycles=0;
  
  // Set the number of seconds from the macros
  key.Form("RxOnSeconds%i",detectorCode);
  PROSPECTMacroInterface.RetrieveValue(key,eRxOnSeconds);
  // Set the number of days from the macros
  key.Form("RxOnDays%i",detectorCode);
  PROSPECTMacroInterface.RetrieveValue(key,eRxOnDays);
  // Set the number of cycles from the macros
  key.Form("RxOnCycles%i",detectorCode);
  PROSPECTMacroInterface.RetrieveValue(key,eRxOnCycles);
  
  // The precedence of assignment is seconds then days then cycles
  // i.e., if seconds are defined in the macro, the RxOn time is defined in seconds(and neglects everything else), otherwise days otherwise Cycles
  if(eRxOnSeconds!=0) printf("Using %f Seconds for reactor-on time\n",eRxOnSeconds );
  else if(eRxOnDays!=0){
    printf("Using %f days for defining reactor-on time\n",eRxOnDays);
    eRxOnSeconds=eRxOnDays*fDayToSeconds;
  }
  else{
    if(eRxOnCycles==0) eRxOnCycles=7;
    printf("Using %f cycles for defining reactor-on time\n",eRxOnCycles);
    eRxOnSeconds=eRxOnCycles*fCycleToSeconds;
  }
  fRxOnExposure[detectorCode]=eRxOnSeconds;
  
  double eRxOffSeconds=0;
  double eRxOffDays=0;
  double eRxOffCycles=0;
  
  // Set the number of seconds from the macros
  key.Form("RxOffSeconds%i",detectorCode);
  PROSPECTMacroInterface.RetrieveValue(key,eRxOffSeconds);
  // Set the number of days from the macros
  key.Form("RxOffDays%i",detectorCode);
  PROSPECTMacroInterface.RetrieveValue(key,eRxOffDays);
  // Set the number of cycles from the macros
  key.Form("RxOffCycles%i",detectorCode);
  PROSPECTMacroInterface.RetrieveValue(key,eRxOffCycles);
  
  // The precedence of assignment is seconds then days then cycles
  // i.e., if seconds are defined in the macro, the RxOff time is defined in seconds(and neglects everything else), otherwise days otherwise Cycles
  if(eRxOffSeconds!=0) printf("Using %f Seconds for reactor-off time \n",eRxOffSeconds );
  else if(eRxOffDays!=0){
    printf("Using %f days for defining reactor-off time\n",eRxOffDays);
    eRxOffSeconds=eRxOffDays*fDayToSeconds;
  }
  else{
    if(eRxOffCycles==0) eRxOffCycles=7;
    printf("Using %f cycles for defining reactor-off time\n",eRxOffCycles);
    eRxOffSeconds=eRxOffCycles*fCycleToSeconds;
  }
  fRxOffExposure[detectorCode]=eRxOffSeconds;
  //  }
  // Calculate the ratio of reactor off to on
  fRxOnOffRatio[detectorCode] =fRxOnExposure[detectorCode]/fRxOffExposure[detectorCode];
}

TExperiment::TExperiment()
{
  fReactor.fNEvents = fNEvents;
  isExperimentDefined=true;  // EDIT: Change this to include at least a detector
}

TExperiment::TExperiment(double deltam2, double sin22theta):fReferenceDeltam2(deltam2),fReferenceSin22theta(sin22theta)
{
  fReactor.fNEvents = fNEvents;
  isExperimentDefined=true;  // EDIT: Change this to include at least a detector
}

TExperiment::~TExperiment()
{
  for(auto it:fDetectors){
    printf("Deleting the detector %i\n",it->fDetectorCode);
    delete it;
  }
}

int TExperiment::GetNFidBins(const TDetector &detector)
{
  int detNumber=detector.fDetectorCode;
  //Find minimum and maximum fiducial segment numbers using the bincontent
  int minFidBin = hPositionToSegment.at(detNumber)->FindFirstBinAbove(0., 2);
  int maxFidBin = hPositionToSegment.at(detNumber)->FindLastBinAbove(0., 2) + 1;
  return maxFidBin - minFidBin;
}

void TExperiment::GenerateFidBinning(const TDetector &detector,double* fidBinning, int arraySize)
{
  int detNumber=detector.fDetectorCode;// Get detetector code
  if(arraySize!=nBinsFid.at(detNumber)+1) {
    printf("In fiducial binning, arraySize(%i) !=nBinsFid(%i)\n,",arraySize,nBinsFid.at(detNumber)+1);
    exit(1);
  }
  int minFidBin = hPositionToSegment.at(detNumber)->FindFirstBinAbove(0., 2);// Get bin number for minimum fiducial bin
  int maxFidBin = hPositionToSegment.at(detNumber)->FindLastBinAbove(0., 2) + 1;// Get bin number for maximum fiducial bin
  double minFiducial = hPositionToSegment.at(detNumber)->GetYaxis()->GetBinLowEdge(minFidBin);
  double maxFiducial = hPositionToSegment.at(detNumber)->GetYaxis()->GetBinLowEdge(maxFidBin);
  
  // Make array of Bins
  for (int i = 0; i < nBinsFid.at(detNumber)+1; i++) {
    fidBinning[i] = minFiducial + i*(maxFiducial-minFiducial)/nBinsFid.at(detNumber);
  }
}

void TExperiment::GenerateSegBinning(const TDetector &detector,double* segBinning, int arraySize)
{
  if(arraySize!=detector.fNSeg+1) {
    printf("In segment binning, ArraySize(%i) !=fNSeg(%i)\n,",arraySize,detector.fNSeg+1);
    exit(1);
  }
  for (int i = 0; i <= detector.fNSeg; i++) {
    segBinning[i] = -0.5 + i;
  }
}

void TExperiment::GenerateRecoEBinning(double* recoEBinning, int arraySize)
{
  if(arraySize!=nRecoEBins+1) {
    printf("In reconstructed binning, ArraySize(%i) !=nRecoEBins(%i)\n,",arraySize,nRecoEBins+1);
    exit(1);
  }
  for (int i = 0; i <= nRecoEBins; i++) {
    recoEBinning[i] = recoMinEEnergy + i * recoEBinWidth;
  }
}

void TExperiment::GenerateTrueEBinning(double* trueEBinning, int arraySize)
{
  if(arraySize!=nTrueEBins+1) {
    printf("In true binning, ArraySize(%i) !=nTrueEBins(%i)\n,",arraySize,nTrueEBins+1);
    exit(1);
  }
  for (int i = 0; i <= nTrueEBins; i++) {
    trueEBinning[i] = trueMinEEnergy + i * trueEBinWidth;
  }
}



void TExperiment::AddDetector(int detType, int detPosition, double nCycles)
{
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  // Set events from the macros
  PROSPECTMacroInterface.RetrieveValue("Events",fNEvents);
  // Set the number of segments from macro
  TString key="IsReferenceData";
  TString isData;
  PROSPECTMacroInterface.RetrieveValue(key,isData);
  
  // Check whether data is being used or no
  isReferenceData=(isData.CompareTo("YES",TString::kIgnoreCase)==0)?true:false;
  
  fDetectors.push_back(new TDetector(detType,detPosition));// Add a detector
  fDetectors.back()->fNEvents = fNEvents;// Add number of events
  int detNumber=fDetectors.back()->fDetectorCode; //Obtain detector number
  
  // If Data is is being used get the data root file name
  if(isData){
    key.Form("DataRootFile%i",detNumber);
    PROSPECTMacroInterface.RetrieveValue(key,dataFileName);
  }
  
  // Setup default values for the reactor-on and reactor-off exposure times
  fRxOnExposure[fDetectors.back()->fDetectorCode]=1.0*nCycles/7.0;
  fRxOffExposure[fDetectors.back()->fDetectorCode]=1.0*nCycles/7.0;
  // Make sure to grab the reactor-on and reactor-off times from the macro file if they are defined there.
  // The defualt values should not be assigned once this function is called.
  ReadMacroInputs(fDetectors.back()->fDetectorCode);
}

void TExperiment::CheckIfExperimentDefined()
{
  if(!isExperimentDefined){
    printf("TExperiment: Experiment not defined!\n");
    exit(1);
  }
}

void TExperiment::ReadEnergySpectrum(TString energySpectrumFileName)
{
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  // Read input energy spectrum file from macro
  PROSPECTMacroInterface.RetrieveValue("EnergySpectrumFile",energySpectrumFileName);
  
  CheckIfExperimentDefined();
  std::ifstream energySpectrumFile;
  energySpectrumFile.open(energySpectrumFileName.Data());
  
  if (!energySpectrumFile.is_open()){
    printf("Energy spectrum file %s unavailable, check location and try again\n",energySpectrumFileName.Data());
    exit(1);
  }
  printf("Reading spectrum and cross-section information from %s\n",energySpectrumFileName.Data());
  
  // Clear vectors to prepare to fill with energy spectrum information
  fAntiNuEnergy.clear();
  fAntiNuFlux.clear();
  fIBDCrossSection.clear();
  
  std::string line;
  
  // The vector will have a first element of the mean energy
  // of the energy bin in MeV, the second element will be the cross
  // section in the energy bin in cm^2 and the last element will be
  // the number of neutrinos per unit energy (Jx10^?).
  // The values come from Vogel and Engel: Neutrino electromagnetic
  // form factors.  Equation 4 using parameters from table 1 gives
  // the spectrum measurement.
  while(std::getline(energySpectrumFile, line)){
    if (line.find(" ")==0) continue;//skip if the line is empty
    if (line.find("#")==0) continue;// skip if the line starts with #
    
    double energy = 0.0;
    double flux = 0.0;
    double xSection = 0.0;
    
    TString inputLine(line);
    if (inputLine.IsWhitespace()) continue;//skip if the line is empty
    if (inputLine.BeginsWith("#")) continue;// skip if the line starts with #
    
    TObjArray *objArr = inputLine.Tokenize(" ");//Divide the line delimited by =
    assert(objArr->GetEntries() == 3);//Make sure there are 3 entries
    
    energy = ((TObjString*)objArr->At(0))->GetString().Atof();
    xSection = ((TObjString*)objArr->At(1))->GetString().Atof();
    flux = ((TObjString*)objArr->At(2))->GetString().Atof();
    
    fAntiNuEnergy.push_back(energy);
    fAntiNuFlux.push_back(flux / 100.0 / 100.0); // Converting cm^2 into m^2
    fIBDCrossSection.push_back(xSection);
  }
  energySpectrumFile.close();
  // Total number of energy points
  //printf("The total number of Energy Bins used is %lu\n",fAntiNuEnergy.size());
}

void TExperiment::WriteEnergySpectrum(TFile& outputFile)
{
  CheckIfExperimentDefined();
  if(!outputFile.IsOpen() || outputFile.IsZombie())
  {
    printf("File not open");
    exit(1);
  }
  outputFile.cd();
  int count=0;
  for(std::vector<double>::iterator it = fAntiNuEnergy.begin() ; it != fAntiNuEnergy.end(); ++it)
  {
    hAntiNuSpectrum->Fill(*it,fAntiNuFlux.at(count));
    hIBDCrossSection->Fill(*it,fIBDCrossSection.at(count));
    hObsSpectrum->Fill(*it,fAntiNuFlux.at(count)*fIBDCrossSection.at(count));
    count++;
  }
  
  // Normalize the histograms to 1
  hAntiNuSpectrum->Scale(1.0/hAntiNuSpectrum->Integral());
  hIBDCrossSection->Scale(1.0/hIBDCrossSection->Integral());
  hObsSpectrum->Scale(1.0/hObsSpectrum->Integral());
  
  // Writing the histograms to the output file
  hAntiNuSpectrum->Write();
  hIBDCrossSection->Write();
  hObsSpectrum->Write();
  printf("Energy spectrum written \n");
}

void TExperiment::CreateEventTree()
{
  CheckIfExperimentDefined();
  // Instantiate event tree/
  // Change it to it& ?
  for(auto it:fDetectors)
  {
    int detNumber = it->fDetectorCode;
    TString eventTreeName;
    eventTreeName.Form("EventTree%i",detNumber);
    fEventTree[detNumber] = new TTree(eventTreeName.Data(),"Flat Tree of the Interaction Baselines");
    
    // Create branches to the event tree
    fEventTree.at(detNumber)->Branch("ReactorPos", "TVector3", &reactorPos);
    fEventTree.at(detNumber)->Branch("DetectorPos", "TVector3", &detectorPos);
    fEventTree.at(detNumber)->Branch("SegmentId", &segmentId, "SegmentId/I");
    fEventTree.at(detNumber)->Branch("SegmentIdX", &segmentIdX, "SegmentIdX/I");
    fEventTree.at(detNumber)->Branch("SegmentIdZ", &segmentIdZ, "SegmentIdZ/I");
    fEventTree.at(detNumber)->Branch("IsFiducial", &isFiducial, "IsFiducial/O");
    fEventTree.at(detNumber)->Branch("TrueBaseline", &trueBaseline, "TrueBaseline/D");
    fEventTree.at(detNumber)->Branch("Weight", &weight, "Weight/D");
    printf("Event tree for %i created\n",detNumber);
  }
}

void TExperiment::CreateTrees()
{
  CheckIfExperimentDefined();
  CreateEventTree();
  fReactor.CreateReactorEventTree();
  for(auto it:fDetectors)
  {
    it->CreatePositionTree();// Create a segment position and baseline TTree
    it->CreateDetectorEventTree(); // Create a TTree that defines the event local(to detector) and global positions.
  }
  areTreesCreated=true;
}

void TExperiment::ResetEventTreeValues()
{
  CheckIfExperimentDefined();
  reactorPos->SetXYZ(0., 0.,-5.);
  detectorPos->SetXYZ(0., 0., -5.);
  segmentId = -1;
  segmentIdX = -1;
  segmentIdZ = -1;
  isFiducial = true;
  trueBaseline = -1.0;
  weight = 0.0;
}

void TExperiment::FillEventTree()
{
  CheckIfExperimentDefined();
  
  for(auto it:fDetectors)
  {
    // Find the total detector volume.
    double detVolume = it->fNSegX* it->fSegWidthX * it->fSegLength * it->fNSegZ * it->fSegWidthZ;
    
    // Calculate normalization: applies to all events regardless
    // of energy or position
    // This, combined with the position weight, the energy weight
    // and the exposure, gives a weight for the event.
    //    double normalization = (it->fProtonDensity) * (it->fDetectorEfficiency) * (fReactor.fReactorPower) * fRxOnExposure[it->fDetectorCode];
    double normalization = (it->fProtonDensity)  * (fReactor.fReactorPower) * fRxOnExposure[it->fDetectorCode];
    
    // Initialize variables to read trees into
    double reactorX = 0.;
    double reactorY = 0.;
    double reactorZ = 0.;
    
    double detectorX = 0.;
    double detectorY = 0.;
    double detectorZ = 0.;
    
    int segIndex = -1;
    int segIndexX = -1;
    int segIndexZ = -1;
    
    bool detectorFiducial = true;
    
    fReactor.ReadReactorEventTree(reactorX, reactorY, reactorZ);
    
    it->ReadDetectorEventTree(detectorX, detectorY, detectorZ);
    it->ReadDetectorEventTree(segIndex, segIndexX, segIndexZ, detectorFiducial);
    
    // Create all events
    for (int i = 0; i < fNEvents; i++) {
      ResetEventTreeValues();
      
      // Access Reactor event information here
      fReactor.fReactorEventTree->GetEntry(i);
      
      reactorPos->SetXYZ(reactorX, reactorY, reactorZ);
      
      // Access Detector event information here
      it->fDetectorEventTree->GetEntry(i);
      
      detectorPos->SetXYZ(detectorX, detectorY, detectorZ);
      
      // Save the segment information
      segmentId = segIndex;
      segmentIdX = segIndexX;
      segmentIdZ = segIndexZ;
      
      // Save fiducial volume
      isFiducial = detectorFiducial;
      
      // Calculate the True Baseline
      TVector3 baselineVector(detectorPos->X() - reactorPos->X(),
                              detectorPos->Y() - reactorPos->Y(),
                              detectorPos->Z() - reactorPos->Z());
      trueBaseline = baselineVector.Mag();
      
      // Calculate the position weight
      double posWeight = detVolume / fNEvents * 1.0 / (4 * TMath::Pi() * trueBaseline * trueBaseline);
      
      // This weight is ONLY missing the weighting effect from the energy and the Exposure time.
      // This weight is equivalent to one year of running.,.
      weight = normalization * posWeight;
      
      fEventTree.at(it->fDetectorCode)->Fill();
      
    } // End of Creating Event Position
    printf("Event tree for detector %i filled\n",it->fDetectorCode);
  }
}


void TExperiment::FillTrees()
{
  CheckIfExperimentDefined();
  if(!areTreesCreated){
    printf("Create event and position trees before filling them!\n");
    exit(1);
  }
  fReactor.FillReactorEventTree();
  for(auto it:fDetectors)
  {
    it->FillPositionTree();
    it->FillDetectorEventTree();
  }
  FillEventTree();
}

/*
 void TExperiment::ReadEventTree(TString inputFileName)
 {
 CheckIfExperimentDefined();
 if(!areTreesCreated){
 printf("Create event and position trees before filling them!\n");
 exit(1);
 }
 TFile *inputFile = TFile::Open(inputFileName.Data());
 if(!inputFile->IsOpen() || inputFile->IsZombie()){
 printf("The file is not open or is a zombie\n");
 exit(1);
 }
 
 fEventTree = (TTree*)inputFile->Get("EventTree");
 
 // Check that the tree exists
 if (!fEventTree) {
 printf("Incomplete file, tree doesn't exist\n");
 exit(1);
 }
 
 // Set the branch addresses
 fEventTree->SetBranchAddress("ReactorPos", &reactorPos);
 fEventTree->SetBranchAddress("DetectorPos", &detectorPos);
 fEventTree->SetBranchAddress("TrueBaseline", &trueBaseline);
 fEventTree->SetBranchAddress("Weight", &weight);
 
 }
 */


void TExperiment::SetupSignalHists()
{
  CheckIfExperimentDefined();
  int nInputBins=120;
  double minInputBin = 0.0;
  double maxInputBin = 12.0;
  hAntiNuSpectrum = new TH1D("hAntiNuSpectrum","hAntiNuSpectrum",nInputBins,minInputBin,maxInputBin);
  hIBDCrossSection = new TH1D("hIBDCrossSection","hIBDCrossSection",nInputBins,minInputBin,maxInputBin);
  hObsSpectrum= new TH1D("hObsSpectrum","hObsSpectrum",nInputBins,minInputBin,maxInputBin);
  
  for (auto it:fDetectors)
  {
    int detNumber = it->fDetectorCode;
    
    TString histName;
    histName.Form("BaselineToSegment%i",detNumber);
    // Maps the true baseline to the segment
    hBaselineToSegment[detNumber] = new TH2D(histName.Data(),
                                             "Map of Segments to true Baselines; Segment Index; True Baseline (m)",
                                             it->fNSeg, -0.5, it->fNSeg-0.5, it->fNBinsL*5, it->fMinBaseline, it->fMaxBaseline);
    
    // Maps the true baseline to the position
    histName.Form("PositionToSegment%i",detNumber);
    hPositionToSegment[detNumber] = new TH2D(histName.Data(),
                                             "Map of Segments to Positions; Segment Index; Position (m)",
                                             it->fNSeg, -0.5, it->fNSeg-0.5, it->fNBinsL, it->fMinBaseline, it->fMaxBaseline);
    
    // Event mapping of the detector
    histName.Form("PseudoDetector%i",detNumber);
    hPseudoDetector[detNumber] = new TH3D(histName.Data(), "Fake Detector; X (m); Y(m); Z(m)",
                                          it->fNSegX, it->fDetLocation.fDistX, it->fDetLocation.fDistX + it->fNSegX*it->fSegWidth,
                                          106, it->fDetLocation.fDistY - it->fSegLength/2.0, it->fDetLocation.fDistY + it->fSegLength/2.0,
                                          it->fNSegZ, it->fDetLocation.fDistZ, it->fDetLocation.fDistZ + it->fNSegZ*it->fSegWidth);
    printf("histograms for the detector %i setup\n",detNumber);
  }
}

void TExperiment::SetupBackgroundHistograms(TString BKGInputFile)
{
  CheckIfExperimentDefined();
  
  printf("Setting up background from the file %s\n",BKGInputFile.Data());
  // The background comes from Michael Mendenhall and includes detector response.
  
  // Pull Information from file (comes from Michael Mendenhall)
  // The background is in mHz/MeV and normalized over the fiducialized near detector.
  std::ifstream ifs (BKGInputFile.Data());
  
  int segmentIndex = 0;   // The index of the segment
  int segmentIndexX = 0;  // The x index of the segment
  int segmentIndexZ = 0;  // The z index of the segment
  double baseline = 0.0;  // The "measured" baseline to the segment
  
  double recoEBinning[nRecoEBins+1];
  GenerateRecoEBinning(recoEBinning,(sizeof(recoEBinning)/sizeof(*recoEBinning)));
  
  TString histName;
  
  for (auto it:fDetectors)
  {
    int detNumber = it->fDetectorCode;
    
    // Segment bins
    double segBinning[it->fNSeg+1];
    GenerateSegBinning(*it,segBinning,(sizeof(segBinning)/sizeof(*segBinning)));
    
    // Fiducial segment bins
    double fidBinning[nBinsFid.at(detNumber)+1];
    GenerateFidBinning(*it,fidBinning,(sizeof(fidBinning)/sizeof(*fidBinning)));
    
    
    histName.Form("BkgRxOn1D%i",detNumber);
    // This describes the cosmic background PER SEGMENT! in the fiducial volume.
    hBkgRxon1D[detNumber] = new TH1D(histName, "Background; Energy (MeV); Entries (mHz/MeV)",
                                     nRecoEBins, recoEBinning);
    
    histName.Form("BkgRxOnSegvsE%i",detNumber);
    // Create the background histogram for fiducialized volume.
    hBkgRxonSeg[detNumber] = new TH2D(histName,
                                      "Simulated Background; Energy (MeV); Segment Number",
                                      nRecoEBins, recoEBinning, it->fNSeg, segBinning);
    
    histName.Form("BkgRxOnLvsE%i",detNumber);
    hBkgRxOnLvsE[detNumber] = new TH2D(histName,
                                       "Simulated Reactor-on Background; Energy (MeV); Baseline (m)",
                                       nRecoEBins, recoEBinning, nBinsFid.at(detNumber), fidBinning);
    
    histName.Form("BkgRxOffLvsE%i",detNumber);
    hBkgRxOffLvsE[detNumber] = new TH2D(histName,
                                        "Simulated Reactor-off Background; Energy (MeV); Baseline (m)",
                                        nRecoEBins, recoEBinning, nBinsFid.at(detNumber), fidBinning);
  }
  
  //  double elow;
  //  double ehi;
  //  double bkg;
  //  double bkgerr;
  //  double sig;
  //  double sigerr;
  
  std::vector< std::vector<double> > BackgroundValues; // Input is mHz/MeV
  
  if (ifs.is_open()) {
    std::string line;
    
    while (std::getline(ifs, line)){
      TString inputLine(line); // put the lines in Tstring form, makes it easier to manipulate
      
      if (inputLine.IsWhitespace()) continue;//skip if the line is empty
      if (inputLine.BeginsWith("#")) continue;// skip if the line starts with #
      
      TObjArray *objArr = inputLine.Tokenize(" ");//Divide the line delimited by =
      
      TIterator *iter=objArr->MakeIterator(); // define an iterator to iteratre over each word in the line
      TObjString *valueRead;
      
      std::vector<double> temp;
      int columnCount=0;
      
      while((valueRead=(TObjString*)iter->Next())) // While there is another string in the line
      {
        columnCount++;
        if(columnCount>3) continue;
        TString value=valueRead->GetString();
        temp.push_back(value.Atof());
      }
      BackgroundValues.push_back(temp);
    } // end of while loop
  }
  else {
    printf("Background File not available\n");
    exit(1);
  }
  
  for (auto it:fDetectors)
  {
    int detNumber = it->fDetectorCode;
    for (unsigned int i = 0; i < BackgroundValues.size(); i++) {
      double pullE = (BackgroundValues[i][0]+BackgroundValues[i][1])/2.0;
      double eWidth = BackgroundValues[i][1]-BackgroundValues[i][0];
      double bkgWeight = BackgroundValues[i][2] /
      1000.0 * // extra factor converts mHz into 1/yr
      fRxOffExposure[detNumber] * // Include the amount of time exposed.
      eWidth / // Include the bin width in MeV
      108.0; // Include the average over the number of fiducial segments.
      hBkgRxon1D.at(detNumber)->Fill(pullE, bkgWeight);
    }
    
    // Fill the background histogram for each segment.
    it->ReadPositionTree(segmentIndex, segmentIndexX, segmentIndexZ, baseline);
    for (int i = 0; i < it->fNSeg; i++) {
      // Get the measured position
      it->fPositionTree->GetEntry(i);
      double pullL = baseline;
      double pullS = segmentIndex;
      
      // Use fiducialized volume
      if (segmentIndexX == 0 || segmentIndexX == it->fNSegX-1) continue;
      if (segmentIndexZ == 0 || segmentIndexZ == it->fNSegZ-1) continue;
      
      for (int i = 0; i < hBkgRxon1D.at(detNumber)->GetNbinsX(); i++) {
        double pullE = hBkgRxon1D.at(detNumber)->GetBinCenter(i+1);
        double bkgWeight = hBkgRxon1D.at(detNumber)->GetBinContent(i+1);
        hBkgRxOnLvsE.at(detNumber)->Fill(pullE, pullL, bkgWeight);
        hBkgRxOffLvsE.at(detNumber)->Fill(pullE, pullL, bkgWeight);
        hBkgRxonSeg.at(detNumber)->Fill(pullE, pullS, bkgWeight);
      }
    } // End loop through Bkg bins
  }
}

void TExperiment::SetupExperiment()
{
  CheckIfExperimentDefined();// Make sure the experiment is defined.
  ReadEnergySpectrum();//Read energy spectrum from the user-defined file
  CreateTrees();//Create position and event trees
  FillTrees();//Fill position and event trees
  SetupSignalHists();// Setup signal histogram
  
  // Create a true neutrino energy to reconstructed energy map
  // The bins are 0.2 MeV bins between 0 .8 and 7.2 MeV
  // This gives a total of 32 Energy Bins.
  // Construct the true energy binning
  
  double trueEBinning[nTrueEBins+1];
  double recoEBinning[nRecoEBins+1];
  GenerateRecoEBinning(recoEBinning,(sizeof(recoEBinning)/sizeof(*recoEBinning)));
  GenerateTrueEBinning(trueEBinning,(sizeof(trueEBinning)/sizeof(*trueEBinning)));
  
  int segmentIndex = 0;   // The index of the segment
  int segmentIndexX = 0;  // The x index of the segment
  int segmentIndexZ = 0;  // The z index of the segment
  double baseline = 0.0;  // The "measured" baseline to the segment
  
  TString histName;
  for (auto it:fDetectors)
  {
    int detNumber = it->fDetectorCode;
    it->ReadPositionTree(segmentIndex, segmentIndexX, segmentIndexZ, baseline);
    
    for (int i = 0; i < it->fNSeg; i++) {
      it->ResetPositionTreeValues();
      it->fPositionTree->GetEntry(i);
      
      // Use fiducialized volume
      if (segmentIndexX == 0 || segmentIndexX == it->fNSegX-1) continue;
      if (segmentIndexZ == 0 || segmentIndexZ == it->fNSegZ-1) continue;
      
      hPositionToSegment.at(detNumber)->Fill(segmentIndex, baseline);
    }
    
    nBinsFid[detNumber]=GetNFidBins(*it);
    double fidBinning[nBinsFid.at(detNumber)+1];
    GenerateFidBinning(*it,fidBinning,(sizeof(fidBinning)/sizeof(*fidBinning)));
    
    histName.Form("PositionToSegmentFid%i",detNumber);
    hPositionToSegmentFid[detNumber] = new TH2D(histName,
                                                "Map of Segments to Positions; Segment Index; Position (m)",
                                                it->fNSeg, -0.5, it->fNSeg-0.5, nBinsFid.at(detNumber),fidBinning);
    
    // Fill the position to segment map
    // First get the total position weight
    for (int i = 0; i < it->fNSeg; i++) {
      it->fPositionTree->GetEntry(i);
      
      // Use fiducialized volume
      if (segmentIndexX == 0 || segmentIndexX == it->fNSegX-1) continue;
      if (segmentIndexZ == 0 || segmentIndexZ == it->fNSegZ-1) continue;
      
      hPositionToSegmentFid.at(detNumber)->Fill(segmentIndex, baseline);
    }
    
    // Fill the PseudoDetector Histogram with the positions of the generated events.
    for (int ev = 0; ev < fNEvents; ev++) {
      fEventTree.at(detNumber)->GetEvent(ev);
      
      hPseudoDetector.at(detNumber)->Fill(detectorPos->X(), detectorPos->Y(), detectorPos->Z(), weight);
      
      //Calculate Segment
      int binX = 0;
      int binY = 0;
      int binZ = 0;
      int globalBin = hPseudoDetector.at(detNumber)->FindBin(detectorPos->X(), detectorPos->Y(), detectorPos->Z());
      hPseudoDetector.at(detNumber)->GetBinXYZ(globalBin, binX, binY, binZ);
      //int segment = (binX-1) + (binZ-1) * it->fNSegX;
      
      // Perform Fiducial Cuts
      if (!isFiducial) continue;
      //      if (binX == 1 || binX == it->fNSegX) continue;
      //      if (binZ == 1 || binZ == it->fNSegZ) continue;
      //      if (TMath::Abs(detectorPos->Y()) > it->fSegLength/2.0-it->fSegWidth) continue;
      hBaselineToSegment.at(detNumber)->Fill(segmentId, trueBaseline, weight);
    }
    //hSimulatedNueEnergy
    histName.Form("SimulatedNueEnergy%i",detNumber);
    hSimulatedNueEnergy[detNumber] = new TH1D(histName, "Simulated #nu Energy; Simulated Energy (MeV); Counts",nTrueEBins, trueEBinning);
    
    // Fill the fine binned histogram, we are doing small steps to more accurately portray the smearing.
    for (unsigned int i = 0; i < fAntiNuEnergy.size(); i++) {
      double trueEnergy = fAntiNuEnergy[i];
      double eneWeight = fAntiNuFlux[i]*fIBDCrossSection[i];
      hSimulatedNueEnergy.at(detNumber)->Fill(trueEnergy,eneWeight);
    }
    
    TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
    // Setup background hists at this point
    TString bkgFile("./inputs/SimulatedBackground.txt");
    
    // Set background file from the macro
    PROSPECTMacroInterface.RetrieveValue("BKGInputFile",bkgFile);
    SetupBackgroundHistograms(bkgFile.Data());
    
    for (auto it:fDetectors)
    {
      int detNumber = it->fDetectorCode;
      it->ReadPositionTree(segmentIndex, segmentIndexX, segmentIndexZ, baseline);
      // Create a segment vs E plot used in covariance matrix calculations.
      double segBinning[it->fNSeg+1];
      GenerateSegBinning(*it,segBinning,(sizeof(segBinning)/sizeof(*segBinning)));
      
      double fidBinning[nBinsFid.at(detNumber)+1];
      GenerateFidBinning(*it,fidBinning,(sizeof(fidBinning)/sizeof(*fidBinning)));
      
      histName.Form("ETrue%i",detNumber);
      hETrue[detNumber] = new TH1D(histName,"Absolute unoscillated MC #nu spectrum; Neutrino Energy (MeV); Counts",nTrueEBins, trueEBinning);
      histName.Form("SegvsETrue%i",detNumber);
      hSegvsETrue[detNumber] = new TH2D(histName,"Unoscillated MC #nu spectrum for each segment; Neutrino Energy (MeV); Segment Index",nTrueEBins, trueEBinning, it->fNSeg, segBinning);
      histName.Form("ENull%i",detNumber);
      hENull[detNumber] = new TH1D(histName,"Absolute unoscillated MC deposited spectrum; Prompt Energy (MeV); Counts",nRecoEBins, recoEBinning);
      histName.Form("SegvsENull%i",detNumber);
      hSegvsENull[detNumber] = new TH2D(histName,"Unoscillated MC deposited spectrum for each segment; Prompt Energy (MeV); Segment Index",nRecoEBins, recoEBinning, it->fNSeg, segBinning);
      histName.Form("LvsENull%i",detNumber);
      hLvsENull[detNumber] = new TH2D(histName, "Unoscillated MC #nu spectrum for each baseline; Prompt Energy (MeV); Measured Baseline(m)",nRecoEBins, recoEBinning, nBinsFid.at(detNumber), fidBinning);
      histName.Form("LvsERelative%i",detNumber);
      hLvsERelative[detNumber] = new TH2D(histName, "Unoscillated scaled MC deposited spectrum for each baseline; Prompt Energy (MeV); Measured Baseline(m)",nRecoEBins, recoEBinning, nBinsFid.at(detNumber), fidBinning);
      
      histName.Form("EOsc%i",detNumber);
      hEOsc[detNumber] = new TH1D(histName,"Absolute oscillated MC #nu spectrum; Neutrino Energy (MeV); Counts",nTrueEBins, trueEBinning);
      histName.Form("SegvsEOsc%i",detNumber);
      hSegvsEOsc[detNumber] = new TH2D(histName,"Oscillated MC #nu spectrum for each segment; Neutrino Energy (MeV); Segment Index",nTrueEBins, trueEBinning, it->fNSeg, segBinning);
      
      histName.Form("ERef%i",detNumber);
      hERef[detNumber] = new TH1D(histName,"Absolute unoscillated reconstructed spectrum; Prompt Energy (MeV); Counts",nRecoEBins, recoEBinning);
      histName.Form("SegvsERef%i",detNumber);
      hSegvsERef[detNumber] = new TH2D(histName,"Absolute reconstructed spectrum at each segment; Prompt Energy (MeV); Segment Index",nRecoEBins, recoEBinning, it->fNSeg, segBinning);
      histName.Form("LvsERef%i",detNumber);
      hLvsERef[detNumber] = new TH2D(histName, "Absolute reconstructed spectrum at each baseline; Prompt Energy (MeV); Measured Baseline(m)", nRecoEBins, recoEBinning, nBinsFid.at(detNumber), fidBinning);
      histName.Form("LvsERelativeRef%i",detNumber);
      hLvsERelativeRef[detNumber] = new TH2D(histName, "Absolute scaled reconstructed spectrum at each baseline; Prompt Energy (MeV); Measured Baseline(m)", nRecoEBins, recoEBinning, nBinsFid.at(detNumber), fidBinning);
    }
    
    // First get the total position weight
    for (int i = 0; i < it->fNSeg; i++) {
      double posWeight = hBaselineToSegment.at(detNumber)->Integral(i+1, i+1, 0, -1);
      it->ResetPositionTreeValues();
      it->fPositionTree->GetEntry(i);
      //        double position = baseline;
      if (posWeight == 0) continue;
      // Now get the total reconstructed energy weight
      for (int j = 0; j < nTrueEBins; j++) {
        double trueEnergy= hSimulatedNueEnergy.at(detNumber)->GetBinCenter(j+1);
        double eneWeight= hSimulatedNueEnergy.at(detNumber)->GetBinContent(j+1);
        
        double totWeight = posWeight*eneWeight;
        hSegvsENull.at(detNumber)->Fill(trueEnergy, segmentIndex, totWeight);
        hSegvsETrue.at(detNumber)->Fill(trueEnergy, segmentIndex, totWeight);
        hETrue.at(detNumber)->Fill(trueEnergy, totWeight);
      }
    }
    printf("Applying detector response to the simulated histograms\n");
    TDetectorResponseInterface fDetRespInt(detNumber);
    // modification by xiaobin lu 2018.2.16

    // fDetRespInt.ApplyDetectorResponse(*hETrue.at(detNumber),hENull.at(detNumber));
    // fDetRespInt.ApplyDetectorResponse(*hSegvsETrue.at(detNumber),hSegvsENull.at(detNumber));
    fDetRespInt.ApplyDetectorResponse(*hSegvsETrue.at(detNumber),hSegvsENull.at(detNumber),hENull.at(detNumber));




    //Fill L vs E histogram from the Seg vs E histogram
    hLvsENull.at(detNumber)->Scale(0);// Make sure there is nothing in the LvsENull histograms
    for (int i = 0; i < it->fNSeg; i++) {
      it->ResetPositionTreeValues();
      it->fPositionTree->GetEntry(i);
      // Get the measured position
      it->fPositionTree->GetEntry(i);
      
      for (int j = 0; j < hSegvsENull.at(detNumber)->GetNbinsX(); j++) {
        double energy = hSegvsENull.at(detNumber)->GetXaxis()->GetBinCenter(j+1);
        hLvsENull.at(detNumber)->Fill(energy, baseline, hSegvsENull.at(detNumber)->GetBinContent(j+1,i+1));
      }
    }
    
    RelativizeHistogram(*hLvsENull.at(detNumber),*hENull.at(detNumber),*hLvsENull.at(detNumber),*hENull.at(detNumber),*hLvsERelative.at(detNumber));
    
    // If data is used for oscillation analysis the signal and background histograms have to come from data
    if(isReferenceData){
      hLvsERef.at(detNumber)->Scale(0.0);
      TString histName;
      TFile *dataFile=TFile::Open(dataFileName.Data());
      
      // Obtain the data histogram from the data files
      TH1D *tempHist1D=(TH1D*)dataFile->Get("E_IBD");
      TH2D *tempHist2D=(TH2D*)dataFile->Get("LvsE_IBD");
      
      //Check to make sure RxOn IBD data hists are consistent with the MC hists
      assert(CheckConsistency(*tempHist2D,*hLvsENull.at(detNumber)));
      assert(CheckConsistency(*tempHist1D,*hENull.at(detNumber)));
      
      tempHist1D->Copy(*hERef.at(detNumber));
      tempHist2D->Copy(*hLvsERef.at(detNumber));
      
      tempHist2D=(TH2D*)dataFile->Get("SegvsE_IBD");
      tempHist2D->Copy(*hSegvsERef.at(detNumber));
      //Clear the exisiting hists
      hBkgRxon1D.at(detNumber)->Scale(0);
      hBkgRxOffLvsE.at(detNumber)->Scale(0);
      hBkgRxOffLvsE.at(detNumber)->Scale(0);
      
      tempHist1D=(TH1D*)dataFile->Get("E_BG");
      tempHist2D=(TH2D*)dataFile->Get("LvsE_BG");
      assert(CheckConsistency(*tempHist2D,*hLvsENull.at(detNumber)));
      assert(CheckConsistency(*tempHist1D,*hENull.at(detNumber)));
      
      tempHist1D->Copy(*hBkgRxon1D.at(detNumber));
      tempHist2D->Copy(*hBkgRxOffLvsE.at(detNumber));
      tempHist2D->Copy(*hBkgRxOnLvsE.at(detNumber));
      
      histName.Form("ERef%i",detNumber);
      hERef.at(detNumber)->SetName(histName);
      histName.Form("LvsERef%i",detNumber);
      hLvsERef.at(detNumber)->SetName(histName);
      histName.Form("BkgRxOffLvsE%i",detNumber);
      hBkgRxOffLvsE.at(detNumber)->SetName(histName);
      histName.Form("BkgRxOnLvsE%i",detNumber);
      hBkgRxOnLvsE.at(detNumber)->SetName(histName);
      
      tempHist2D=(TH2D*)dataFile->Get("SegvsE_BG");
      
      // Subtract backgrounds, this obviously has to be done in a much cleaner way
      hLvsERef.at(detNumber)->Add(hBkgRxOnLvsE.at(detNumber),-1);
      hERef.at(detNumber)->Add(hBkgRxon1D.at(detNumber),-1);
      hSegvsERef.at(detNumber)->Add(tempHist2D,-1);
      
      RelativizeHistogram(*hLvsENull.at(detNumber),*hENull.at(detNumber),*hLvsERef.at(detNumber),*hERef.at(detNumber),*hLvsERelativeRef.at(detNumber));
    }// end of if loop
    
    // Else the histograms are obtained from the histograms populated by the simulated TTrees
    else{
      it->ResetPositionTreeValues();
      for (int i = 0; i < it->fNSeg; i++) {
        // Get position weight of the tree
        double posWeight = hBaselineToSegment.at(detNumber)->Integral(i+1, i+1, 0, -1);
        it->fPositionTree->GetEntry(i);
        //        double position = baseline;
        if (posWeight == 0) continue;
        // Now get the total reconstructed energy weight
        for (int j = 0; j < nTrueEBins; j++) {
          double trueEnergy= hSimulatedNueEnergy.at(detNumber)->GetBinCenter(j+1);
          double eneWeight= hSimulatedNueEnergy.at(detNumber)->GetBinContent(j+1);
          double totWeight = posWeight*eneWeight;
          
          // Calculate the Delta m2 contribution
          double dm2Term = fOscillationSimulator.GetDeltam2Term(fReferenceDeltam2,trueEnergy,baseline);
          hSegvsEOsc.at(detNumber)->Fill(trueEnergy, segmentIndex, totWeight*dm2Term);
          hEOsc.at(detNumber)->Fill(trueEnergy, totWeight*dm2Term);
        }
      }
      
      // Apply detector response
      fDetRespInt.ApplyDetectorResponse(*hSegvsEOsc.at(detNumber),hSegvsERef.at(detNumber));
      fDetRespInt.ApplyDetectorResponse(*hEOsc.at(detNumber),hERef.at(detNumber));
      // Make sure to scale the spectrum with fReferenceSin22theta, this gives disappeared number of neutrios
      hSegvsERef.at(detNumber)->Scale(-fReferenceSin22theta);
      hERef.at(detNumber)->Scale(-fReferenceSin22theta);
      // Subtract from total neutrinos to get surviving neutrinos
      hERef.at(detNumber)->Add(hENull.at(detNumber));
      hSegvsERef.at(detNumber)->Add(hSegvsENull.at(detNumber));
      
      
      //Fill L vs E histogram from the Seg vs E histogram
      for (int i = 0; i < it->fNSeg; i++) {
        it->ResetPositionTreeValues();
        it->fPositionTree->GetEntry(i);
        // Get the measured position
        it->fPositionTree->GetEntry(i);
        
        for (int j = 0; j < hSegvsERef.at(detNumber)->GetNbinsX(); j++) {
          double energy = hSegvsERef.at(detNumber)->GetXaxis()->GetBinCenter(j+1);
          hLvsERef.at(detNumber)->Fill(energy, baseline, hSegvsERef.at(detNumber)->GetBinContent(j+1,i+1));
        }
      }
      
      RelativizeHistogram(*hLvsENull.at(detNumber),*hENull.at(detNumber),*hLvsERef.at(detNumber),*hERef.at(detNumber),*hLvsERelativeRef.at(detNumber));
      
      // Set background on to background off before scaling
      if(fRxOnOffRatio.at(detNumber)==0){
        printf("Ratio of reactor off to on cannot be %f \n",fRxOnOffRatio.at(detNumber));
        exit(1);
      }
      hBkgRxOffLvsE.at(detNumber)->Scale(1/fRxOnOffRatio.at(detNumber));
    }// end of else loop
    
    // Check the signal to background
    // Work between 0.8 MeV and 7.2 MeV
    double lowEBin = 1;
    double highEBin = 32;
    double lowPBin = hLvsENull.at(detNumber)->FindFirstBinAbove(10,2);
    double highPBin = hLvsENull.at(detNumber)->FindLastBinAbove(10,2);
    
    
    //EDIT: Signal should not be null histogram
    // Calculate the normalization of each histogram.
    double nEventsNull = hLvsENull.at(detNumber)->Integral(lowEBin, highEBin, lowPBin, highPBin);
    // In mHz/MeV with 0.2 MeV bins.
    double nEventsBkg = hBkgRxOnLvsE.at(detNumber)->Integral(lowEBin,highEBin,lowPBin, highPBin);
    printf(" Number of Background: %f\n",nEventsBkg);
    // Ensure far detector has same S:B as near.
    if (it->fDetLocation.fDetectorType == 2) {
      double SignalToBackground = 3.25;
      hBkgRxOnLvsE.at(detNumber)->Scale(nEventsNull / SignalToBackground / nEventsBkg);
      hBkgRxonSeg.at(detNumber)->Scale(nEventsNull / SignalToBackground / nEventsBkg);
      nEventsBkg = hBkgRxOnLvsE.at(detNumber)->Integral(lowEBin, highEBin, lowPBin, highPBin);
    }
    printf("Details of detector %i setup.\n",it->fDetectorCode);
    
    printf(" Min E: %f\n",hBkgRxOnLvsE.at(detNumber)->GetXaxis()->GetBinLowEdge(lowEBin));
    printf(" Max E: %f\n",hBkgRxOnLvsE.at(detNumber)->GetXaxis()->GetBinLowEdge(highEBin) + hBkgRxOnLvsE.at(detNumber)->GetXaxis()->GetBinWidth(highEBin));
    printf(" Number of Signal: %8.3f\n",nEventsNull);
    printf(" Number of Background: %f\n",nEventsBkg);
    printf(" double check S/B : %8.3f\n",nEventsNull / nEventsBkg);
  }
}


void TExperiment::WriteEventTree(TFile& outputFile)
{
  CheckIfExperimentDefined();
  outputFile.cd();
  for(auto it:fDetectors)
  {
    fEventTree.at(it->fDetectorCode)->Write();
  }
  printf("Event trees written to %s \n",outputFile.GetName());
}

void TExperiment::WriteTrees(TFile& outputFile)
{
  CheckIfExperimentDefined();
  WriteEventTree(outputFile);
  fReactor.WriteReactorEventTree(outputFile);
  for(auto it:fDetectors)
  {
    it->WritePositionTree(outputFile);
    it->WriteDetectorEventTree(outputFile);
  }
}

void TExperiment::WriteHistograms(TFile& outputFile)
{
  CheckIfExperimentDefined();
  outputFile.cd();
  // setup all the histograms
  WriteEnergySpectrum(outputFile);
  for (auto it:fDetectors)
  {
    int detNumber = it->fDetectorCode;
    
    hPositionToSegment.at(detNumber)->Write();
    hPositionToSegmentFid.at(detNumber)->Write();
    hBaselineToSegment.at(detNumber)->Write();
    hPseudoDetector.at(detNumber)->Write();
    
    hETrue.at(detNumber)->Write();
    hENull.at(detNumber)->Write();
    hSegvsETrue.at(detNumber)->Write();
    hSegvsENull.at(detNumber)->Write();
    hLvsENull.at(detNumber)->Write();
    hLvsERelative.at(detNumber)->Write();
    
    hERef.at(detNumber)->Write();
    hEOsc.at(detNumber)->Write();
    hSegvsERef.at(detNumber)->Write();
    hSegvsEOsc.at(detNumber)->Write();
    hLvsERef.at(detNumber)->Write();
    hLvsERelativeRef.at(detNumber)->Write();
    
    hBkgRxon1D.at(detNumber)->Write();
    hBkgRxonSeg.at(detNumber)->Write();
    hBkgRxOnLvsE.at(detNumber)->Write();
    hBkgRxOffLvsE.at(detNumber)->Write();
  }
}

void TExperiment::Print()
{
  CheckIfExperimentDefined();
  printf("Details of the experiment used:\n");
  fReactor.Print();
  for (auto it:fDetectors)
  {
    it->Print();
    printf("Experiment : Exposure of the experiment %f\n",fRxOnExposure[it->fDetectorCode]);
  }
  printf("Experiment : Number of events created %i\n",fNEvents);
}
