#include "TDataExtractor.hh"

TDataExtractor::TDataExtractor(TString extName, TExperiment &experiment, TDetector &detector):fExperiment(experiment),fDetector(detector),extractorName(extName.Data())
{
  fDetCode=fDetector.fDetectorCode;
  /// TChain for IBD-like events
  /// Assumes that the IBD ttree name is Tibd and it is in a TDirectoryFile P2kIBDPlugin
  fInputDataTChain= new TChain("P2kIBDPlugin/Tibd");
  
  printf("For %s extractor, using detector with code %i for extracting baselines\n",extractorName.Data(), fDetCode);
  
  isTDetectorDefined=true;
}

void TDataExtractor::AddTTree(TString dataDirectoryName)
{
  if(isFinal){
    printf("Cannot add a TTree after extracting data, skipping files in %s",dataDirectoryName.Data());
  }
  // Open the system directory for the ROOT file
  dataDirectoryName.Append("/");
  TSystemDirectory dir("IBDDataDir",dataDirectoryName);
  // List all the files in the directory
  TCollection *listOfFiles=(TCollection*)dir.GetListOfFiles();
  
  printf("Adding TTrees from %s \n",dataDirectoryName.Data());
  if(listOfFiles){
    
    TSystemFile *file;
    // Make an iterator of the list of TSystemFile
    TIter next(listOfFiles);
    while ((file=(TSystemFile*)next())) {
      TString fname = file->GetName();
      // continue if the file is not a directory or if it is . or .. . Also make sure it ends with .root
      
      // Search inside subdirectories
      if(file->IsDirectory()&& !(fname.CompareTo(".")==0) && !(fname.CompareTo("..")==0)){
        
        TString fsubdirname(dataDirectoryName);
        fsubdirname.Append(fname);
        fsubdirname.Append("/");
        
        TSystemDirectory subdir("IBDDataSubDir",fsubdirname);
        TCollection *listOfsubFiles=(TCollection*)subdir.GetListOfFiles();
        
        if(listOfsubFiles){
          //std::cout << "List of subfiles found!" << std::endl;
          
          TSystemFile *subfile;
          // Make an iterator of the list of TSystemFile
          TIter next2(listOfsubFiles);
          while ((subfile=(TSystemFile*)next2())) {
            TString fsubname = subfile->GetName();
            
            //std::cout << "Subfile name = " << fsubname.Data() << std::endl;
            // continue if the file is not a directory or if it is . or .. . Also make sure it ends with .root
            if(!(fsubname.EndsWith("P2kAnalyzer.root"))||subfile->IsDirectory()||fsubname.CompareTo(".")==0 ||fsubname.CompareTo("..")==0 ||fsubname.EndsWith(".xml")) continue;
            // Location of the datafile
            TString fileLocation(fsubdirname);
            // Name of the file
            fileLocation.Append(fsubname);
            
            // Add the root file to the TTree
            fInputDataTChain->Add(fileLocation.Data());
            // Weird thing you have to do for SetBranchAddress to work properly later on
            // Commmented out later since this is not needed apparently, will leave it here in case ROOT acts weird again
            //fInputDataTChain->SetBranchStatus("*",1);
            
          }
        }
      }
      
      else{
        if(!(fname.EndsWith("P2kAnalyzer.root"))||file->IsDirectory()||fname.CompareTo(".")==0 ||fname.CompareTo("..")==0 ||fname.EndsWith(".xml")) continue;
        // Location of the datafile
        TString fileLocation(dataDirectoryName);
        // Name of the file
        fileLocation.Append(fname);
        std::cout << "File location = " << fileLocation.Data() << std::endl;
        
        // Add the root file to the TTree
        fInputDataTChain->Add(fileLocation.Data());
        // Weird thing you have to do for SetBranchAddress to work properly later on
        // Commmented out later since this is not needed apparently, will leave it here in case ROOT acts weird again
        //fInputDataTChain->SetBranchStatus("*",1);
      }
    }
  }
}

void TDataExtractor::FinalizeTTreeAddition()
{
  isFinal=true;
  // Assign TChain to the TTree
  fInputDataTree=fInputDataTChain;
  nFiles=fInputDataTChain->GetNtrees();
}

void TDataExtractor::ReadDataTree()
{
  // Set the branch addresses
  
  fInputDataTree->SetBranchAddress("evt", &evt);
  fInputDataTree->SetBranchAddress("t_abs", &t_abs);
  fInputDataTree->SetBranchAddress("E", &E);
  fInputDataTree->SetBranchAddress("maxseg", &maxseg);
  fInputDataTree->SetBranchAddress("xyz", &xyz);
  fInputDataTree->SetBranchAddress("E_maxseg", &E_maxseg);
  fInputDataTree->SetBranchAddress("E_adjacent", &E_adjacent);
  fInputDataTree->SetBranchAddress("segmult", &segmult);
  fInputDataTree->SetBranchAddress("diameter", &diameter);
  fInputDataTree->SetBranchAddress("tspread", &tspread);
  fInputDataTree->SetBranchAddress("ncapt_dt", &ncapt_dt);
  fInputDataTree->SetBranchAddress("n_seg", &n_seg);
  fInputDataTree->SetBranchAddress("n_xyz", &n_xyz);
  fInputDataTree->SetBranchAddress("veto_t", &veto_t);
  fInputDataTree->SetBranchAddress("cut", &cut);
  fInputDataTree->SetBranchAddress("detgeom", &detgeom);
  fInputDataTree->SetBranchAddress("rxpwr", &rxpwr);
  //bLineBranch = fInputDataTree->Branch("bLine", &bLine, "bLine/D");
}

void TDataExtractor::CreateHistograms(){
  TString histName;
  // If the baselines come from TDetector
  if(isTDetectorDefined){
    histName.Form("LvsE_%s",extractorName.Data());
    hLvsE=(TH2D*)(fExperiment.hLvsENull.at(fDetCode)->Clone(histName));
    hLvsE->Scale(0);
    
    histName.Form("SegvsE_%s",extractorName.Data());
    hSegvsE=(TH2D*)(fExperiment.hSegvsENull.at(fDetCode)->Clone(histName));
    hSegvsE->Scale(0);
    
    histName.Form("E_%s",extractorName.Data());
    hE=(TH1D*)(fExperiment.hENull.at(fDetCode)->Clone(histName));
    hE->Scale(0);
    
    for(int i =0;i<fDetector.fNSeg;i++){
      int segX = i%fDetector.fNSegX;
      int segZ = std::floor(i/(fDetector.fNSegX));
      if(segX<fDetector.fNOuterXLeft || segX>=fDetector.fNSegX-fDetector.fNOuterXRight) continue;
      if(segZ<fDetector.fNOuterZBottom || segZ>=fDetector.fNSegZ-fDetector.fNOuterZTop) continue;
      histName.Form("SegE%i_%s",i,extractorName.Data());
      hSegE[i]=(TH1D*)(fExperiment.hENull.at(fDetCode)->Clone(histName));
      hSegE[i]->Scale(0);
    }
  }
}

void TDataExtractor::ExtractData()
{
  
  FinalizeTTreeAddition();
  CreateHistograms();
  ReadDataTree();
  fInputDataTree->GetEntry(0);
  
  int segmentIndex = 0;   // The index of the segment
  int segmentIndexX = 0;  // The x index of the segment
  int segmentIndexZ = 0;  // The z index of the segment
  double baseline = 0.0;  // The "measured" baseline to the
  
  fDetector.ReadPositionTree(segmentIndex, segmentIndexX, segmentIndexZ, baseline);
  
  printf("Extracting data from the extractor %s\n",extractorName.Data());
  for (int i=0;i<fInputDataTree->GetEntries();i++){
    fInputDataTree->GetEntry(i);
    
    if(cut!=3) continue;
    // conversion from P2X numbering scheme to OscSens framework scheme
    int segX = (fDetector.fNSegX-1-(maxseg)%fDetector.fNSegX)%fDetector.fNSegX;
    int segZ = std::floor(maxseg/(fDetector.fNSegX));
    
    int trueSeg=segZ*fDetector.fNSegX+segX;
    
    // Get access to position tree.
    // The segment numbering is different between P2X, PG4 and OscSens
    fDetector.fPositionTree->GetEntry(trueSeg);
    
    // Only fill fiducialized bins
    if(segX<fDetector.fNOuterXLeft || segX>=fDetector.fNSegX-fDetector.fNOuterXRight) continue;
    if(segZ<fDetector.fNOuterZBottom || segZ>=fDetector.fNSegZ-fDetector.fNOuterZTop) continue;
    
    hLvsE->Fill(E,baseline);
    hE->Fill(E);
    hSegE.at(trueSeg)->Fill(E);
    hSegvsE->Fill(E,trueSeg);
  }
  // Weird conversion from reactor files to seconds
  //This has to change for each detector and RxOn/RxOff files
  // This has to be done in a whole different
  fT=nFiles*3600;
  isDataExtracted=true;
  printf("Data successfully extracted for %s\n",extractorName.Data());
}

void TDataExtractor::WriteHistograms(TFile& outputFile)
{
  outputFile.cd();
  hLvsE->Write();
  hE->Write();
  hSegvsE->Write();
  for (auto& it:hSegE) {
    it.second->Write();
  }
}
