#include "TDataExtractor.hh"
#include "TMacroInterface.hh"

void usage(int arguments)
{
  printf("Cannot use %i arguments\n",arguments);
  printf("Usage:\n");
  printf("$ ./test \n");
}

int main(int argc,char** argv)
{
  TString inputFile("");
  if(argc==1)inputFile.Append("mac/extract.mac");
  else if(argc==2)inputFile=argv[1];
  else usage(argc);
  
  // Create singleton instance of Macro interface for extracting information from macros
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  
  // This step is very important to be able to use the macro files
  PROSPECTMacroInterface.Initialize(inputFile);
  
  
  // Go through the inputs from macro and check for input and output file names and locations
  
  TString outputFileName;
  TString outputFileName_mac("outputFileName");
  PROSPECTMacroInterface.RetrieveValue(outputFileName_mac,outputFileName);
  
  TString IBDDataFileLocation;
  TString IBDDataFileLocation_mac("IBDDataFileLocation");
  PROSPECTMacroInterface.RetrieveValue(IBDDataFileLocation_mac,IBDDataFileLocation);
  
  TString BGDataFileLocation;
  TString BGDataFileLocation_mac("BGDataFileLocation");
  PROSPECTMacroInterface.RetrieveValue(BGDataFileLocation_mac,BGDataFileLocation);
  
  TExperiment experiment(0,0);
  
  // Go through the inputs from macro and check for what detectors to use for the Experiment
  for(std::pair<TString,TString> keyValue:(TMacroExtractor::Instance().GetKeyValueMap())){
    TString key = keyValue.first;//Get the keys of the PROSPECTKeyValuePairs
    TString value = keyValue.second;//Get the values of the PROSPECTKeyValuePairs
    if(!key.Contains("DETECTOR")) continue;
    if(value.CompareTo("YES",TString::kIgnoreCase)!=0) continue; // Continue only if the detector is set to yes
    key.ReplaceAll("DETECTOR","");
    int detCode = key.Atoi(); // Obtain the detector number from the key
    int detLocation = detCode%100; // Obtain the detector location from the value
    int detNumber = detCode/100;
    printf("Adding detector %i\n", detCode);
    experiment.AddDetector(detNumber,detLocation); // Add a detector to the experiment
  }
  
  experiment.SetupExperiment();
  
  
  TFile *outputFile=new TFile(outputFileName,"RECREATE");
  for(auto it:experiment.fDetectors){
    // BG and Signal files are different and so need two extractors to extarct IBD or signal data
    TDataExtractor *fIBDDataExtractor=new TDataExtractor("IBD",experiment,*it);
    TDataExtractor *fBGDataExtractor=new TDataExtractor("BG",experiment,*it);
    
    fIBDDataExtractor->AddTTree(IBDDataFileLocation); // Call the Add TTree function in the TDataExtractor
    fBGDataExtractor->AddTTree(BGDataFileLocation);// Call the Add TTree function in the TDataExtractor
    
    //TSystemDirectory BGDir("BGDataParentDir",BGDataFileLocation);
    fIBDDataExtractor->ExtractData();
    fBGDataExtractor->ExtractData();
    
    fIBDDataExtractor->WriteHistograms(*outputFile);
    fBGDataExtractor->WriteHistograms(*outputFile);
  }
  
  outputFile->Close();
  return 0;
}


