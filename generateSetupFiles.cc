///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Jan, 2017
// Executable to generate setup root files 
////////////////////////////////////////////////////////////////////////

#include "TDetector.hh"
#include "TReactor.hh"
#include "TExperiment.hh"
#include "TOscillationModelBuilder.hh"
#include "TMinimization.hh"
#include "TMacroInterface.hh"

//Please make sure to follow the steps in the order shown below
//TODO: Need to have better error-checking if the execution do not follow these steps

void usage(int arguments)
{
  printf("Cannot use %i arguments\n",arguments);
  printf("Usage:\n");
  printf("$ ./generateSetupFiles ./macro.mac\n");
}

int main(int argc,char** argv)
{
  TString inputFile("");
  if(argc==1)inputFile.Append("mac/default.mac");
  else if(argc==2)inputFile = (argv[1]);
  else usage(argc);
  
	printf("Using %s macro file.\n",inputFile.Data());
	
	TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  PROSPECTMacroInterface.Initialize(inputFile);
  
  TString outputSetupFileName("set.root");
  TString outputOscillationFileName("Osc.root");
  TString outputMinimizationFileName("Min.root");
  
  PROSPECTMacroInterface.RetrieveValue("SetupFileName",outputSetupFileName);
  
  TExperiment experiment(0,0);
  
  std::map<int,TH2D*> hRef;
  
  // Go through the inputs from macro and check for what detectors to use for the Experiment
  for(std::pair<TString,TString> keyValue:(TMacroExtractor::Instance().GetKeyValueMap())){
    TString key = keyValue.first;//Get the keys of the PROSPECTKeyValuePairs
    TString value = keyValue.second;//Get the values of the PROSPECTKeyValuePairs
    if(key.Contains("DETECTOR")){
      if(value.CompareTo("YES",TString::kIgnoreCase)!=0) continue; // Continue only if the detector is set to yes
      key.ReplaceAll("DETECTOR","");
      int detNumber = key.Atoi(); // Obtain the detector number from the key
      int detLocation = detNumber%100; // Obtain the detector location from the value
      detNumber = detNumber/100;
      experiment.AddDetector(detNumber,detLocation); // Add a detector to the experiment
    }
  }
  // Make sure there is at least one detector with non-zero run-time
  if((experiment.fDetectors).size()==0){
    printf("No detectors added\n");
    printf("Add 'AddDetector = 101' for adding near detector at near position \n");
    exit(1);
  }
  
  TFile* outputSetupFile = new TFile(outputSetupFileName,"RECREATE"); 
  experiment.SetupExperiment();
		
	experiment.fReactor.Print();
  experiment.Print();
  for(auto it:experiment.fDetectors){
		//it->Print();
		it->fDetLocation.Print();
	}

  experiment.WriteTrees(*outputSetupFile);//Write position and event trees to the outputSetupFile
  experiment.WriteHistograms(*outputSetupFile);//Write histograms to the outputSetupFile
  
  outputSetupFile->Close();
  return 0;
}

