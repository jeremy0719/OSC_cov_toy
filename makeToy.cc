///////////////////////////////////////////////////////////////////////
// Author: K Gilje adapted from P T Surukuchi
// Date: Oct 20, 2016
// Executable to throw a single toy.
// Will be expanded on to do Feldman Cousins in the future.
////////////////////////////////////////////////////////////////////////

#include "TDetector.hh"
#include "TReactor.hh"
#include "TExperiment.hh"
#include "TMCToyInterface.hh"
#include "TMacroInterface.hh"

//Please make sure to follow the steps in the order shown below
//TODO: Need to have better error-checking if the execution do not follow these steps

void usage(int arguments)
{
  printf("Cannot use %i arguments\n",arguments);
  printf("Usage:\n");
  printf("$ ./makeToy ./macro.mac\n");
}

int main(int argc,char** argv)
{
  TString inputFile("");
  if(argc==1)inputFile.Append("mac/toyThrow.mac");
  else if(argc==2)inputFile = (argv[1]);
  else usage(argc);
  // Create singleton instance of Macro interface for extracting information from macros
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  
  // This step is very important to be able to use the macro files
  PROSPECTMacroInterface.Initialize(inputFile);
  
  // Name of the output toy file
  TString outputToyFileName("TestingToy.root");
  
  // Output file to save the toy
  TFile* outputToyFile = new TFile(outputToyFileName,"RECREATE");
  
  // Create an experiment object
  TExperiment experiment;
  
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
    printf("For example add 'AddDetector = 101' for adding near detector at near position \n");
    exit(1);
  }
  
  // Setup the experiment
  experiment.SetupExperiment();
  
  // Write the histograms to the output file
  // Useful for comparisons to the toys
  experiment.WriteHistograms(*outputToyFile);
  
  // Create a toy interface
  TMCToyInterface singleToy = TMCToyInterface(experiment);
  
  // Set the seed
  singleToy.SetSeed(0);
  
  // Throw a single toy.  This can be done any number of times.
  // For Feldman Cousins, just make a loop to throw the toy n times.
  singleToy.ThrowToy();
  
  const auto & reactorOnHists=singleToy.GetReactorOnHists();
  const auto & reactorOffHists=singleToy.GetReactorOffHists();
  
  outputToyFile->cd();
  for(auto it:reactorOnHists)
  {
    int detNumber=it.first;
    hRef[detNumber]=(TH2D*)((it.second))->Clone();
    hRef.at(detNumber)->Add(reactorOffHists.at(detNumber),-experiment.fRxOnOffRatio.at(detNumber));
    hRef.at(detNumber)->Write();
//    reactorOffHists.at(detNumber)->Write();
  }
  
  // Write the toy histograms to the file.
  //singleToy.WriteHistograms(*outputToyFile);
  
  return 0;
}

