///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Sep 07, 2016
// Executable to estimate sensitivity of the PROSPECT experiment.
// This code shows an example that gets the inputs from macro files
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
  printf("$ ./estimateSensitivity ./macro.mac\n");
}

int main(int argc,char** argv)
{
  TString inputFile("");
  if(argc==1)inputFile.Append("mac/default.mac");
  else if(argc==2)inputFile = (argv[1]);
  else usage(argc);
  // Create singleton instance of Macro interface for extracting information from macros
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  
  // This step is very important to be able to use the macro files
  PROSPECTMacroInterface.Initialize(inputFile);
  
  // Name of the output setup file
  //TString outputSetupFileName("set.root");
  // Name of the output oscillation file
  TString outputOscillationFileName("Osc.root");
  // Name of the output minimization file
  TString outputMinimizationFileName("Min.root");
  
  // Read the output setup file from the macro file
  //PROSPECTMacroInterface.RetrieveValue("SetupFileName",outputSetupFileName);
  // Read the output oscillation file from the macro file
  PROSPECTMacroInterface.RetrieveValue("OscillationFileName",outputOscillationFileName);
  // Read the output minimization file from the macro file
  PROSPECTMacroInterface.RetrieveValue("MinimizationFileName",outputMinimizationFileName);
  
  //TFile* outputSetupFile = new TFile(outputSetupFileName,"RECREATE");//output file for the experimental setup and event and position trees
  TExperiment experiment;// Create an experiment object
  
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
  
  experiment.SetupExperiment();
  
  //experiment.WriteTrees(*outputSetupFile);//Write position and event trees to the outputSetupFile
  //experiment.WriteHistograms(*outputSetupFile);//Write histograms to the outputSetupFile
  
  TOscillationModelBuilder oscillationModelBuilder(experiment);//Create an oscillation simulator object
  oscillationModelBuilder.SetUpModelOscillation();//Setup the OscillationSimulator,takes the default number of bins
  //TFile* outputOscFile = new TFile(outputOscillationFileName,"RECREATE");//output file with information on the oscillation from the created model

  oscillationModelBuilder.SimulateModelOscillation();//
  //oscillation.WriteOscHistograms(*outputOscFile);//Write histograms to the outputOscFile
  
  TFile* outputMinimizationFile = new TFile(outputMinimizationFileName,"RECREATE"); //
  
  // Sensitivity object
  TMinimization *minimization=new TMinimization(oscillationModelBuilder);
  
  //TSensitivity sensitivity(oscillation);
  minimization->SetupMinimizer();//Setup covariance matrcies and histograms for the minimization process
  minimization->Minimize();//Perform minimization process
  minimization->WriteHistograms(*outputMinimizationFile);//Write histograms to the outputMinimizationFile

  //  Close all the files, make sure to close or it might lead to memory leaks
  //outputSetupFile->Close();
  //outputOscFile->Close();
  outputMinimizationFile->Close();
  return 0;
}

