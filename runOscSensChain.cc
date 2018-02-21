///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Sep 07, 2016
// Executable to run through the sensitivity/oscillation chain .
// This code shows an example that gets the inputs from macro files
////////////////////////////////////////////////////////////////////////

#include <time.h>

#include "TDetector.hh"
#include "TReactor.hh"
#include "TExperiment.hh"
#include "TOscillationModelBuilder.hh"
//#include "TSensitivity.hh"
#include "TMinimization.hh"
//#include "TOscillationFitter.hh"
#include "TMacroInterface.hh"

//Please make sure to follow the steps in the order shown below
//TODO: Need to have better error-checking if the execution do not follow these steps

void usage(int arguments)
{
  printf("Cannot use %i arguments\n",arguments);
  printf("Usage:\n");
  printf("$ ./runOscSensChain ./macro.mac\n");
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
  
  
  double refDeltam2=0.0;
  double refSin22theta=0.0;
  
  TString value;
  PROSPECTMacroInterface.RetrieveValue("IsReferenceOscillated",value);
  if(value.CompareTo("YES",TString::kIgnoreCase)==0) {
    PROSPECTMacroInterface.RetrieveValue("ReferenceDeltam2",refDeltam2);
    PROSPECTMacroInterface.RetrieveValue("ReferenceSin22theta",refSin22theta);
    printf("Using %f, %f for ReferenceDeltam2 and ReferenceSin22theta\n",refDeltam2,refSin22theta);
    if(refSin22theta==0.0 || refDeltam2==0.0){
      printf("--WARNING--\n");
      printf("--Using 0.0, 0.0 for ReferenceDeltam2 and ReferenceSin22theta--\n");
      printf("--Set ReferenceDeltam2 and ReferenceSin22theta in your macro properly\n--");
    }
  }
  TExperiment experiment(refDeltam2,refSin22theta);
  
  
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
      std::cout << "Adding detector " << detNumber <<std::endl;
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
//  experiment.WriteHistograms(*outputMinimizationFile);//Write histograms to the outputSetupFile
  
  TOscillationModelBuilder oscillationModelBuilder(experiment);//Create an oscillationModelBuilder simulator object
  oscillationModelBuilder.SetUpModelOscillation();//Setup the OscillationSimulator,takes the default number of bins
  //TFile* outputOscFile = new TFile(outputOscillationFileName,"RECREATE");//output file with information on the oscillation from the created model
  
  
  TFile* outputMinimizationFile = new TFile(outputMinimizationFileName,"RECREATE"); //
  oscillationModelBuilder.SimulateModelOscillation();//
  oscillationModelBuilder.WriteHistograms(*outputMinimizationFile);//Write histograms to the outputOscFile
  
  /*TString value;
   PROSPECTMacroInterface.RetrieveValue("IsReferenceOscillated",value);
   
   if(value.CompareTo("YES",TString::kIgnoreCase)==0) minimization=new TOscillationFitter(oscillation);
   else minimization= new TSensitivity(oscillation);
   
   //TSensitivity sensitivity(oscillation);
   */
  
  TMinimization minimization(oscillationModelBuilder);
  minimization.SetupMinimizer();//Setup covariance matrcies and histograms for the minimization process
  
  minimization.Minimize();//Perform minimization process
  minimization.WriteHistograms(*outputMinimizationFile);
  
  std::cout << outputMinimizationFile->GetName() <<std::endl;
  //  Close all the files, make sure to close or it might lead to memory leaks
  outputMinimizationFile->Close();
  return 0;
}

