///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: Nov, 2016
// Executable to use toys to generate sensitivity.
////////////////////////////////////////////////////////////////////////

#include "TDetector.hh"
#include "TReactor.hh"
#include "TExperiment.hh"
#include "TOscillationModelBuilder.hh"
#include "TMCToyInterface.hh"
#include "TMacroInterface.hh"
#include "TMinimization.hh"
#include "time.h"

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
  if(argc==1)inputFile.Append("./mac/default.mac");
  else if(argc==2)inputFile = (argv[1]);
  else usage(argc);
  // Create singleton instance of Macro interface for extracting information from macros
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  
  // This step is very important to be able to use the macro files
  PROSPECTMacroInterface.Initialize(inputFile);
  
  // Name of the output toy file
  TString outputFileName("toySensitivity.root");
  
  // Create an experiment object
//  TExperiment experiment;
  
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
  
  std::map<int,TH2D> hRef;
  
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
  
  
  // Read the output minimization file from the macro file
  PROSPECTMacroInterface.RetrieveValue("MinimizationFileName",outputFileName);
  TFile* outputMinimizationFile = new TFile(outputFileName,"RECREATE");
  experiment.WriteHistograms(*outputMinimizationFile);
  
  TOscillationModelBuilder oscillationModelBuilder(experiment);//Create an oscillationModelBuilder simulator object
  oscillationModelBuilder.SetUpModelOscillation();//Setup the OscillationSimulator,takes the default number of bins
  
  oscillationModelBuilder.SimulateModelOscillation();
  
  TMinimization minimization(oscillationModelBuilder);
  minimization.SetupMinimizer();//Setup covariance matrcies and histograms for the minimization process
  
  
//  TString value;
  PROSPECTMacroInterface.RetrieveValue("UseToys",value);
  
  if(value.CompareTo("YES",TString::kIgnoreCase)==0) {
    
    int nToys=1;
    PROSPECTMacroInterface.RetrieveValue("NToys",nToys);
    
    
    // Create a toy interface
    TMCToyInterface toyInterface = TMCToyInterface(experiment);
    
    // Set the seed
    toyInterface.SetSeed(10);
    
    for(int i=0;i<nToys;i++)
    {
      
      outputMinimizationFile->cd();
      toyInterface.ThrowToy();
      const auto & reactorOnHists=toyInterface.GetReactorOnHists();
      const auto & reactorOffHists=toyInterface.GetReactorOffHists();
      
      // Use system time to check running time for each toy
      time_t currentTimer=time(NULL);
      time_t prevTimer=time(NULL);
      
      for(auto it:reactorOnHists)
      {
        int detNumber=it.first;
        (it.second)->Copy((hRef[detNumber]));
        hRef.at(detNumber).Add(reactorOffHists.at(detNumber),-experiment.fRxOnOffRatio.at(detNumber));
        outputMinimizationFile->cd();
        hRef.at(detNumber).Write();
        
        ((minimization.hReferenceLvsE).at(detNumber))=(TH2D*)hRef.at(detNumber).Clone();
      }
      
      minimization.Minimize();//Perform minimization process
      auto minChi=minimization.GetCumulativeChiSquareMap();
      TString chi2MapName;
      chi2MapName.Form("CumulativeChiMap%i",i);
      minimization.SetCumulativeChiSquareMapName(chi2MapName);
      minChi->Write();
      
      time(&currentTimer);
      double dt=difftime(currentTimer,prevTimer);
			if(i%10==0) printf("Running toy %i, dt=%f seconds\n",i,dt);
      time(&prevTimer);
    }
    // Write the toy histograms to the file.
    toyInterface.WriteHistograms(*outputMinimizationFile);
  }
  outputMinimizationFile->Close();
  
  return 0;
}

