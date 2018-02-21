#include "TDetectorLocation.hh"

void TDetectorLocation::ReadMacroInputs()
{
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  TString key;
  
  // Set the detector-reactor distances from the macro
  // Mostly will not need to be changed, but it is good to have the ability to set
  
  int detectorCode=fDetectorType*100+fDetectorPosition;
  key.Form("Xdistance%i",detectorCode);
  PROSPECTMacroInterface.RetrieveValue(key,fDistX);
  key.Form("Ydistance%i",detectorCode);
  PROSPECTMacroInterface.RetrieveValue(key,fDistY);
  key.Form("Zdistance%i",detectorCode);
  PROSPECTMacroInterface.RetrieveValue(key,fDistZ);
  key.Form("AngleOffset%i",detectorCode);
  PROSPECTMacroInterface.RetrieveValue(key,fAngleOffset);
}

TDetectorLocation::TDetectorLocation():fDetectorType(0),fDetectorPosition(0)
{
}

void TDetectorLocation::ComputeDistances()
{
  // Make sure the detector is defined
  if(fDetectorType==0||fDetectorPosition==0){
    printf("You need to assign detector type and detector position before computing distances\n");
    printf("Either call the non-default constructor or call DefineDetectorLocation method\n");
    exit(1);
  }
  if(fDetectorType==1){
    // The distance from the reactor to the XY center of the near detector with the Z coordinate corresponding to the floor (front position)
    // Originally based on reactor center to detector center (14x10) and 14.6 cm segments engineering drawing.
    // Edited on Jan 11, 2017 based on docdb.wlab.yale.edu/0015/001528/001/Final_Shield.pdf
    fDistX = 5.904;
    fDistY = 0.0;
    fDistZ = 5.16;
    
    ReadMacroInputs();
    
    fDistX = fDistX + (fDetectorPosition-1)*1.5;
    
    // Doesn't need fixing here
    // Set the angular offset of the detector
    // fAngleOffset = 0.;
  }
  
  else if(fDetectorType==2){
    // The distance from the reactor to the XY center of the far detector with the Z coordinate correstponding to the floor
    fDistX = 15.8 + 14 * 0.146;
    fDistY = 0.0;
    fDistZ = 0.9;
    
    ReadMacroInputs();
    
    // Doesn't need fixing here
    // Set the angular offset of the detector
    // fAngleOffset = 0.;
  }
}

TDetectorLocation::TDetectorLocation(int detType, int detPosition):fDetectorType(detType),fDetectorPosition(detPosition)
{
  ComputeDistances();
}

void TDetectorLocation::DefineDetectorLocation(int detType, int detPosition){
  // Only define a type and location if it hasn't been defined yet
  if(!(fDetectorType==0)||!(fDetectorPosition==0)){
    printf("Cannot reassign the detector type and position\n");
    exit(1);
  }
  fDetectorType = detType;
  fDetectorPosition = detPosition;
  ComputeDistances();
}

void TDetectorLocation::Print()
{
  // Make sure the detector is defined
  if(fDetectorType==0||fDetectorPosition==0){
    printf("You need to assign detector type and detector position\n");
    printf("Either call the non-default constructor or call DefineDetectorLocation method\n");
    exit(1);
  }
  printf("PROSPECT Detector Location : Type is %i\n",fDetectorType);
  printf("PROSPECT Detector Location : Position identifier is %i\n",fDetectorPosition);
  printf("PROSPECT Detector Location : Distance to first segment in x, y and z are %f, %f and %f respectively\n",fDistX,fDistY,fDistZ);
  printf("PROSPECT Detector Location : Detector angular offset is %f\n",fAngleOffset);
}
