#include "TDetectorResponseInterface.hh"

TDetectorResponseInterface::TDetectorResponseInterface()
{
  printf("WARNING: Using the default constructor for TDetectorResponseInterface \n");
  printf("WARNING: Please make sure to call the constructor with detctor code argument \n");
  
  // Call the singleton instance of the Macro interface
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  TString key;
  int detType=detectorCode/100;// Get detector type. i.e, near, middle, far etc
  key.Form("EnergyResolution%i",detType);
  PROSPECTMacroInterface.RetrieveValue(key,fEnergyResolution);// Store energy resolution from macro file
}

TDetectorResponseInterface::TDetectorResponseInterface(int detCode):detectorCode(detCode)
{
  // Call the singleton instance of the Macro interface
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  TString key,value;
  key.Form("DetectorResponseFile%i",detectorCode);
  // If the value exists.
  if(!PROSPECTMacroInterface.RetrieveValue(key,value)){
    printf("WARNING: No input for the key %s exists in the macro file\n",key.Data());
    detRespFileName="./inputs/DetectorResponse.root";
    printf("Trying to use %s file\n",detRespFileName.Data());
  }
  else{
    detRespFileName=value;
  }
  
  detRespFile=TFile::Open(detRespFileName.Data());
  if(detRespFile->IsZombie() || !(detRespFile->IsOpen())){
    printf("The file %s was not loaded properly\n",detRespFileName.Data());
    exit(1);
  }
  isDefaultConstructor=false;
}

TDetectorResponseInterface::~TDetectorResponseInterface()
{
  detRespFile->Close();
  delete detRespFile;
}

/// Apply detector response to the full detector
/*
 This function takes as the input the 1D histogram of a full detector absolute spectrum
 */
void TDetectorResponseInterface::ApplyDetectorResponse(const TH1 &inputHistogram,TH1 *outputHistogram)
{
  if(isDefaultConstructor) detResp.ApplyResponse(inputHistogram,outputHistogram,fEnergyResolution);
  else{
    TString detHistName;
    // TODO: This has to be changed to the correct name of the absolute detector response matrix
    detHistName="Abs";
    TH2D* segDetResponse=(TH2D*)detRespFile->Get(detHistName.Data());
    detResp.ApplyResponse(inputHistogram,outputHistogram,segDetResponse);
  }
}

/// Apply response segment by segment
/*
 This function takes as the input the 2D histogram spectrum with y bins being the baseline bins
 */
void TDetectorResponseInterface::ApplyDetectorResponse(const TH2 &inputHistogram,TH2 *outputHistogram)
{
  if(isDefaultConstructor) detResp.ApplyResponse(inputHistogram,outputHistogram,fEnergyResolution);
  else{
    std::map<int,TH2*> segDetResponse;
    for(int i =0;i<inputHistogram.GetNbinsY();i++){
      
      // Assuming the detector is 14x11, will need to use the TDetector object if this has to be changed.
      int segX = i%14;
      int segY = std::floor(i/(14));
      // Don't bother with the fiducial segments
      if(segX<=0 || segX>=13) continue;
      if(segY<=0 || segY>=10) continue;
      TString detHistName;
      detHistName.Form("Segment%i",i);
      // Assume that the detector response file has detector response matrix histograms even for the non-fiducial segments.
      segDetResponse[i] = (TH2D*)detRespFile->Get(detHistName.Data());// Create the segment response 2D histogram
    }
    detResp.ApplyResponse(inputHistogram,outputHistogram,segDetResponse);
  }
}


void TDetectorResponseInterface::ApplyDetectorResponse(const TH2 &inputHistogram, TH2* outputHistogram, TH1* & outputHistogramFull )
{
  const double neutrino_conversion_energy = 0.800 ; 
  // setup const, bin low edge, bin high edge, number of bins.
	double ENull_xlow = 0 ;
	double ENull_xup  =	10 ;
  double BWidth = 0.100000;
  int ENull_NBin = ( ENull_xup - ENull_xlow ) / BWidth ;


  TH1D* SegENull = new TH1D("SegENull" , "ETrue after energy shift ",  ENull_NBin  , ENull_xlow , ENull_xup  );
  TH1D*    ENull = new TH1D("   ENull" , "ETrue after energy shift ",  ENull_NBin  , ENull_xlow , ENull_xup  );
  // start computing hSegvsENull from hSegvsETrue for each non-fiducial segments.
  // loop through all segments ;
  for (int i = 0; i < inputHistogram.GetNbinsY(); i++ )
    // Assuming the detector is 14x11, will need to use the TDetector object if this has to be changed.
    int segX = i%14;
    int segY = std::floor(i/(14));
    // Don't bother with the fiducial segments
    if(segX<=0 || segX>=13) continue;
    if(segY<=0 || segY>=10) continue;
    // scale to 0 for each segment calculation.
    SegENull->Scale(0);
    for (int j = 1; j <= SegENull->GetNbinsX(); j++ )
    {

    }
    // done with each segment calculation.
    // get full spectrum using the projection function 
    // might need to change the projection range.
    outputHistogramFull = ( TH1D* )SegENull->ProjectionX("",9,40);

}