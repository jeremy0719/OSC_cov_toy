#include "TMinimization.hh"

TMinimization::TMinimization(TOscillationModelBuilder & oscModelBuilder):fOscillationModelBuilder(oscModelBuilder),fExperiment(oscModelBuilder.fExperiment),  fCovMatInt(fExperiment,false)
{
  TMacroInterface &PROSPECTMacroInterface = TMacroInterface::Instance();
  TString key,value;
  key.Form("doRelativeMinimization");
  PROSPECTMacroInterface.RetrieveValue(key,value);
  doRelativeMinimization=(value.CompareTo("YES",TString::kIgnoreCase)==0)?true:false;
  /*TODO:: Throw exception when the reactorOnFraction + reactorOffFraction > 1.0, cannot exit(1) in a constructor
   if (reactorOnFraction + reactorOffFraction > 1.0) {
   printf("The fraction of Reactor On Time and Cosmic Reactor Off Time is more than 1\n");
   printf("This is impossible ;) Returning... \n");
   }*/
}

TMinimization::~TMinimization()
{
  //  delete hCumulativeChiSquareMap;
}

void TMinimization::SetupMinimizer()
{
  // Obtain the number of bins in deltam2 and sin22theta
  fNDeltam2 = fOscillationModelBuilder.fNDeltam2;
  fNSinSq2Theta = fOscillationModelBuilder.fNSinSq2Theta;
  
  TString histName = "CumulativeChiSquareMap";
  hCumulativeChiSquareMap = new TH2D(histName.Data(),"#chi^{2} Values; sin^{2} 2#theta_{14}; #Delta m^{2}_{14}", fNSinSq2Theta, fOscillationModelBuilder.fSinSq2ThetaBins, fNDeltam2, fOscillationModelBuilder.fDeltam2Bins);
  
  //    hLvsENull.at(detNumber)->Scale();
  //    hLvsEBackgroundRxOn.at(detNumber)->Scale();
  //    hLvsEBackgroundRxOff.at(detNumber)->Scale();
  
  
  for(auto it:fExperiment.fDetectors){
    int detNumber = it->fDetectorCode;
    histName.Form("LvsE_%i",detNumber);
    hLvsENull[detNumber] = (TH2D*)(fExperiment.hLvsENull.at(detNumber)->Clone(histName));
    // Calculate the number of energy bins (the same for every position!)
    if(fEneBins==0) fEneBins = hLvsENull.at(detNumber)->GetNbinsX();
    double posBins = hLvsENull.at(detNumber)->GetNbinsY();
    fCovMatInt.fDetCovMatBins[detNumber] = fEneBins * posBins;
    fCovMatInt.fCumulativeCovMatBins += fCovMatInt.fDetCovMatBins.at(detNumber);
  }
  
  fCovMatInt.SetupCovarianceMatrices(fEneBins);
  
  for(auto it:fExperiment.fDetectors){
    int detNumber = it->fDetectorCode;
    
    for (int j = 0; j < fNDeltam2; j++) {
      histName.Form("LvsE_%i_%i",j,detNumber);
      int deltam2Number = 1000*j + detNumber;
      hLvsEDeltam2[deltam2Number] = (TH2D*)(fOscillationModelBuilder.hLvsEDeltam2.at(deltam2Number)->Clone(histName));
      //      LvsETemp->Scale();
      //      hLvsEDeltam2.push_back(LvsETemp);
    }
    if(doRelativeMinimization){
      histName.Form("LvsERefRelative_%i",detNumber);
      hReferenceLvsE[detNumber]=(TH2D*)fExperiment.hLvsERelativeRef.at(detNumber)->Clone(histName);
    }
    else{
      histName.Form("LvsERef_%i",detNumber);
      hReferenceLvsE[detNumber]=(TH2D*)fExperiment.hLvsERef.at(detNumber)->Clone(histName);
    }
    
    histName.Form("ChiSquareMap%i",detNumber);
    hChiSquareMap[detNumber] = new TH2D(histName,"#chi^{2} Values; sin^{2} 2#theta_{14}; #Delta m^{2}_{14}", fNSinSq2Theta, fOscillationModelBuilder.fSinSq2ThetaBins, fNDeltam2, fOscillationModelBuilder.fDeltam2Bins);
    
    printf("The total number of covariance bins for detector %i is: %i\n",detNumber, fCovMatInt.fDetCovMatBins.at(detNumber));
  } //End of loop through detectors
  
  fCovMatInt.BuildCovarianceMatrices();
}

double TMinimization::CalculateChiSquareValue(TDetector& detector, double deltam2, double sinSq2Theta)
{
  int detNumber= detector.fDetectorCode;
  TVectorD XVector(fCovMatInt.fDetCovMatBins.at(detNumber));
  
  
  TH2D OscHistTemp;
  hLvsENull.at(detNumber)->Copy(OscHistTemp);
  TH2D lvsEDeltam2;
  hLvsENull.at(detNumber)->Copy(lvsEDeltam2);
  lvsEDeltam2.Scale(0.0);
  fOscillationModelBuilder.SimulateModelOscillation(fExperiment,detector,lvsEDeltam2,deltam2);
  
  OscHistTemp.Add(&lvsEDeltam2, -1.0 * sinSq2Theta);
  
  TH2D OscHistTempRelative;
  // 1D histogram, projection of temporary oscillated scaled histogram
  TH1D OscTotHistTemp=*(TH1D*)OscHistTemp.ProjectionX("",1,OscHistTemp.GetNbinsY());
  OscHistTemp.Copy(OscHistTempRelative);
  RelativizeHistogram(*hLvsENull.at(detNumber),*fExperiment.hENull.at(detNumber),OscHistTemp,OscTotHistTemp,OscHistTempRelative);
  OscHistTempRelative.Copy(OscHistTemp);
  
  for (int i = 0; i < fCovMatInt.fDetCovMatBins.at(detNumber); i++) {
    // Access the position bin and energy bin indices.
    int posIndexI = TMath::Floor((double) i / fEneBins);
    int eneIndexI = i - posIndexI * fEneBins;
    
    // Get the total statistics in bin I
    double obsSigRxOnI= hReferenceLvsE.at(detNumber)->GetBinContent(eneIndexI+1, posIndexI+1);
    double  predSigRxOnI = OscHistTemp.GetBinContent(eneIndexI+1, posIndexI+1);
    
    // Fill X vector
    XVector(i) = obsSigRxOnI-predSigRxOnI;
    
  } // End of loop through i bins.
  
  //Check and make sure the covariance matrix and the vector have the same number of terms
  if(fCovMatInt.fCovMatrixTot.at(detNumber).GetNrows()!=XVector.GetNrows()){
    printf("The covariance matrix size is %i where it has to be %i\n",
           fCovMatInt.fCovMatrixTot.at(detNumber).GetNrows(),XVector.GetNrows());
    printf("Exiting\n");
    exit(1);
  }
  
  // Find Transpose
  TVectorD XTVector = XVector;
  
  // Perform X^T * C^{-1} * X
  if(doRelativeMinimization)XTVector *= fCovMatInt.fCovMatrixTotRelative.at(detNumber);
  else XTVector *= fCovMatInt.fCovMatrixTot.at(detNumber);
  double chiSquareValue = XTVector * XVector;
  return chiSquareValue;
}


void TMinimization::BuildChiSquareMap()
{
  //  hCumulativeChiSquareMap->Reset();
  //  for(auto it:fExperiment.fDetectors)hChiSquareMap.at(it->fDetectorCode)->Reset();
  
  // Loop through the delta m values.
  for (int m = 0; m < fNDeltam2; m++) {
    
    //printf("Evaluating Delta m2 = %f\n",fOscillationModelBuilder.fDeltam2[m] );
    
    // Loop through sin2 2theta values.
    for (int s = 0; s < fNSinSq2Theta; s++) {
      
      TVectorD XCumulativeVector(fCovMatInt.fCumulativeCovMatBins);
      int detCovMatBins=0;
      
      for(auto it:fExperiment.fDetectors){
        int detNumber = it->fDetectorCode;
        // Create the "X" vector
        TVectorD XVector(fCovMatInt.fDetCovMatBins.at(detNumber));
        
        // Build temporary histogram
        // The oscillation histogram corresponding to the oscillation parameters m and s
        TH2D OscHistTemp; // Temporary oscillation histogram
        hLvsENull.at(detNumber)->Copy(OscHistTemp);// Copy LvsENull to this hist
        int deltam2Number = 1000*m + detNumber;
        // Subtract the LvsEDeltam2 histograms
        OscHistTemp.Add(hLvsEDeltam2.at(deltam2Number), -1.0 * fOscillationModelBuilder.fSinSq2Theta[s]);
        
        // Only enter this loop if the LvsE histogram has to be scaled.
        if(doRelativeMinimization){
          OscHistTemp.Copy(OscHistTemp);
          // 1D histogram, projection of temporary oscillated scaled histogram
          TH1D OscTotHistTemp=*(TH1D*)OscHistTemp.ProjectionX("",1,OscHistTemp.GetNbinsY());
          TH2D OscHistTempRelative;
          OscHistTemp.Copy(OscHistTempRelative);
          RelativizeHistogram(*hLvsENull.at(detNumber),*fExperiment.hENull.at(detNumber),OscHistTemp,OscTotHistTemp,OscHistTempRelative);
          OscHistTemp.Scale(0.0);
          OscHistTempRelative.Copy(OscHistTemp);
        }
        
        for (int i = 0; i < fCovMatInt.fDetCovMatBins.at(detNumber); i++) {
          // Access the position bin and energy bin indices.
          int posIndexI = TMath::Floor((double) i / fEneBins);
          int eneIndexI = i - posIndexI * fEneBins;
          
          // Get the total statistics in bin I
          double obsSigRxOnI = hReferenceLvsE.at(detNumber)->GetBinContent(eneIndexI+1, posIndexI+1);
          double predSigRxOnI  = OscHistTemp.GetBinContent(eneIndexI+1, posIndexI+1);
          
          // Fill X vector
          XVector(i) = obsSigRxOnI-predSigRxOnI;
          XCumulativeVector(detCovMatBins+i) = obsSigRxOnI-predSigRxOnI;
          
        } // End of loop through i bins.
        // Find Transpose
        TVectorD XTVector = XVector;
        
        //Check and make sure the covariance matrix and the vector have the same number of terms
        if(fCovMatInt.fCovMatrixTot.at(detNumber).GetNrows()!=XVector.GetNrows()){
          printf("The covariance matrix size is %i where it has to be %i\n",
                 fCovMatInt.fCovMatrixTot.at(detNumber).GetNrows(),XVector.GetNrows());
          printf("Exiting\n");
          exit(1);
        }
        // Perform X^T * C^{-1} * X
        if(doRelativeMinimization){
          XTVector *= fCovMatInt.fCovMatrixTotRelative.at(detNumber);
        }
        else {
          XTVector *= fCovMatInt.fCovMatrixTot.at(detNumber);
        }
        double CovarianceChiSquare = XTVector * XVector;
        // Fill Chi Square Histogram
        hChiSquareMap.at(detNumber)->Fill(fOscillationModelBuilder.fSinSq2Theta[s], fOscillationModelBuilder.fDeltam2[m], CovarianceChiSquare);
        //hChiSquareMap.at(detNumber)->SetBinContent(s,m, CovarianceChiSquare);
        detCovMatBins += fCovMatInt.fDetCovMatBins.at(detNumber);
      }
      //Check and make sure the covariance matrix and the vector have the same number of terms
      if(fCovMatInt.fCumulativeCovMatrix.GetNrows()!=XCumulativeVector.GetNrows()){
        printf("The covariance matrix size is %i where it has to be %i\n",
               fCovMatInt.fCumulativeCovMatrix.GetNrows(),XCumulativeVector.GetNrows());
        printf("Exiting\n");
        exit(1);
      }
      TVectorD XTCumulativeVector = XCumulativeVector;
      if(doRelativeMinimization) XTCumulativeVector *= fCovMatInt.fCumulativeCovMatrixRelative;
      else XTCumulativeVector *= fCovMatInt.fCumulativeCovMatrix;
      double CovarianceChiSquare = XTCumulativeVector * XCumulativeVector;
      hCumulativeChiSquareMap->Fill(fOscillationModelBuilder.fSinSq2Theta[s], fOscillationModelBuilder.fDeltam2[m], CovarianceChiSquare);
      //      hCumulativeChiSquareMap->SetBinContent(s,m, CovarianceChiSquare);
    } // End of loop through sin22theta values.
  } // End of loop through delta m2 values.
}


void TMinimization::Minimize()
{
  // Clear histograms before applying minimization(mainly applicable for toys)
  ClearChi2Histograms();
  BuildChiSquareMap();
}

void TMinimization::ClearChi2Histograms()
{
  for(auto it:fExperiment.fDetectors){
    int detNumber = it->fDetectorCode;
    hChiSquareMap.at(detNumber)->Reset();}
  hCumulativeChiSquareMap->Reset();
}

void TMinimization::WriteHistograms(TFile& outputFile)
{
  printf("Writing minimization histograms \n");
  if(!outputFile.IsOpen() || outputFile.IsZombie()){
    printf("File not open\n");
    return;
  }
  outputFile.cd();
  for(auto it:fExperiment.fDetectors){
    int detNumber = it->fDetectorCode;
    fCovMatInt.hCumulativeCovMatrix->Write();
    hChiSquareMap.at(detNumber)->Write();
    hCumulativeChiSquareMap->Write();
    
    hLvsENull.at(detNumber)->Write();
    hReferenceLvsE.at(detNumber)->Write();
    
    fExperiment.hENull.at(detNumber)->Write();
    fExperiment.hETrue.at(detNumber)->Write();
    fExperiment.hSegvsENull.at(detNumber)->Write();
    fExperiment.hSegvsETrue.at(detNumber)->Write();
    fExperiment.hLvsERelative.at(detNumber)->Write();
    
    fExperiment.hERef.at(detNumber)->Write();
    fExperiment.hSegvsERef.at(detNumber)->Write();
    fExperiment.hLvsERef.at(detNumber)->Write();
    fExperiment.hLvsERelativeRef.at(detNumber)->Write();
    
//    fExperiment.hEfficiency->Write();
//    fExperiment.hRxOn->Write();
//    fExperiment.hRxOff->Write();
    fExperiment.hBkgRxOnLvsE.at(detNumber)->Write();
    fExperiment.hBkgRxOffLvsE.at(detNumber)->Write();
  }
  fCovMatInt.WriteHistograms(outputFile);
}

