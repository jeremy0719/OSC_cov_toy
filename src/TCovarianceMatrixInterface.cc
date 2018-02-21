#include "TCovarianceMatrixInterface.hh"

TCovarianceMatrixInterface::TCovarianceMatrixInterface(TExperiment &experiment,bool isToyInt):fExperiment(experiment),isToyInterface(isToyInt)
{}

TCovarianceMatrixInterface::~TCovarianceMatrixInterface()
{}

void TCovarianceMatrixInterface::TrackCovMatrixList(int detNumber)
{
  // Call the singleton instance of the Macro interface
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();
  TString key;
  
  // based on what this covaraicnce matrix interface is for, the approriate key has to assigned
  isToyInterface?key.Form("ToyUncertainties%i",detNumber):key.Form("MinimizationUncertainties%i",detNumber);
  
  
  // iterate over the list of the enumerator
  for (int covMatrixList = SSTAT; covMatrixList != SB2B; covMatrixList++ ) {
    // cast to CovarianceMatrices enum
    CovarianceMatrices covMatrixType = static_cast<CovarianceMatrices>(covMatrixList);
    TString covType = GetString(covMatrixType); // Assign the type of covariance matrix to string
    // Check from the macro file if this covariance matrix should be used and if so store in the bitset
    if(PROSPECTMacroInterface.CheckListForValue(key,covType))CovMatList[detNumber].set(covMatrixType);
  }
}

void TCovarianceMatrixInterface::ExtractCovarianceMatrix(int detNumber,int covFileType,TH2D* hRCovMatrix)
{
  CovarianceMatrices covMatrixType = static_cast<CovarianceMatrices>(covFileType);    // cast to CovarianceMatrices enum
  TString covType = GetString(covMatrixType);// Assign the type of covariance matrix to string
  //  covType=(covType.BeginsWith("T")?covType.Strip(TString::kLeading,'T'):covType.Strip(TString::kLeading,'S'));
  
  
  TMacroInterface& PROSPECTMacroInterface = TMacroInterface::Instance();  // Call the singleton instance of the Macro interface
  TString key="CovFile";  // Initialize key for now append and prepend other terms later
  TString value;
  TString iSSignal=(covType.BeginsWith("S")?"Signal":"Background");// Check if the uncertainty is Signal or Background
  covType.Remove(0,1); // Remove the leading term in covType string
  key.Prepend(covType);// Prepend the type of cov matrix to key
  isToyInterface?key.Prepend("Toy"):key.Prepend(iSSignal); // prepend toy or signal based on the inerface type
  key.Append(std::to_string(detNumber)); //Append detector number
  key.ToUpper();// Convert the whole sting to upper level/ done just to maintain consistency
  // Making sure that the code flows without break even if the macro file doesn't dfine what covariance matrices are to be used.
  if(CovMatList.find(detNumber) == CovMatList.end()) return;
  // Check if covariance file type has a matrix associated with it
  if(CovMatList.at(detNumber)[covFileType]){
    if(!PROSPECTMacroInterface.RetrieveValue(key,value)){
      double signalNormSigma;
      double bkgNormSigma;
      double signalShapeSigma;
      double bkgShapeSigma;
      // HARD-CODED NUMBERS
      // For toys, the uncertainty is obtained from the RAA deficit value, whereas for minimization it has to be a 100% deficit
      signalNormSigma=isToyInterface?0.06:1;
      bkgNormSigma=isToyInterface?0:0.02;
      signalShapeSigma=isToyInterface?0.1:0.1;
      bkgShapeSigma=isToyInterface?0:0.01;
      // If the covariance matrix file is not defined but the uncertainty is signal/bkg norm, use hard-coded uncertainties instead and prepare covariance matrices
      if(covFileType>=SNORM && covFileType<=BSHAPE){
        hTempMat=(TH2D*)hRCovMatrix->Clone("ReducedCovMatrix");
        hTempMat->Scale(0.0);
        switch (covFileType) {
          case SNORM:
          {
            for(int i=1;i<=hTempMat->GetXaxis()->GetNbins();i++){
              for(int j=1;j<=hTempMat->GetYaxis()->GetNbins();j++){
                hTempMat->SetBinContent(i,j,signalNormSigma*signalNormSigma);
              }
            }
          }
            break;
            
          case BNORM:
          {
            for(int i=1;i<=hTempMat->GetXaxis()->GetNbins();i++){
              for(int j=1;j<=hTempMat->GetYaxis()->GetNbins();j++){
                hTempMat->SetBinContent(i,j,bkgNormSigma*bkgNormSigma);
              }
            }
          }
            break;
            
          case SSHAPE:
          {
            for(int i=1;i<=hTempMat->GetXaxis()->GetNbins();i++){
              for(int j=1;j<=hTempMat->GetYaxis()->GetNbins();j++){
                double tempSigmaSquare=signalShapeSigma*signalShapeSigma;
                // Fill in the uncertainties only in the bins corresponding to the same energy bin for each position bin
                if(!((std::abs(i-j))%fEneBins))hTempMat->SetBinContent(i,j,tempSigmaSquare);
              }
            }
          }
            break;
            
          case BSHAPE:
          {
            for(int i=1;i<=hTempMat->GetXaxis()->GetNbins();i++){
              for(int j=1;j<=hTempMat->GetYaxis()->GetNbins();j++){
                double tempSigmaSquare=bkgShapeSigma*bkgShapeSigma;
                if(!((std::abs(i-j))%fEneBins))hTempMat->SetBinContent(i,j,tempSigmaSquare);
              }
            }
          }
            break;
            
          default:
            break;
        }
      }
      else displayMacroError(key);
    }
    else{
      TFile *tempFile= new TFile(value.Data());// Open the specified file
      if(!tempFile->IsOpen() || tempFile->IsZombie()){
        displayFileError(value);
      }
      hTempMat=(TH2D*)tempFile->Get("ReducedCovMatrix");// Get the temporary histogram
    }
    hRCovMatrix->Add(hTempMat,1);// Add the temporary histogram to specified covariance matrix
  }
}


void TCovarianceMatrixInterface::ExtractScalingVector(){
  TVectorD TotScaleVector;
  ExtractInvertedVector(fExperiment.hLvsENull,TotScaleVector);
  TotSignalScalingMatrix.ResizeTo(TotScaleVector.GetNoElements(),TotScaleVector.GetNoElements());
  TotSignalScalingMatrix=OuterProduct(TotScaleVector,TotScaleVector);
}

void TCovarianceMatrixInterface::ExtractScalingVectors(int detNumber){
  TVectorD DetScaleVector;
  ExtractInvertedVector(fExperiment.hLvsENull.at(detNumber),DetScaleVector);
  SignalScalingMatrix[detNumber].ResizeTo(DetScaleVector.GetNoElements(),DetScaleVector.GetNoElements());
  SignalScalingMatrix.at(detNumber)=OuterProduct(DetScaleVector,DetScaleVector);
}


void TCovarianceMatrixInterface::ExtractCovarianceMatrices(int detNumber)
{
  TString histName;
  /// Initialize all the covariance histograms before filling them
  histName.Form("hRCovUCSigMatrix%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hRCovUCSigMatrix[detNumber]=new TH2D(histName.Data(),"hRCovUCSigMatrix",fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("hRCovFCSigMatrix%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hRCovFCSigMatrix[detNumber]=new TH2D(histName.Data(),"hRCovFCSigMatrix",fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("hRCovUCBkgMatrix%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hRCovUCBkgMatrix[detNumber]=new TH2D(histName.Data(),"hRCovUCBkgMatrix",fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("hRCovFCBkgMatrix%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hRCovFCBkgMatrix[detNumber]=new TH2D(histName.Data(),"hRCovFCBkgMatrix",fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  // add signal normalization matrix to the signal FC covariance matrix
  ExtractCovarianceMatrix(detNumber,SNORM,hRCovFCSigMatrix[detNumber]);
  // add background normalization matrix to the bkg FC covariance matrix
  ExtractCovarianceMatrix(detNumber,BNORM,hRCovFCBkgMatrix[detNumber]);
  // If/when the backgrounds are actually established at all the locations, this will be UC
  // add signal shape matrix to the signal FC covariance matrix
  ExtractCovarianceMatrix(detNumber,SSHAPE,hRCovFCSigMatrix[detNumber]);
  // add bkg shape matrix to the bkg UC covariance matrix
  ExtractCovarianceMatrix(detNumber,BSHAPE,hRCovFCSigMatrix[detNumber]);
  // add signal Energy scale matrix to the signal UC covariance matrix
  ExtractCovarianceMatrix(detNumber,SESCALE,hRCovUCSigMatrix[detNumber]);
  // add signal bin to bin matrix to the signal UC covariance matrix
  ExtractCovarianceMatrix(detNumber,SB2B,hRCovUCSigMatrix[detNumber]);
  
  /*-----This list will increase, need to implement a better way if it gets out of hand***/
  
  // Create the Covariance Matrices needed.
  // FC Matrices are Fully Correlated
  // UC Matrices are UnCorrelated
  // PC Matrices are Partially Correlated, needed in future
  fCovMatrixSigFC[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixSigUC[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixBkgOnFC[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixBkgOnUC[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixBkgOffFC[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixBkgOffUC[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixSigStat[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixBkgOnStat[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixBkgOffStat[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixBkgStat[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixTot[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixTotRelative[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
  fCovMatrixTotHist[detNumber].ResizeTo(fDetCovMatBins.at(detNumber), fDetCovMatBins.at(detNumber));
}


void TCovarianceMatrixInterface::SetupCovarianceMatrixHists(int detNumber)
{
  TString histName;
  // Final Covariance matrices used in the analysis.
  histName.Form("SigCovMatrixFC%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hSigCovMatrixFC[detNumber] = new TH2D(histName, "Fully Correlated Signal Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("SigCovMatrixUC%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hSigCovMatrixUC[detNumber] = new TH2D(histName, "Uncorrelated Signal Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("BkgOnCovMatrixFC%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hBkgOnCovMatrixFC[detNumber] = new TH2D(histName, "Fully Correlated Reactor-On Background Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("BkgOffCovMatrixFC%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hBkgOffCovMatrixFC[detNumber] = new TH2D(histName, "Fully Correlated Reactor-Off Background Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("BkgOnCovMatrixUC%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hBkgOnCovMatrixUC[detNumber] = new TH2D(histName, "Uncorrelated Reactor-On Background Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("BkgOffCovMatrixUC%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hBkgOffCovMatrixUC[detNumber] = new TH2D(histName, "Uncorrelated Reator-Off Background Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("SigStatCovMatrix%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hSigStatCovMatrix[detNumber] = new TH2D(histName, "Full Statistical Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("BkgOnStatCovMatrix%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hBkgOnStatCovMatrix[detNumber] = new TH2D(histName, "Background Reactor-On Statistical Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("BkgOffStatCovMatrix%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hBkgOffStatCovMatrix[detNumber] = new TH2D(histName, "Background Reactor-Off Statistical Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("BkgStatCovMatrix%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hBkgStatCovMatrix[detNumber] = new TH2D(histName, "Background Statistical Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("TotalCovMatrix%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hTotalCovMatrix[detNumber] = new TH2D(histName, "Full Total Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
  
  histName.Form("TotalCovMatrixRelative%i",detNumber);
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hTotalCovMatrixRelative[detNumber] = new TH2D(histName, "Full Total Covariance Matrix", fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5,fDetCovMatBins.at(detNumber), -0.5, fDetCovMatBins.at(detNumber)-0.5);
}

void TCovarianceMatrixInterface::SetupDetCovarianceMatrices(int detNumber)
{
  SetupCovarianceMatrixHists(detNumber);
  TrackCovMatrixList(detNumber);
  ExtractCovarianceMatrices(detNumber);
  ExtractScalingVectors(detNumber);
}

void TCovarianceMatrixInterface::SetupCovarianceMatrices(int eneBins)
{
  fEneBins=eneBins;
  TString histName;
  // Loop through the detectors
  for(auto it:fExperiment.fDetectors){
    int detNumber = it->fDetectorCode;
    // Call setup covariance matrix for each individual detector
    SetupDetCovarianceMatrices(detNumber);
  }
  ExtractScalingVector();
  // Cumulative covariance matrix
  fCumulativeCovMatrix.ResizeTo(fCumulativeCovMatBins, fCumulativeCovMatBins);
  fCumulativeCovMatrixRelative.ResizeTo(fCumulativeCovMatBins, fCumulativeCovMatBins);
  fCumulativeCovMatrixHist.ResizeTo(fCumulativeCovMatBins, fCumulativeCovMatBins);
  
  fReactorOnSigCovMatrix.ResizeTo(fCumulativeCovMatBins,fCumulativeCovMatBins);
  fReactorOnBkgCovMatrix.ResizeTo(fCumulativeCovMatBins,fCumulativeCovMatBins);
  fReactorOffBkgCovMatrix.ResizeTo(fCumulativeCovMatBins,fCumulativeCovMatBins);
  
  histName = "CumulativeCovMatrix";
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hCumulativeCovMatrix = new TH2D(histName.Data(),"",fCumulativeCovMatBins, -0.5, fCumulativeCovMatBins-0.5,fCumulativeCovMatBins, -0.5, fCumulativeCovMatBins-0.5);
  histName = "CumulativeCovMatrixRelative";
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hCumulativeCovMatrixRelative = new TH2D(histName.Data(),"",fCumulativeCovMatBins, -0.5, fCumulativeCovMatBins-0.5,fCumulativeCovMatBins, -0.5, fCumulativeCovMatBins-0.5);
  
  histName = "hReactorOnSigCovMatrix";
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hReactorOnSigCovMatrix = new TH2D(histName.Data(),"",fCumulativeCovMatBins, -0.5, fCumulativeCovMatBins-0.5,fCumulativeCovMatBins, -0.5, fCumulativeCovMatBins-0.5);
  histName = "hReactorOnBkgCovMatrix";
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hReactorOnBkgCovMatrix = new TH2D(histName.Data(),"",fCumulativeCovMatBins, -0.5, fCumulativeCovMatBins-0.5,fCumulativeCovMatBins, -0.5, fCumulativeCovMatBins-0.5);
  histName = "hReactorOffBkgCovMatrix";
  isToyInterface?histName.Prepend("Toy"):histName.Prepend("Min");
  hReactorOffBkgCovMatrix = new TH2D(histName.Data(),"",fCumulativeCovMatBins, -0.5, fCumulativeCovMatBins-0.5,fCumulativeCovMatBins, -0.5, fCumulativeCovMatBins-0.5);
}

void TCovarianceMatrixInterface::BuildCovarianceMatrices()
{
  
  int detCount=0;//keep track of number of detectors parsed through
  int detCovMatBins =0;//keep track of number of toal bins parsed through
  
  printf("Building covariance matrices\n");
  for(auto it:fExperiment.fDetectors){
    
    int detNumber = it->fDetectorCode;
    for (int i = 0; i < fDetCovMatBins.at(detNumber); i++) {
      
      // Access the position bin and energy bin indices.
      int posIndexI = TMath::Floor((double) (i) / fEneBins);
      int eneIndexI = (i) - posIndexI * fEneBins;
      
      // Get the total statistics in bin I Reactor On
      double predSigRxOnI = fExperiment.hLvsENull.at(detNumber)->GetBinContent(eneIndexI+1, posIndexI+1);
      double predBkgRxOnI = fExperiment.hBkgRxOnLvsE.at(detNumber)->GetBinContent(eneIndexI+1, posIndexI+1);
      double predTotRxOnI = predSigRxOnI + predBkgRxOnI;
      
      // Get the total statistics in bin I Reactor Off
      double predBkgRxOffI = fExperiment.hBkgRxOffLvsE.at(detNumber)->GetBinContent(eneIndexI+1, posIndexI+1);
      
      for (int j = 0; j < fDetCovMatBins.at(detNumber); j++) {
        
        // Access the position bin and energy bin indices.
        int posIndexJ = TMath::Floor((double) (j) / fEneBins);
        int eneIndexJ = (j) - posIndexJ * fEneBins;
        
        // Get the total statistics in bin J Reactor On
        double predSigRxOnJ = fExperiment.hLvsENull.at(detNumber)->GetBinContent(eneIndexJ+1, posIndexJ+1);
        double predBkgRxOnJ = fExperiment.hBkgRxOnLvsE.at(detNumber)->GetBinContent(eneIndexJ+1, posIndexJ+1);
        double predTotRxOnJ = predSigRxOnJ + predBkgRxOnJ;
        
        // Get the total statistics in bin J Reactor Off
        double predBkgRxOffJ = fExperiment.hBkgRxOffLvsE.at(detNumber)->GetBinContent(eneIndexJ+1, posIndexJ+1);
        
        // Build up stats matrix msking sure that the stat matrix for this detector has to be built
        // The error on S+B is Sqrt(S+B).
        
        fCovMatrixSigStat.at(detNumber)(i, j) = TMath::Sqrt(predTotRxOnI*predTotRxOnJ) * (i==j);
        // Build up signal systematic matrix
        // Need to convert the reduced covariance matrix back into a full covariance matrix
        // Uncorrelated Errors
        fCovMatrixSigUC.at(detNumber)(i, j) = hRCovUCSigMatrix.at(detNumber)->GetBinContent(i+1, j+1) * predSigRxOnI * predSigRxOnJ;
        
        // Correlated Errors
        fCovMatrixSigFC.at(detNumber)(i, j) = hRCovFCSigMatrix.at(detNumber)->GetBinContent(eneIndexI+1, eneIndexJ+1) * predSigRxOnI * predSigRxOnJ;
        
        // Build up background systematic matrix
        // Need to convert the reduced covariance matrix back into a full covariance matrix
        // Uncorrelated Errors
        fCovMatrixBkgOnUC.at(detNumber)(i, j) = hRCovUCBkgMatrix.at(detNumber)->GetBinContent(i+1, j+1) * predBkgRxOnI * predBkgRxOnJ;
        fCovMatrixBkgOffUC.at(detNumber)(i, j) = hRCovUCBkgMatrix.at(detNumber)->GetBinContent(i+1, j+1) * predBkgRxOffI * predBkgRxOffJ;
        
        // Covariance matrix corresponding to the reactor-on and reactor-off background statistics
        // These matrices will be used for toy geneation
        fCovMatrixBkgOnStat.at(detNumber)(i, j) = TMath::Sqrt(predBkgRxOnI * predBkgRxOnJ) * (i==j);
        fCovMatrixBkgOffStat.at(detNumber)(i, j) = TMath::Sqrt(predBkgRxOffI * predBkgRxOffI) * (i==j);
        
        // Since the background is "measured" when the reactor is off,
        // the statistical error on the background should be determined
        // by the statistics of the reactor off period.
        fCovMatrixBkgStat.at(detNumber)(i, j) = TMath::Power(fExperiment.fRxOnOffRatio.at(detNumber),2)*fCovMatrixBkgOffStat.at(detNumber)(i, j);
        
        // Correlated Errors
        fCovMatrixBkgOnFC.at(detNumber)(i, j) = hRCovFCBkgMatrix.at(detNumber)->GetBinContent(eneIndexI+1, eneIndexJ+1) * predBkgRxOnI * predBkgRxOnJ;
        fCovMatrixBkgOffFC.at(detNumber)(i, j) = hRCovFCBkgMatrix.at(detNumber)->GetBinContent(eneIndexI+1, eneIndexJ+1) * predBkgRxOffI * predBkgRxOffJ;
        
      } // End of loop through j bins.
    } // End of loop through i bins.
    
    // Add Matrices
    // Options: CovMatrixSigStat, CovMatrixBkgStat, CovMatrixSigFC, CovMatrixSigUC, CovMatrixBkgFC, CovMatrixBkgUC
    fCovMatrixTot.at(detNumber) = fCovMatrixSigStat.at(detNumber) + fCovMatrixBkgStat.at(detNumber) + fCovMatrixSigFC.at(detNumber)+fCovMatrixBkgOnFC.at(detNumber) + fCovMatrixSigUC.at(detNumber) +fCovMatrixBkgOnFC.at(detNumber);
    
    // Loop through bins of total covariance matrix and fill the associated histogram.
    // Do this before inversion.
    for (int i = 0; i < fDetCovMatBins.at(detNumber); i++) {
      for (int j = 0; j < fDetCovMatBins.at(detNumber); j++) {
        fCumulativeCovMatrix(detCovMatBins+i, detCovMatBins+j) = fCovMatrixTot.at(detNumber)(i,j);
        (fReactorOnSigCovMatrix)(detCovMatBins+i, detCovMatBins+j) = fCovMatrixSigStat.at(detNumber)(i,j) + fCovMatrixSigFC.at(detNumber)(i,j) + fCovMatrixSigUC.at(detNumber)(i,j);
        (fReactorOnBkgCovMatrix)(detCovMatBins+i, detCovMatBins+j) = fCovMatrixBkgOnStat.at(detNumber)(i,j) + fCovMatrixBkgOnFC.at(detNumber)(i,j) + fCovMatrixBkgOnUC.at(detNumber)(i,j) ;
        (fReactorOffBkgCovMatrix)(detCovMatBins+i, detCovMatBins+j) = fCovMatrixBkgOffStat.at(detNumber)(i,j)  + fCovMatrixBkgOffFC.at(detNumber)(i,j)+ fCovMatrixBkgOffUC.at(detNumber)(i,j);
      }
    }
    fCovMatrixTotRelative.at(detNumber)=fCovMatrixTot.at(detNumber);
    fCovMatrixTotHist.at(detNumber)=fCovMatrixTot.at(detNumber);
    ElementMult(fCovMatrixTotRelative.at(detNumber),SignalScalingMatrix.at(detNumber));
    printf("fCovMatrixTot Covariance Matrix for the detector %i constructed.\n",it->fDetectorCode );
    // Invert matrix
    if(fCovMatrixTot.at(detNumber).Invert()==0){
      printf("fCovMatrixTot Matrix couldn't be inverted\n");
      exit(1);
    }
    if(fCovMatrixTotRelative.at(detNumber).Invert()==0){
      printf("fCovMatrixTot Matrix couldn't be inverted\n");
      exit(1);
    }
    printf("fCovMatrixTot Covariance matrix for the detector %i inverted \n",it->fDetectorCode );
    
    detCount += 1;
    detCovMatBins += fDetCovMatBins.at(detNumber);
  }
  
  
  // The block matrices corresponding to each detector are filled, but the detector level correlations have not been filled and that has to be done at this point
  for (int i = 0; i < fCumulativeCovMatBins; i++){
    // Access the detector index
    // It is important to see that fDetectorType is same between same detector at different locations
    int detNumberI;
    int maxCovBinsI=0;
    int minCovBinsI=0;
    
    for(auto it:fExperiment.fDetectors){
      detNumberI = it->fDetectorCode;
      maxCovBinsI+= fDetCovMatBins.at(detNumberI);
      if(maxCovBinsI>=i) break;
      minCovBinsI += fDetCovMatBins.at(detNumberI);
    }
    // Access the position bin and energy bin indices.
    int posIndexI = TMath::Floor((double) (i -minCovBinsI) / fEneBins);
    int eneIndexI = (i - minCovBinsI) - posIndexI * fEneBins;
    double predSigRxOnI = fExperiment.hLvsENull.at(detNumberI)->GetBinContent(eneIndexI+1, posIndexI+1);
    double predBkgRxOnI = fExperiment.hBkgRxOnLvsE.at(detNumberI)->GetBinContent(eneIndexI+1, posIndexI+1);
    double predBkgRxOffI = fExperiment.hBkgRxOffLvsE.at(detNumberI)->GetBinContent(eneIndexI+1, posIndexI+1);
    
    for (int j = 0; j < fCumulativeCovMatBins; j++){
      
      int detNumberJ;
      int maxCovBinsJ=0;
      int minCovBinsJ=0;
      
      // Access the detector index
      for(auto it:fExperiment.fDetectors){
        detNumberJ = it->fDetectorCode;
        maxCovBinsJ+= fDetCovMatBins.at(detNumberJ);
        if(maxCovBinsJ>=j) break;
        minCovBinsJ += fDetCovMatBins.at(detNumberJ);
      }
      // Access the position bin and energy bin indices.
      int posIndexJ = TMath::Floor((double) (j-minCovBinsJ) / fEneBins);
      int eneIndexJ = (j-minCovBinsJ) - posIndexJ * fEneBins;
      double predSigRxOnJ = fExperiment.hLvsENull.at(detNumberJ)->GetBinContent(eneIndexJ+1, posIndexJ+1);
      double predBkgRxOnJ = fExperiment.hBkgRxOnLvsE.at(detNumberJ)->GetBinContent(eneIndexJ+1, posIndexJ+1);
      double predBkgRxOffJ = fExperiment.hBkgRxOffLvsE.at(detNumberJ)->GetBinContent(eneIndexJ+1, posIndexJ+1);
      
      if(detNumberI==detNumberJ) continue;
      // Only true of matrices fully correlated b/w detectors
      // because these values have already been added
      fCumulativeCovMatrix(i, j) = hRCovFCSigMatrix.at(detNumberJ)->GetBinContent(eneIndexI+1, eneIndexJ+1) * predSigRxOnI * predSigRxOnJ;
      (fReactorOnSigCovMatrix)(i,j) = hRCovFCSigMatrix.at(detNumberJ)->GetBinContent(eneIndexI+1, eneIndexJ+1) * predSigRxOnI * predSigRxOnJ;
      (fReactorOnBkgCovMatrix)(i,j) = hRCovFCBkgMatrix.at(detNumberJ)->GetBinContent(eneIndexI+1, eneIndexJ+1) * predBkgRxOnI * predBkgRxOnJ;
      (fReactorOffBkgCovMatrix)(i,j) = hRCovFCBkgMatrix.at(detNumberJ)->GetBinContent(eneIndexI+1, eneIndexJ+1) * predBkgRxOffI * predBkgRxOffJ;
    }
  }
  fCumulativeCovMatrixHist=fCumulativeCovMatrix;
  fCumulativeCovMatrixRelative=fCumulativeCovMatrix;
  ElementMult(fCumulativeCovMatrixRelative,TotSignalScalingMatrix);
  printf("Cumulative Covariance Matrix constructed.\n");
  
  if(fCumulativeCovMatrix.Invert()==0){
    printf("fCumulativeCovMatrix Matrix couldn't be inverted\n");
    exit(1);
  }
  if(fCumulativeCovMatrixRelative.Invert()==0){
    printf("fCumulativeCovMatrixRelative Matrix couldn't be inverted\n");
    exit(1);
  }
  printf("Cumulative Covariance Matrix inverted.\n");
}

void TCovarianceMatrixInterface::FillCovHistograms()
{
  int detCount=0;//keep track of number of detectors parsed through
  int detCovMatBins =0;//keep track of number of toal bins parsed through
  
  for(auto it:fExperiment.fDetectors){
    
    int detNumber = it->fDetectorCode;
    for (int i = 0; i < fDetCovMatBins.at(detNumber); i++) {
      for (int j = 0; j < fDetCovMatBins.at(detNumber); j++) {
        hSigStatCovMatrix.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixSigStat.at(detNumber)(i, j));
        hSigCovMatrixUC.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixSigUC.at(detNumber)(i, j));
        hSigCovMatrixFC.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixSigFC.at(detNumber)(i, j));
        hBkgOnCovMatrixUC.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixBkgOnUC.at(detNumber)(i, j));
        hBkgOffCovMatrixUC.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixBkgOffUC.at(detNumber)(i, j));
        hBkgOnStatCovMatrix.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixBkgOnStat.at(detNumber)(i, j));
        hBkgOffStatCovMatrix.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixBkgOffStat.at(detNumber)(i, j));
        hBkgStatCovMatrix.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixBkgStat.at(detNumber)(i, j));
        hBkgOnCovMatrixFC.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixBkgOnFC.at(detNumber)(i, j));
        hBkgOffCovMatrixFC.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixBkgOffFC.at(detNumber)(i, j));
      }
    }
    for (int i = 0; i < fDetCovMatBins.at(detNumber); i++) {
      for (int j = 0; j < fDetCovMatBins.at(detNumber); j++) {
        //        hCumulativeCovMatrix->SetBinContent(detCovMatBins+i+1, detCovMatBins+j+1, fCumulativeCovMatrix(detCovMatBins+i, detCovMatBins+j));
        hTotalCovMatrix.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixTotHist.at(detNumber)(i,j));
        hTotalCovMatrixRelative.at(detNumber)->SetBinContent(i+1, j+1, fCovMatrixTotRelative.at(detNumber)(i,j));
      }
    }
    printf("Covariance matrix histograms for the detctor %i filled.\n",it->fDetectorCode );
    detCount += 1;
    detCovMatBins += fDetCovMatBins.at(detNumber);
  }
  for (int i = 0; i < fCumulativeCovMatBins; i++){
    for (int j = 0; j < fCumulativeCovMatBins; j++){
      hCumulativeCovMatrix->SetBinContent(i+1, j+1, fCumulativeCovMatrixHist(i, j));
      hCumulativeCovMatrixRelative->SetBinContent(i+1, j+1, fCumulativeCovMatrixRelative(i, j));
      hReactorOnSigCovMatrix->SetBinContent(i+1, j+1, (fReactorOnSigCovMatrix)(i, j));
      hReactorOnBkgCovMatrix->SetBinContent(i+1, j+1,  (fReactorOnBkgCovMatrix)(i, j));
      hReactorOffBkgCovMatrix->SetBinContent(i+1, j+1, (fReactorOffBkgCovMatrix)(i, j));
    }
  }
  printf("Cumulative covariance matrix histograms filled.\n");
}


void TCovarianceMatrixInterface::WriteHistograms(TFile& outputFile)
{
  FillCovHistograms();
  printf("Writing Covariance Matrix histograms \n");
  if(!outputFile.IsOpen() || outputFile.IsZombie()){
    printf("File not open\n");
    return;
  }
  outputFile.cd();
  for(auto it:fExperiment.fDetectors){
    int detNumber = it->fDetectorCode;
    // Write all the covariance matrix 2D histograms
    hSigCovMatrixFC.at(detNumber)->Write();
    hSigCovMatrixUC.at(detNumber)->Write();
    hBkgOnCovMatrixFC.at(detNumber)->Write();
    hBkgOffCovMatrixFC.at(detNumber)->Write();
    hBkgOnCovMatrixUC.at(detNumber)->Write();
    hBkgOffCovMatrixUC.at(detNumber)->Write();
    hSigStatCovMatrix.at(detNumber)->Write();
    hBkgOnStatCovMatrix.at(detNumber)->Write();
    hBkgOffStatCovMatrix.at(detNumber)->Write();
    hBkgStatCovMatrix.at(detNumber)->Write();
    hTotalCovMatrix.at(detNumber)->Write();
    hTotalCovMatrixRelative.at(detNumber)->Write();
  }
  hCumulativeCovMatrix->Write();
  hCumulativeCovMatrixRelative->Write();
  hReactorOnSigCovMatrix->Write();
  hReactorOnBkgCovMatrix->Write();
  hReactorOffBkgCovMatrix->Write();
}
