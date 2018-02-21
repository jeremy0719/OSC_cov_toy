#include "TOscSensUtils.hh"

int Get2DHistBins(const TH2D *inputHist){
  int totalBins=inputHist->GetSize();
  totalBins-=2*(inputHist->GetNbinsX()+inputHist->GetNbinsY())+4;
  return totalBins;
}

// Turn a TH2D to a vector
void ExtractVector(TH2D *inputHist,TVectorD& outputVector) {
  
  int nBinsTot = Get2DHistBins(inputHist);
  outputVector.ResizeTo(nBinsTot);
  
  int iElement=0;
  
  for (int j = 0; j < inputHist->GetNbinsY(); j++) {
    for (int i = 0; i < inputHist->GetNbinsX(); i++) {
      outputVector(iElement) = inputHist->GetBinContent(i+1, j+1);
      iElement++;
    } // End of loop through Y bins
  } // End of loop through X bins
}

// Turn map of TH2Ds to a vector
void ExtractVector(std::map<int, TH2D*> inputMap,TVectorD& outputVector) {
  
  int nBinsTot = 0;
  
  // Loop through the map to get the size of the output vector
  std::map<int, TH2D*>::iterator it;
  for (it = inputMap.begin(); it != inputMap.end(); it++) {
    nBinsTot += Get2DHistBins(it->second);
    
  } // End of loop through detectors
  
  outputVector.ResizeTo(nBinsTot);
  
  // Index of position in output vector
  int iElement = 0;
  
  // Loop through the map again, filling the vector as you go.
  for (it = inputMap.begin(); it != inputMap.end(); it++) {
    
    // Calculate the number of bins.
    int nBinsX = (it->second)->GetNbinsX();
    int nBinsY = (it->second)->GetNbinsY();
    
    for (int j = 0; j < nBinsY; j++) {
      for (int i = 0; i < nBinsX; i++) {
        outputVector(iElement) = (it->second)->GetBinContent(i+1, j+1);
        iElement++;
      } // End of loop through Y bins
    } // End of loop through X bins
    
  } // End of loop through detectors
}


// Turn TH2D to a vector but with inverted elements
void ExtractInvertedVector(TH2D* inputHist,TVectorD &outputVector) {
  TVectorD directVector;
  ExtractVector(inputHist,directVector);
  outputVector.ResizeTo(directVector.GetNoElements());
  int i=0;
  while(i<directVector.GetNoElements()){
    if(directVector(i)==0){
      printf("Something wrong, there are elemets with values '0'\n");
      exit(1);
    }
    else outputVector(i)=1/directVector(i);
    ++i;
  }
}


// Turn map of TH2Ds to a vector but with inverted elements
void ExtractInvertedVector(std::map<int, TH2D*> inputMap, TVectorD &outputVector) {
  int prevBins=0;
  int totBins=0;
  std::map<int, TH2D*>::iterator it;
  // Loop through the map again, filling the vector as you go.
  for (it = inputMap.begin(); it != inputMap.end(); it++) {
    TVectorD detInvVector;
    ExtractInvertedVector(it->second,detInvVector);
    int tempBins=detInvVector.GetNoElements();
    totBins+=tempBins;
    outputVector.ResizeTo(totBins);
    int i=0;
    while(i<tempBins){
      outputVector(prevBins+i)=detInvVector(i);
      ++i;
    }
    prevBins+=tempBins;
  }
}


void RelativizeHistogram(const TH2D& hLvsE, const TH1D& hE,const TH2D& hLvsEInput, const TH1D& hEInput, TH2D& hLvsEOutput){
//  hLvsEInput.Copy(hLvsEOutput);
  hLvsEOutput.Scale(0);
  for(int j=1;j<=hLvsEInput.GetNbinsY();j++){
    for (int i = 1; i <= hLvsEInput.GetNbinsX(); i++) {
      // Get scaling factor, this term essentially tells us what the bin contents of an
      // absolute reference histogram scaled to each baseline would look like
      double scalingFactor=hEInput.GetBinContent(i)*hLvsE.GetBinContent(i,j)/hE.GetBinContent(i);
      // This is essentially LvsE(L)/E
      double scaledBinContent=hLvsEInput.GetBinContent(i,j)/scalingFactor;
      hLvsEOutput.SetBinContent(i,j,scaledBinContent);
    }
  }
}

bool CheckConsistency(const TH1 &h1,const TH1 &h2){
  if(h1.GetNbinsX()!=h2.GetNbinsX()) return false;
  if(h1.GetXaxis()->GetBinLowEdge(1)!=h2.GetXaxis()->GetBinLowEdge(1)) return false;
  if(h1.GetXaxis()->GetBinLowEdge(h1.GetNbinsX())!=h2.GetXaxis()->GetBinLowEdge(h1.GetNbinsX())) return false;
  for(int i=1;i<=h1.GetNbinsX();i++){
    if(h1.GetXaxis()->GetBinCenter(i)!=h1.GetXaxis()->GetBinCenter(i)) return false;
  }
  return true;
}

bool CheckConsistency(const TH2 &h1,const TH2 &h2){
  if(h1.GetNbinsY()!=h2.GetNbinsY()) return false;
  if(h1.GetNbinsY()!=h2.GetNbinsY()) return false;
  if(h1.GetYaxis()->GetBinLowEdge(1)!=h2.GetYaxis()->GetBinLowEdge(1)) return false;
  if(h1.GetYaxis()->GetBinLowEdge(h1.GetNbinsY())!=h2.GetYaxis()->GetBinLowEdge(h1.GetNbinsY())) return false;
  for(int i=1;i<=h1.GetNbinsY();i++){
    if(h1.GetYaxis()->GetBinCenter(i)!=h1.GetYaxis()->GetBinCenter(i)) return false;
  }
  return true;
}
