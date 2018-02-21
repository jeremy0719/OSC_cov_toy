#include "TDetectorResponse.hh"
#include "TF1.h"
#include "TMath.h"

// Smearing from Xianyi, 2017/08/30
void TDetectorResponse::ApplyResponse(const TH1 &inputHist,TH1* outputHist, double res){
//  if (!inputHist) return;
  
  int nBins = inputHist.GetNbinsX();
  
  for (int i = 0; i < nBins; i++){
    double x = inputHist.GetBinCenter(i+1);	// The X value of the bin to be smeared
    double sigma = res*TMath::Sqrt(x);	// Using a fixed resolution. Change "x" to "TMath::Sqrt(x)" to apply energy dependent resolution.
    int bandMin=1;
    int bandMax=inputHist.GetNbinsX();
    double thisBinContent = inputHist.GetBinContent(i+1);	// Get the content of this bin.
				double thisBinWidth	= inputHist.GetBinWidth(i+1);	// Get the width of the this bin.
				assert(thisBinWidth == inputHist.GetBinWidth(i+1));	//Bin widths of inputHist and inputHist must equal.
    
				// Now, loop through the bins in 5*sigma range to fill contents.
    for (int j = bandMin; j < bandMax+1; j++){
      double E = inputHist.GetBinCenter(j);
						double lowE = inputHist.GetBinLowEdge(j);
						double highE = lowE + thisBinWidth;
      double lowEGaus = thisBinContent*TMath::Gaus(lowE, x, sigma, false)/TMath::Sqrt(2.*TMath::Pi())/sigma;
						double centerGaus = thisBinContent*TMath::Gaus(E, x, sigma, false)/TMath::Sqrt(2.*TMath::Pi())/sigma;
      double highEGaus = thisBinContent*TMath::Gaus(highE, x, sigma, false)/TMath::Sqrt(2.*TMath::Pi())/sigma;
						
						double spreadBinContent = 0;	// The content needs to be filled.
						// If not the center bin, fill single trapezoid integral.
						//if (j != i+1) spreadBinContent = (lowEGaus+highEGaus)*thisBinWidth/2;
      
						// If is center bin but 2*sigma is wider than half bin, fill two trapezoid on each side of the center.
      if (thisBinWidth < 4*sigma) spreadBinContent = (lowEGaus+centerGaus)*(thisBinWidth/2)/2 + (highEGaus+centerGaus)*(thisBinWidth/2)/2;
      
						// If is center bin and half bin is wider than 2*sigma, fill the trapezoid value on each sigma.
      if (thisBinWidth > 4*sigma) {
        int k = 1;
        double prevE, currE, prevEGaus, currEGaus;
        while (thisBinWidth > 2*k*sigma){
          if (k > 5) break; // If the bin width > 5*sigma, integrate up to 5*sigma.
          prevE = E - sigma*(k-1);
          currE = E - sigma*k;
          prevEGaus = thisBinContent*TMath::Gaus(prevE, x, sigma, false)/TMath::Sqrt(2.*TMath::Pi())/sigma;
          currEGaus = thisBinContent*TMath::Gaus(currE, x, sigma, false)/TMath::Sqrt(2.*TMath::Pi())/sigma;
          spreadBinContent += 2*(prevEGaus+currEGaus)*(sigma)/2;
          k++;
        }
        if (k <= 5) spreadBinContent += 2*(lowEGaus+currEGaus)*(currE-lowE)/2;
        // There is a simplified way, instead of using the while loop:
        // currE = E - sigma*2;
        // currEGaus = thisBinContent*TMath::Gaus(currE, x, sigma, false)/TMath::Sqrt(2.*TMath::Pi())/sigma;
        // speadBinContent = 2*(lowEGaus+currEGaus)*(currE-lowE)/2 + 2*(centerGaus+currEGaus)*(sigma)/2;
      }
      
      outputHist->AddBinContent(j, spreadBinContent);
    }
  }
}

void TDetectorResponse::ApplyResponse(const TH2 &inputHist,TH2* outputHist, double smearing){
  //
  //  // temporary inputHist for copying input inputHist to be smeared
  //  TH2D *tempResHist = (TH2D*)inputHist.Clone("tempResHist");
  //  inputHist.Scale(0.0);
  
  // Project 2D histogram onto X-axis to 1D histogram
  TH1D *tempProjectedHist;
  TH1D *tempOutProjectedHist;
  
  // Loop through the y bins
  for(int j=1; j<=inputHist.GetNbinsY();j++){
    // Make sure to empty the histogram
    TString histName;
    histName.Form("%i",j);// Set a unique name for the histogram
    // Project the histogram corresponsing to the particular y bin
    tempProjectedHist = (TH1D*)inputHist.ProjectionX(histName,j,j);
    histName.Form("out%i",j);// Set a unique name for the histogram
    tempOutProjectedHist = (TH1D*)inputHist.ProjectionX(histName);
    tempOutProjectedHist->Scale(0);
    ApplyResponse(*tempProjectedHist,tempOutProjectedHist,smearing);// Apply smearing to the projected histogram
    for(int i=1; i<=inputHist.GetNbinsX();i++){
      // Fill the input histogram with the smeared content
      outputHist->SetBinContent(i,j,tempProjectedHist->GetBinContent(i));
    }
    tempProjectedHist->Scale(0.0);
    tempProjectedHist->SetDirectory(0);// Decouple histogram from the directory
  }
  delete tempOutProjectedHist;
  delete tempProjectedHist;
}

void TDetectorResponse::ApplyResponse(const TH1 &inputHist, TH1* outputHist, TH2* detResponseMatrix){
  // temporary inputHist for copying input inputHist to the response histogram
  //  TH1D *tempResHist = (TH1D*)inputHist.Clone("tempResHist");
  //Check if the input histogram and the detector response matrix have the same number of bins
  if(inputHist.GetNbinsX()!=detResponseMatrix->GetNbinsX()){
    printf("WARNING:The input histogram and the response matrix do not have the same number of bins\n");
    printf("WARNING:Not able to apply detector response\n");
    return;
  }
  outputHist->Scale(0);
  // Run through detector response matrix true energy bins
  for(int i =1;i<=detResponseMatrix->GetNbinsX();i++){
    //Figure out true energy bin content
    double trueBinContent=inputHist.GetBinContent(i);
    // Run through the reconstructed bin contents
    for(int j =1;j<=detResponseMatrix->GetNbinsY();j++){
      double iJContribution=0.0;
      // Contribution of true bin to the reconstructed bin content
      iJContribution=detResponseMatrix->GetBinContent(i,j);
      // Fill in the spectrum with detector response applied
      outputHist->Fill(detResponseMatrix->GetYaxis()->GetBinCenter(j),iJContribution*trueBinContent);
    }
  }
  // Delete temporary histogram
  //  delete tempResHist;
}

void TDetectorResponse::ApplyResponse(const TH2 &inputHist,TH2* outputHist,TH2* detResponseMatrix){
  
  // temporary inputHist for copying input inputHist to be smeared
  //  TH2D *tempResHist = (TH2D*)inputHist.Clone("tempResHist");
  outputHist->Scale(0);
  
  if((inputHist.GetNbinsX()!=detResponseMatrix->GetNbinsX())||(detResponseMatrix->GetNbinsX()!=detResponseMatrix->GetNbinsY())){
    printf("There is a binnning issue in applying detector response matrix\n");
    printf("Please make sure that the number of energy bins are consistent b/w input histogram and detector response matrix\n");
    printf("Also make sure that the detector response matrix is a square matrix\n");
    exit(1);
  }
  // Project 2D histogram onto X-axis to 1D histogram
  TH1D *tempProjectedHist;
  TH1D *tempOutProjectedHist;
  
  // Loop through the y bins
  for(int j=1; j<=inputHist.GetNbinsY();j++){
    // Make sure to empty the histogram
    TString histName;
    histName.Form("%i",j);// Set a unique name for the histogram
    // Project the histogram corresponsing to the particular y bin
    tempProjectedHist = (TH1D*)inputHist.ProjectionX(histName,j,j);
    histName.Form("out%i",j+1);// Set a unique name for the histogram
    tempOutProjectedHist = (TH1D*)detResponseMatrix->ProjectionY(histName);
    tempOutProjectedHist->Scale(0);
    ApplyResponse(*tempProjectedHist,tempOutProjectedHist,detResponseMatrix);// Apply smearing to the projected histogram
    for(int i=1; i<=inputHist.GetNbinsX();i++){
      // Fill the input histogram with the smeared content
      outputHist->SetBinContent(i,j,tempOutProjectedHist->GetBinContent(i));
    }
    tempProjectedHist->Scale(0.0);
    tempProjectedHist->SetDirectory(0);// Decouple histogram from the directory
  }
  delete tempProjectedHist;
  delete tempOutProjectedHist;
}



void TDetectorResponse::ApplyResponse(const TH2 &inputHist,TH2* outputHist,const std::map<int,TH2*> &detResponseMatrices){
  // temporary inputHist for copying input inputHist to be smeared
  //  TH2D *tempResHist = (TH2D*)inputHist.Clone("tempResHist");
  //  inputHist.Scale(0.0);Keep in mind the true energy have have wider range than reconstructed energy
  outputHist->Scale(0.0);
  // Project 2D histogram onto X-axis to 1D histogram
  TH1D *tempProjectedHist;
  TH1D *tempOutProjectedHist;
  
  // Loop through the y bins
  for(int j=0; j<inputHist.GetNbinsY();j++){
    
    // Assuming the detector is 14x11, will need to use the TDetector object if this has to be changed.
    int segX = j%14;
    int segY = std::floor(j/(14));
    // Don't bother with the fiducial segments
    if(segX<=0 || segX>=13) continue;
    if(segY<=0 || segY>=10) continue;
    
    
    if((inputHist.GetNbinsX()!=detResponseMatrices.at(j)->GetNbinsX())){
      printf("There is a binnning issue in applying detector response matrix\n");
      printf("Please make sure that the number of energy bins are consistent b/w input histogram and detector response matrices\n");
      exit(1);
    }
    
    TString histName;
    histName.Form("%i",j+1);// Set a unique name for the histogram
    // Project the histogram corresponsing to the particular y bin
    tempProjectedHist = (TH1D*)inputHist.ProjectionX(histName,j+1,j+1);
    histName.Form("out%i",j+1);// Set a unique name for the histogram
    tempOutProjectedHist = (TH1D*)detResponseMatrices.at(j)->ProjectionY(histName);
    tempOutProjectedHist->Scale(0);
    ApplyResponse(*tempProjectedHist,tempOutProjectedHist,detResponseMatrices.at(j));// Apply response to the projected histogram
    
    for(int i=1; i<=tempOutProjectedHist->GetNbinsX();i++){
      // Fill the input histogram with the smeared content
      outputHist->Fill(tempOutProjectedHist->GetBinCenter(i),inputHist.GetYaxis()->GetBinCenter(j+1),tempOutProjectedHist->GetBinContent(i));
    }
    tempProjectedHist->Scale(0.0);
    tempProjectedHist->Scale(0.0);
    tempOutProjectedHist->SetDirectory(0);// Decouple histogram from the directory
    tempOutProjectedHist->SetDirectory(0);// Decouple histogram from the directory
  }
  delete tempProjectedHist;
  delete tempOutProjectedHist;
  //  delete tempResHist;
}


