#include "TCovUtilities.hh"

TH2D* rebinTH2(TH2* inputHistogram)
{
    int nBins=inputHistogram->GetXaxis()->GetNbins();
    double* energyBins= new double[nBins+1];
    
    energyBins[0]=0.0;
    for(int i=0;i<nBins-1;i++)
    {
        double tempNumber=1.0+0.2*i;
        energyBins[i+1]=tempNumber;
    }
    energyBins[nBins]=12.0;
    
    TH2D *h = new TH2D("hMeanCovMatrixRebin",inputHistogram->GetTitle(),nBins,energyBins,nBins,energyBins);
    TAxis *xaxis = inputHistogram->GetXaxis();
    TAxis *yaxis = inputHistogram->GetYaxis();
    for (int j=1; j<=yaxis->GetNbins();j++)
    {
        for (int i=1; i<=xaxis->GetNbins();i++)
        {
            h->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),inputHistogram->GetBinContent(i,j));
        }
    }
    delete [] energyBins;
    return h;
}

/// Functions borrowed from Karin Gilje

// This function takes an input of a number of different throws in each pull.
std::vector<TH1D*> CreateThrowHistograms(int nThrowTypes, double sigma) {
    std::vector<TH1D*> outputVector;
    outputVector.clear();
    
    for (int i = 0; i < nThrowTypes; i++) {
        TString title = Form("Throw %d; Throw (%%)", i+1);
        TString name = Form("Throw%d", i+1);
        TH1D* tempThrows = new TH1D(name, title, 50, -5.0*sigma, 5.0*sigma);
        outputVector.push_back(tempThrows);
    }
    
    return outputVector;
}

// This function takes a vector of histograms and clones their structure.
std::vector<TH1D*> CreateToyHistograms(std::vector<TH1D*> inputHistograms) {
    std::vector<TH1D*> outputVector;
    outputVector.clear();
    
    for (unsigned int i = 0; i < inputHistograms.size(); i++) {
        TString name = Form("tempBin%d", i+1);
        TH1D* tempHist = (TH1D*)inputHistograms[i]->Clone(name);
        tempHist->Scale(0.0);
        outputVector.push_back(tempHist);
    }
    
    return outputVector;
}
