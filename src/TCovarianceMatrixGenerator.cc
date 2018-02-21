#include "TCovarianceMatrixGenerator.hh"

TCovMatGenerator::TCovMatGenerator(int x)
{
    nBins=x;
    outputFile=new TFile("covariance.root","RECREATE");
    
    energyBins= new double[nBins+1];
    energyBins[0]=0.0;
    for(int i=0;i<nBins-1;i++)
    {
        double tempNumber=1.0+0.2*i;
        energyBins[i+1]=tempNumber;
    }
    energyBins[nBins]=12.0;
}

TCovMatGenerator::TCovMatGenerator(int x, TString outputRoot)
{
    nBins=x;
    outputFile=new TFile(outputRoot,"RECREATE");
}

TCovMatGenerator::TCovMatGenerator(int x, TFile& userDefinedOutput)
{
    nBins=x;
    outputFile=&userDefinedOutput;
}

// RefCovMatGenerator starts here
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
RefCovMatGenerator::RefCovMatGenerator(const vector<double> &v):TCovMatGenerator(v.size())
{
    initialize(v);
}

RefCovMatGenerator::RefCovMatGenerator(const vector<double> &v,TString outputRoot):TCovMatGenerator(v.size(),outputRoot)
{
    initialize(v);
}

RefCovMatGenerator::RefCovMatGenerator(const vector<double> &v,TFile &userDefinedOutput):TCovMatGenerator(v.size(),userDefinedOutput)
{
    initialize(v);
}

int RefCovMatGenerator::initialize(const vector<double>& v)
{
    // Copies contents from v to refVector, might want reduce the cost by not copying
    refVector=v;
    nBins=refVector.size();
    hCovMatrix=new TH2D("CovMatrix","Covariance matrix",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
    hCorrMatrix=new TH2D("CorrMatrix","Correlation matrix",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
    hRedCovMatrix=new TH2D("ReducedCovMatrix","Reduced covariance matrix",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
    hSigmaValues=new TH1D("SigmaValues","Reduced Covariance Matrix; Energy Position Bin; Energy Position Bin",nBins,-0.5,nBins-0.5);
    hMeanFracSample=new TH1D("hMeanFracSample","hMeanFracSample",nBins,0,nBins-0.5);
    
    for(int i=0;i<nHists;i++)
    {
        TString histName = Form("ToyCovMatrix%i",i);
        TH2D* tempHist= new TH2D(histName,"Toy Covariance matrix",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
        hToyCovMatrix.push_back(tempHist);
        tempHist->Reset();
        histName =Form("ToyRedCovMatrix%i",i);
        tempHist= new TH2D(histName,"Toy Covariance matrix",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
        hRedToyCovMatrix.push_back(tempHist);
    }
    return 1;
}

int RefCovMatGenerator::addSample(const vector<double>& sampleVector)
{
    if(sampleVector.size()!=refVector.size())
    {
        std::cerr << "The length of reference vector "<< refVector.size() <<" for sample " << nSamples <<" does not match the length of provided sample/toy " <<
        sampleVector.size()<<std::endl;
        return 0;
    }
    
    for(int i=0;i<nBins;++i)
    {
        if(refVector[i]!=0.0)
        {
            for(int j=0;j<nBins;++j)
            {
                if(refVector[j]!=0.0)
                {
                    double tempCovariance=(sampleVector[i]-refVector[i]);
                    tempCovariance*=(sampleVector[j]-refVector[j]);
                    // PTS: Changes
                    hCovMatrix->Fill(i,j,tempCovariance);
                    if(nSamples<nHists)
                    {
                        hToyCovMatrix[nSamples]->Fill(i,j,tempCovariance);
                        double tempRedCovariance = tempCovariance;
                        double statI= refVector[i];
                        double statJ= refVector[j];
                        if(statI!=0 && statJ!=0)
                        {
                        tempRedCovariance=tempRedCovariance/(statI*statJ);
                        hRedToyCovMatrix[nSamples]->Fill(i,j,tempRedCovariance);
                        }
                    }
                }
            }
        }
        
        double fracSample;
        if(refVector[i]!=0.0)
        {
            fracSample=(sampleVector[i]-refVector[i]);
            fracSample/=refVector[i];
        }
        else fracSample=0.0;
        hMeanFracSample->Fill(i,fracSample);
    }
    
    //Appending a number to nSamples
    ++nSamples;
    return 1;
}

int RefCovMatGenerator::finalize()
{
    hCovMatrix->Scale(1.0/nSamples);
    
    // TODO: Move this to a new utility function
    for(int i = 0; i< nBins ; i++)
    {
        double diagElementI = std::sqrt(hCovMatrix->GetBinContent(i+1,i+1));
        double statI= refVector[i];
        if(diagElementI!=0.0 && statI!=0.0)hSigmaValues->Fill(i, diagElementI/statI);
        for(int j = 0; j< nBins ; j++)
        {
            double diagElementJ = std::sqrt(hCovMatrix->GetBinContent(j+1,j+1));
            double statJ= refVector[j];
            if(diagElementI!=0.0 && diagElementJ!=0.0)
            {
                double tempCorrelation = hCovMatrix->GetBinContent(i+1,j+1);
                tempCorrelation/=(diagElementI*diagElementJ);
                hCorrMatrix->Fill(i,j,tempCorrelation);
            }
            
            if(statI!=0 && statJ!=0)
            {
                double tempRedCovariance = hCovMatrix->GetBinContent(i+1,j+1);
                tempRedCovariance/=(statI*statJ);
                hRedCovMatrix->Fill(i,j,tempRedCovariance);
            }
        }
    }
    hSigmaValues->Write();
    hCovMatrix->Write();
    hCorrMatrix->Write();
    hRedCovMatrix->Write();
    hMeanFracSample->Write();
    // Filling the histograms
    for(int i=0;i<nHists;i++)
    {
        //hToyCovMatrix[i]->Write();
        //hRedToyCovMatrix[i]->Write();
    }
    outputFile->Close();
    return 1;
}

// RefCovMatGenerator ends here
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////


// MeanTCovMatGenerator starts here
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
MeanTCovMatGenerator::MeanTCovMatGenerator(int x):TCovMatGenerator(x)
{
    initialize();
}

MeanTCovMatGenerator::MeanTCovMatGenerator(int x, TString outputRoot):TCovMatGenerator(x,outputRoot)
{
    initialize();
}

int MeanTCovMatGenerator::initialize()
{
    binContent=new double[nBins];
    principal = new TPrincipal(nBins,"D");
    return 1;
}

int MeanTCovMatGenerator::addSample(const vector<double>& sampleVector)
{
    int i=0;
    for(auto it=sampleVector.begin();it!=sampleVector.end();++it)
    {
        binContent[i]=*it;
        ++i;
    }
    principal->AddRow(binContent);
    ++nSamples;
    return 1;
}

int MeanTCovMatGenerator::finalize()
{
    principal->MakePrincipals();
    meanCovMat =principal->GetCovarianceMatrix();
    return 1;
}
