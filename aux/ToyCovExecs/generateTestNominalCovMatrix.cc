// Program to generate a covariance matrix for normalization of

#include "TCovarianceMatrixGenerator.hh"
#include <sstream>
// Program defaults to this point if the input is not correct
int Usage()
{
  // The input root file should be the output from generateSetupFiles
  printf("Make sure to input the reference ROOT file and detector number \n");
  printf("Usage: ./generateSignalNormCovMatrix detector_number\n");
  printf("Usage: ./generateSignalNormCovMatrix detector_number nToys\n");
  return -1;
}


int main(int argc, char** argv)
{
  if(argc!=2 && argc!=3 )
  {
    int errorCode = Usage();
    std::cerr<< "Exiting with error code " << errorCode <<std::endl;
    return errorCode;
  }
  
  TFile* inputFile = new TFile("Nu_Spectrum.root");
  int nToys;
  // If number of toys are input use it, otherwise use 10000 toys
  if(argc==2) {
    nToys = 10000;
  }
  else nToys = atoi(argv[2]);
  
  
  printf("Using %i toys\n",nToys);
  // Define random nymber generator
  TRandom3* random=new TRandom3();
  TString inputFileName = inputFile->GetName();
  // Check if file is attached.
  if (!inputFile->IsOpen() || inputFile->IsZombie())
  {
    std::cerr << "File not attached or object doesn't exist" << std::endl;
    std::cerr << "Returning ..." << std::endl;
    return -1;
  }
  std::clog << "Accessing " << inputFileName <<std::endl;
  
  // Opening default LvsE plot
  TString hName("hObsSpectrum");
  TKey *key = inputFile->FindKey(hName);
  if (key ==0)
  {
    printf("Histogram %s does not exist!!\n",hName.Data());
    return -1;
  }
  
  TH1D* hObsSpectrum = (TH1D*)inputFile->Get(hName);
  
  // Create Output File
  TString outputName("Cov_Norm.root");
  TFile* outputFile = new TFile(outputName,"RECREATE");
  
  double nEneBins = hObsSpectrum->GetXaxis()->GetNbins();
  double covMatBins = nEneBins;
  
  TH1D* EneHist = new TH1D("EneHist", "Energy Histograms; Energy bin number; arb", covMatBins, -0.5, covMatBins-0.5);
  
  vector <double> aInput(covMatBins);
  //Initial set of values for reference
  
  for(int i=0;i<nEneBins;i++)
  {
    double tempBinContent=0.0;
    tempBinContent = hObsSpectrum->GetBinContent(i+1);
    aInput[i]= tempBinContent;
    EneHist->SetBinContent(i,tempBinContent);
  }
  hObsSpectrum->Write();
  
  // Refernce covariabce matrix generator object created. This object takes a nominal histogram initially and the output TFile
  RefCovMatGenerator refCov(aInput,*outputFile);
  
  random->SetSeed(0);
  std::clog << nToys << " toys being created" << std::endl;
  
  double sigma = 0.05;
  // Create a histogram to record the energy scale throw values.
  TH1D* throwHist = new TH1D("Throw histogram","",100,-2.5*sigma,2.5*sigma);
  std::vector<TH1D*> tempToyHist;
  
  std::time_t initialTime = std::time(nullptr);
  std::clog << initialTime << " seconds since the Epoch\n";
  for(int k=0; k <nToys;k++)
  {
    // gaussian random pull
    double pull = sigma * random->Gaus(0.0, 1.0);
    //Filling the throw histogram to be used later
    throwHist->Fill(pull);
    
    if(k<20){
      TString name = Form("tempBin%d", k+1);
      TH1D* tempHist = (TH1D*)hObsSpectrum->Clone(name);
      tempHist->Scale(0.0);
      tempToyHist.push_back(tempHist);
    }
    
    if(k%500 == 0)
    {
      std::clog << "Working on toy: " << k <<std::endl;
      std::time_t presentTime = std::time(nullptr);
      std::clog << "dt: " << presentTime - initialTime << " s" <<std::endl;
    }
    
    for(int i=0;i<nEneBins;i++)
    {
      double tempBinContent = hObsSpectrum->GetBinContent(i+1);
      tempBinContent = tempBinContent* (1+pull);
      if(k<20)tempToyHist[k]->SetBinContent(i+1, tempBinContent);
      // Fill input histogram that will be supplied as a toy to the covariance matrix generator
      aInput[i]= tempBinContent;
    }
    
    
    // Supplying a toy to the covariance matrix generator, returns 0 if toy couldn't be used for some reason
    if(refCov.addSample(aInput)==0)
    {
      std::cerr << "Something went wrong, samples cannot be added" << std::endl;
      return 0;
    }
    
    // Printin the toy histograms for first 20 toys
    if (k < 20) {
      TString name = Form("Toy%d", k);
      tempToyHist[k]->SetName(name);
      tempToyHist[k]->Write();
    }
  }
  EneHist->Write();
  throwHist->Write();
  std::clog<< "Saving the file" <<outputName.Data()<<std::endl;
  
  return 0;
}
