// Variables: Number of samples
// Number of bins(energy bins and position bins) of covariance matrix
//
//
#include "TCovarianceMatrixGenerator.hh"
//#include "PlotEditor.hh"
#include <sstream>
// Program defaults to this point if the input is not correct
int Usage()
{
  // The input root file should be the output from generateSetupFiles
  printf("Make sure to input the reference ROOT file and detector number \n");
  printf("Usage: ./generateSignalNormCovMatrix input_ROOT_file detector_number\n");
  printf("Usage: ./generateSignalNormCovMatrix input_ROOT_file detector_number nToys\n");
  return -1;
}

int main(int argc, char** argv)
{
  if(argc!=3 && argc!=4 )
  {
    int errorCode = Usage();
    std::cerr<< "Exiting with error code " << errorCode <<std::endl;
    return errorCode;
  }
  
  TFile* inputFile = new TFile(argv[1]);
  int nToys;
  // If number of toys are input use it, otherwise use 10000 toys
  if(argc==3) {
      nToys = 10000;
    }
  else nToys = atoi(argv[3]);
  
  int detNumber = atoi(argv[2]);
  
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
  TString LvsEName;
		LvsEName.Form("LvsENull%i",detNumber);
  TKey *key = inputFile->FindKey(LvsEName);
  if (key ==0)
  {
    printf("Histogram %s does not exist!!\n",LvsEName.Data());
    return -1;
  }
  
  TH2F* LvsENull = (TH2F*)inputFile->Get(LvsEName);
  
  // Create Output File
  TString outputName(inputFileName);
  outputName.ReplaceAll("Setup", "Cov_Sys_SpecNorm");
  TFile* outputFile = new TFile(outputName,"RECREATE");
  
  // Create vector of Energy Spectra for varying position
  double nEneBins = LvsENull->GetXaxis()->GetNbins();
  double nPosBins = LvsENull->GetYaxis()->GetNbins();
  
  // Total number of bin in covariance matrix are a product of position and energy bins
  double covMatBins = nEneBins*nPosBins;
  
  TH1D* posEneHist = new TH1D("posEneHist", "Position Energy Histograms; Position-Energy bin number; arb", covMatBins, -0.5, covMatBins-0.5);
  
  // Position energy spectrum histogram
  std::vector<TH1D*> LvsEPosSpectrum;
  TString posBins="";
  
  std::clog << "Energy Bins: " << nEneBins << ", Position Bins: " << nPosBins << ", Covariance matrix size: "<< covMatBins <<std::endl;
  
  vector <double> aInput(covMatBins);
  //Initial set of values for reference
  for(int i=0;i<nPosBins;i++)
  {
    posBins=Form("PositionBin%d", i+1);
    LvsEPosSpectrum.push_back((TH1D*)LvsENull->ProjectionX(posBins,i+1,i+1));
    double baseline = LvsENull->GetYaxis()->GetBinCenter(i+1);
    LvsEPosSpectrum[i]->SetTitle(Form("Baseline %f; Energy (MeV)",baseline));
    for(int j=0;j<nEneBins;j++)
    {
      double tempBinContent=0.0;
      tempBinContent = LvsEPosSpectrum[i]->GetBinContent(j+1);
      aInput[j+i*nEneBins]= tempBinContent;
      posEneHist->SetBinContent(j+i*nEneBins,tempBinContent);
    }
    LvsEPosSpectrum[i]->Write();
  }
  
  // Refernce covariabce matrix generator object created. This object takes a nominal histogram initially and the output TFile
  RefCovMatGenerator refCov(aInput,*outputFile);
  
  random->SetSeed(0);
  std::clog << nToys << " toys being created" << std::endl;
  
  double sigma = 1.0;
  // Create a histogram to record the energy scale throw values.
  TH1D* throwHist = new TH1D("Throw histogram","",100,-2.5*sigma,2.5*sigma);
  std::vector<TH1D*> tempToyHists;
  for(int i=0; i <nPosBins;i++)
  {
    TString name = Form("tempBin%d", i+1);
    TH1D* tempHist = (TH1D*)LvsEPosSpectrum[i]->Clone(name);
    tempHist->Scale(0.0);
    tempToyHists.push_back(tempHist);
  }
  
  std::time_t initialTime = std::time(nullptr);
  std::clog << initialTime << " seconds since the Epoch\n";
  for(int k=0; k <nToys;k++)
  {
    // gaussian random pull
    double pull = sigma * random->Gaus(0.0, 1.0);
    //Filling the throw histogram to be used later
    throwHist->Fill(pull);
    
    if(k%500 == 0)
    {
      std::clog << "Working on toy: " << k <<std::endl;
      std::time_t presentTime = std::time(nullptr);
      std::clog << "dt: " << presentTime - initialTime << " s" <<std::endl;
    }
    
    for(int i=0;i<nPosBins;i++)
    {
      for(int j=0;j<nEneBins;j++)
      {
        
        double tempBinContent = LvsEPosSpectrum[i]->GetBinContent(j+1);
        tempBinContent = tempBinContent* (1+pull);
        tempToyHists[i]->SetBinContent(j+1, tempBinContent);
        // Fill input histogram that will be supplied as a toy to the covariance matrix generator
        aInput[j+i*nEneBins]= tempBinContent;
      }
    }
    
    // Supplying a toy to the covariance matrix generator, returns 0 if toy couldn't be used for some reason
    if(refCov.addSample(aInput)==0)
    {
      std::cerr << "Something went wrong, samples cannot be added" << std::endl;
      return 0;
    }
    
    // Printin the toy histograms for first 20 toys
    if (k < 20) {
      TString name = Form("Toy%dBin%d", k, 4);
      tempToyHists[3]->SetName(name);
      tempToyHists[3]->Write();
    }
  }
  posEneHist->Write();
  LvsENull->Write();
  throwHist->Write();
  std::clog<< "Saving the file" <<outputName.Data()<<std::endl;
  
  return 0;
}
