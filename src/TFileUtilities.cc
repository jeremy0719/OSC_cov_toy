#include "TFileUtilities.hh"

void usage()
{
  std::cout <<"Input has to include the name of the ROOT file including the location and the name of the histogram " <<std::endl;
}

// Extract Data from a 1D histogram with a given name, takes name of the ROOT file and name of the histogram in the ROOT file in that order
void hDataExtractor(TString fileName, TString histogramName)
{
  TFile* file=TFile::Open(fileName.Data());
  if(!file->IsOpen() || !file->IsZombie())
  {
    usage();
    return;
  }
  TH1* histogram=(TH1*)file->Get(histogramName.Data());
  for(int i=1;i<histogram->GetNbinsX();i++)
  {
    std::cout <<i*histogram->GetBinWidth(i) << ","<<histogram->GetBinContent(i)<<std::endl;
  }
}

// Extract Data from a 2D histogram with a given name, takes name of the ROOT file and name of the histogram in the ROOT file in that order
void tdhDataExtractor(TString fileName,TString histogramName)
{
  TFile* file=TFile::Open(fileName.Data());
  if(!file->IsOpen() || !file->IsZombie())
  {
    usage();
    return;
  }
  TH2* histogram=(TH2*)file->Get(histogramName.Data());
  for(int i=1;i<histogram->GetNbinsX();i++)
  {
    for(int j=1;j<histogram->GetNbinsX();j++)
    {
      std::cout <<i*(histogram->GetXaxis())->GetBinWidth(i)<< "," << j*histogram->GetYaxis()->GetBinWidth(j)<< ","<<histogram->GetBinContent(i,j)<<std::endl;
    }
  }
}

TH2D* makeTXTTH2(TString fileName,int XBins, int YBins)
{
  TString histName = fileName;
  histName.ReplaceAll("./","");
  histName.ReplaceAll(".txt","");
  TH2D* tdhOutput = new TH2D(histName.Data(),"Best shift values",XBins,0,XBins,YBins,0,YBins);
  
  std::ifstream in;
  //    int nlines=0;
  in.open(fileName);
  printf("Reading input file '%s'\n",fileName.Data());
  
  int segNo;
  double z;
  
  while (true)
  {
    in >> segNo >> z;
    if (!in.good()) break;
    std::cout << (int)(segNo%XBins) << " " <<(int)(segNo/XBins) <<" " << std::abs(z) <<std::endl;
    tdhOutput->SetBinContent(segNo%XBins+1,(int)(segNo/XBins)+1,std::abs(z));
  }
  return tdhOutput;
}


void makeTXTTGraph(TString fileName, TGraph *tGOutput)
{
  TString histName = fileName;
  histName.ReplaceAll("./","");
  histName.ReplaceAll("/","");
  histName.ReplaceAll(".txt","");
  tGOutput->SetName(histName.Data());
  
  std::ifstream infile;
  std::string lineRead;
  
  infile.open(fileName.Data());
  printf("Reading input file '%s'\n",fileName.Data());
  
  double numberRead=0.0;
  double x=0.0;
  double y=0.0;
  int nlines=0;
  while(infile.good()){
    while(getline(infile, lineRead)){
      std::istringstream streamA(lineRead);
      int columncount =0;
      while(streamA >>numberRead){
        if(columncount >=2){
          printf("More than 2 points in a line\n");
          break;
        }
        (columncount==0)?x=numberRead:y=numberRead;
        columncount++;
      }
      tGOutput->SetPoint(nlines,x,y);
      ++nlines;
    }
  }
  infile.close();
}
