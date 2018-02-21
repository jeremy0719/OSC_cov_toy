#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <map>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TList.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TOscillationModelBuilder.hh"
#include "PROSPECTStyle.hh"


void EstimateBinNumbers(double deltam2Value, double sinSq2ThetaValue, int& deltam2Bin,int& sin22theta2Bin)
{
  int fNDeltam2=57;
  int fNSinSq2Theta=49;
  deltam2Bin = (std::log10(deltam2Value)+1.425)*fNDeltam2/(std::log10(28.2)+1.425);
  sin22theta2Bin = (std::log10(sinSq2ThetaValue)+2.425)*fNSinSq2Theta/(std::log10(2.425)+2.4);
}

void usage(int nInputs)
{
  printf("Are you serious!, you think this program could run with %i inputs\n",nInputs);
  printf("Usage: ./DrawChiPicture list_of_files list_of_legends_files CriticalChi2File.root\n");
}

std::vector<TFile*> UnpackFiles(TString inputList) {
  std::vector<TFile*> outputVector;
  outputVector.clear();
  
  std::ifstream inputFile(inputList.Data());
  if (!inputFile.is_open()) {
    printf("File %s not found!",inputList.Data());
    exit(1);
  }
  
  TString fileName;
  while (!inputFile.eof()) {
    inputFile >> fileName;
    // Need in case of empty line, ensures file name exists.
    if (fileName.IsWhitespace()) continue;
    TFile* temp = new TFile(fileName.Data(), "READ");
    if (!temp->IsOpen() || temp->IsZombie()) exit(1);
    outputVector.push_back(temp);
  }
  
  inputFile.close();
  return outputVector;
}

std::vector<TString> UnpackLegend(TString inputList) {
  std::vector<TString> outputVector;
  outputVector.clear();
  
  std::ifstream inputFile(inputList.Data());
  if (!inputFile.is_open()) {
    printf("File %s not found!",inputList.Data());
    exit(1);
  }
  
  std::string line;
  while (!inputFile.eof()) {
    std::getline(inputFile, line);
    TString legendTitle(line);
    // Need in case of empty line, ensures file name exists.
    if (legendTitle.IsWhitespace()) continue;
    outputVector.push_back(legendTitle);
  }
  
  inputFile.close();
  return outputVector;
}

int main(int argc, char** argv)
{
  SetupProspectStyle();
  if(argc !=4)
  {
    usage(argc);
    return 1;
  }
  TString inputList= argv[1];
  TString legendList= argv[2];
  TString CriticalChi2FileName = argv[3];
  //Load style
  SetupProspectStyle();
  
  // Make Canvas
  TCanvas *c1 = new TCanvas("c1","c1",0,0,700,600);
  c1->SetTickx();
  c1->SetTicky();
  
  // Check that the input file name ends with .list
  // It should be a list of .root files that you want to analyze together
  if (!inputList.EndsWith(".list") || !legendList.EndsWith(".list")) {
    printf("Incorrect argument structure.\n");
    printf(" Argument: %s must end in '.list' \n",inputList.Data());
    printf(" Argument: %s must end in '.list' \n",legendList.Data());
    exit(1);
  }
  
  //Get previously made energy-position graphs
  std::vector<TFile*> inputFiles = UnpackFiles(inputList);
  
  //Get Legend Names
  std::vector<TString> legendTitles = UnpackLegend(legendList);
  
  // Check that the inputs have the same size.
  if (inputFiles.size() != legendTitles.size()) {
    printf("File sizes not the same!\n");
    printf("root files list: %lu\n",inputFiles.size());
    printf("legends list   : %lu\n",legendTitles.size());
    exit(1);
  }
  
  // Get comparison regions.
  std::ifstream infile;
  double Kopp90x[4][2000];  double Kopp90y[4][2000];
  double Kopp95x[4][2000];  double Kopp95y[4][2000];
  std::map<double,double> Kopp95Map[4];
  
  double RAA95x[2][2000];  double RAA95y[2][2000];
  
  ///////Kopp Regions
  int Kpoints90[4];
  int Kpoints95[4];
  for(int j=0;j<4;j++){
    infile.open(Form("inputs/Kopp2013_SBL_%d.txt",j+1));
    int i=0;
    while(!infile.eof()){
      infile >> Kopp95x[j][i] >> Kopp95y[j][i];
      Kopp95y[j][i] = 4*Kopp95y[j][i]*(1-Kopp95y[j][i]);
      (Kopp95Map[j])[Kopp95y[j][i]]=Kopp95x[j][i];
      i++;
    }
    Kpoints95[j] = i;
    infile.close();
  }
  
  printf("Kopp SBL region loaded.\n");
  
  for(int j=0;j<4;j++){
    infile.open(Form("inputs/Kopp2013_nue_%d.txt",j+1));
    int i=0;
    while(!infile.eof()){
      infile >> Kopp90x[j][i] >> Kopp90y[j][i];
      Kopp90y[j][i] = 4*Kopp90y[j][i]*(1-Kopp90y[j][i]);
      i++;
    }
    Kpoints90[j] = i;
    infile.close();
  }
  
  printf("Kopp nue region loaded.\n");
  
  int RAApoints95[2];
  for(int j=0;j<2;j++){
    infile.open(Form("inputs/WhitePaper95_%d.txt",j+1));
    int i=0;
    while(!infile.eof()){
      infile >> RAA95x[j][i] >> RAA95y[j][i];
      i++;
    }
    RAApoints95[j] = i;
    infile.close();
  }
  
  printf("LSN Whitepaper region loaded.\n");
  
  c1->cd();
  
  // Construct the frame
  TH1F *hr = c1->DrawFrame(0.0045,0.04,1,22);
  hr->GetXaxis()->SetTitle("sin^{2}2#theta_{ee}");
  hr->GetYaxis()->SetTitle("#Deltam^{2}_{14} [eV^{2}]");
  
  std::vector<TH2F*> chiSqGraphs;
  chiSqGraphs.clear();
  
  // Define pretty colors
  std::vector<Color_t> colors;
  colors.clear();
  colors.push_back(kBlack);
  colors.push_back(kMagenta+2);
  colors.push_back(kCyan+2);
  colors.push_back(kPink-2);
  colors.push_back(kGreen-2);
  colors.push_back(kOrange+7);
  colors.push_back(kGray+1);
  colors.push_back(kOrange+2);
  colors.push_back(kAzure-2);
  colors.push_back(kViolet+2);
  
  // Set contour level
  double contourSig[1] = {0};
  for (unsigned int i = 0; i < inputFiles.size(); i++) {
    TH2F* temp;
    temp = (TH2F*)inputFiles[i]->Get("CumulativeChiMap5");
    if (i < colors.size()) temp->SetLineColor(colors[i]);
    else temp->SetLineColorAlpha(kRed-i,0.8);
    temp->SetLineWidth(2);
    chiSqGraphs.push_back(temp);
  }
  
  TGraph *K90[4];
  TGraph *K95[4];
  for(int i=0;i<4;i++){
    K90[i] = new TGraph(Kpoints90[i]-1, Kopp90y[i], Kopp90x[i]);
    K95[i] = new TGraph(Kpoints95[i]-1, Kopp95y[i], Kopp95x[i]);
    K90[i]->SetFillColorAlpha(kRed-4, 0.25);
    K95[i]->SetFillColorAlpha(kBlue-7, 0.25);
  }
  
  TGraph *RAA95[2];
  for(int i=0;i<2;i++){
    RAA95[i] = new TGraph(RAApoints95[i]-1, RAA95x[i], RAA95y[i]);
    RAA95[i]->SetFillColorAlpha(kGreen+2, 0.20);
    RAA95[i]->SetLineColor(kGreen+2);
    RAA95[i]->SetLineWidth(1);
  }
  
  c1->Clear();
  TGraph *coverage3SigTemp;
  std::map<int, TGraph*> coverage3Sig;

  
  for(int i=0;i<4;i++){
    K95[i]->Draw("FSame");
    K90[i]->Draw("FSame");
  }
  
  for(int i=0;i<2;i++){
    RAA95[i]->Draw("FSame");
    RAA95[i]->Draw("CSame");
  }
  
  
  double Koppxpoint[1] = {0.09};
  double Koppypoint[1] = {1.78};
  double RAAxpoint[1] = {0.165};
  double RAAypoint[1] = {2.39};
  
  TFile *CriticalChi2File=TFile::Open(CriticalChi2FileName.Data());
  TH2D *CriticalChi2Map=(TH2D*)CriticalChi2File->Get("hCriticalChi2Map");
  
  TGraph* chi2Graphs[chiSqGraphs.size()];
  
  for (unsigned int i = 0; i < chiSqGraphs.size(); i++) {
    
    std::cout << "test " <<std::endl;
    //chiSqGraphs[i]->Draw("COLZ");
    TString histName("temp");
    histName.Append(i);
    //TH2D* tempHist=(TH2D*)chiSqGraphs[i]->Clone(histName);
    TH2D* tempHist=(TH2D*)CriticalChi2Map->Clone(histName);
//    tempHist->Add(CriticalChi2Map,-1);
      tempHist->Add(chiSqGraphs[i],-1);
    for(int j=1;j<=tempHist->GetXaxis()->GetNbins();j++){
      for(int k=1;k<=tempHist->GetYaxis()->GetNbins();k++){
        //if(tempHist->GetBinContent(j,k)>(0)) tempHist->SetBinContent(j,k,0);
      }
    }
    
    tempHist->Draw("cont LIST");
    tempHist->SetContour(1, contourSig);
    tempHist->Draw("CONT3 same");
    std::cout << legendTitles[i] << std::endl;
    
    // Calculate Kopp best-fit sigma values
    //int KoppBinX(0);
    //int KoppBinY(0);
    //EstimateBinNumbers(Koppypoint[0],Koppxpoint[0],KoppBinY,KoppBinX);
    double chiSquare = chiSqGraphs[i]->GetBinContent(28, 35);
    //std::cout << chiSqGraphs[i]->GetXaxis()->GetBinCenter(28)<<std::endl;
    //std::cout << chiSqGraphs[i]->GetYaxis()->GetBinCenter(35)<<std::endl;
    double sigma = TMath::ErfcInverse(TMath::Prob(chiSquare, 2)) * TMath::Sqrt(2);
    printf(" Kopp Best Fit Sigma:    %f\n " ,sigma);
    
    // Calculate RAA best-fit sigma values
    chiSquare = chiSqGraphs[i]->GetBinContent(33, 38);
    //std::cout << chiSqGraphs[i]->GetXaxis()->GetBinCenter(33)<<std::endl;
    //std::cout << chiSqGraphs[i]->GetYaxis()->GetBinCenter(38)<<std::endl;
    sigma = TMath::ErfcInverse(TMath::Prob(chiSquare, 2)) * TMath::Sqrt(2);
    printf(" RAA Best Fit Sigma:    %f\n " ,sigma);
    double minBin = chiSqGraphs[i]->FindFirstBinAbove(11.83, 1);
    double minAngle = chiSqGraphs[i]->GetXaxis()->GetBinCenter(minBin);
    printf(" Lowest sin22theta :    %f\n " ,minAngle);
  }
  
  c1->SetLogx();
  c1->SetLogy();
  
  TGraph *KPoint = new TGraph(1,Koppxpoint,Koppypoint);
  KPoint->SetMarkerStyle(29);
  KPoint->SetMarkerSize(1.5);
  KPoint->SetMarkerColor(kBlack);
  KPoint->Draw("PSame");
  
  TGraph *KPoint2 = new TGraph(1,Koppxpoint,Koppypoint);
  KPoint2->SetMarkerStyle(24);
  KPoint2->SetMarkerSize(1.5);
  KPoint2->SetMarkerColor(kBlack);
  KPoint2->Draw("PSame");
  
  TGraph *RAAPoint = new TGraph(1,RAAxpoint,RAAypoint);
  RAAPoint->SetMarkerStyle(29);
  RAAPoint->SetMarkerSize(1.5);
  RAAPoint->SetMarkerColor(kGreen+2);
  RAAPoint->Draw("PSame");
  
  TGraph *RAAPoint2 = new TGraph(1,RAAxpoint,RAAypoint);
  RAAPoint2->SetMarkerStyle(24);
  RAAPoint2->SetMarkerSize(1.5);
  RAAPoint2->SetMarkerColor(kBlack);
  RAAPoint2->Draw("PSame");
  
  TLegend* leg = new TLegend(0.12,0.12,0.42,0.34);
  leg->SetTextFont(132);
  leg->SetTextSize(0.03);
  leg->SetHeader("#font[22]{Sensitivity: 3#kern[1]{#sigma} CL}");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  for (unsigned int i = 0; i < legendTitles.size(); i++) {
    leg->AddEntry(chiSqGraphs[i], legendTitles[i], "l");
  }
  
  K90[0]->SetMarkerStyle(29);
  K90[0]->SetMarkerSize(1.5);
  K90[0]->SetMarkerColor(kBlack);
  
  RAA95[0]->SetMarkerStyle(29);
  RAA95[0]->SetMarkerSize(1.5);
  RAA95[0]->SetMarkerColor(kGreen+2);
  
  leg->AddEntry(K95[0],"SBL Anomaly (Kopp), 95% CL","f");
  leg->AddEntry(K90[0],"All #nu_{e} Disappearance Exps (Kopp), 95% CL","fp");
  leg->AddEntry(RAA95[0],"SBL + Gallium Anomaly (RAA), 95% CL","fp");
  leg->Draw();
  
  
  TString outputTitle(inputList);
  outputTitle.ReplaceAll(".list","");
  c1->SaveAs(outputTitle + ".png");
  c1->SaveAs(outputTitle + ".eps");
  c1->SaveAs(outputTitle + ".pdf");
  c1->SaveAs(outputTitle + ".C");
  
}
