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
#include "TKey.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "PROSPECTStyle.hh"

void usage(int nInputs)
{
  printf("Are you serious!, you think this program could run with %i inputs\n",nInputs);
  printf("Usage: ./DrawChiPicture inputfile\n");
  exit(0);
}

int main(int argc, char** argv)
{
  if(argc !=2)
  {
    usage(argc);
    return -1;
  }
  
  TString inputFileName(argv[1]);
  TFile* inputFile = TFile::Open(inputFileName);
  
  //Load style
  SetupProspectStyle();
  
  // Make Canvas
  TCanvas *c1 = new TCanvas("c1","c1",0,0,700,600);

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
  
  // Set contour level
  double contourSig[1] = {10.83};
  
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
  
  TList *listOfKeys=inputFile->GetListOfKeys();
  for (int i = 0; i < listOfKeys->GetSize(); i++) {
    TString keyName=((TKey*)listOfKeys->At(i))->GetName();
    if(!(keyName.Contains("CumulativeChiMap"))) continue;
    keyName.ReplaceAll("CumulativeChiMap","");
    TH2F* temp;
    TString HistogramName;
    std::cout << keyName <<std::endl;
    HistogramName.Form("CumulativeChiMap%i",keyName.Atoi());
    temp = (TH2F*)inputFile->Get(HistogramName.Data());
    temp->SetLineWidth(2);
    temp->SetContour(1, contourSig);
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
  for (unsigned int i = 0; i < chiSqGraphs.size(); i++) {
    chiSqGraphs[i]->Draw("cont LIST");
    c1->Update();
    
    TObjArray *contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList *list = (TList*)contours->At(0);
    
    coverage3SigTemp = (TGraph*)list->First();
    coverage3Sig[i] =(TGraph*) coverage3SigTemp->Clone();
  }
  
  for(int i=0;i<4;i++){
    K95[i]->Draw("FSame");
    K90[i]->Draw("FSame");
  }
  
  for(int i=0;i<2;i++){
    RAA95[i]->Draw("FSame");
    RAA95[i]->Draw("CSame");
  }
  
  for (unsigned int i = 0; i < chiSqGraphs.size(); i++) {
    chiSqGraphs[i]->Draw("cont3SAME");
  }
  
  c1->SetLogx();
  c1->SetLogy();
  
  double xpoint[1] = {0.09};
  double ypoint[1] = {1.78};
  TGraph *KPoint = new TGraph(1,xpoint,ypoint);
  KPoint->SetMarkerStyle(29);
  KPoint->SetMarkerSize(1.5);
  KPoint->SetMarkerColor(kBlack);
  KPoint->Draw("PSame");
  
  TGraph *KPoint2 = new TGraph(1,xpoint,ypoint);
  KPoint2->SetMarkerStyle(24);
  KPoint2->SetMarkerSize(1.5);
  KPoint2->SetMarkerColor(kBlack);
  KPoint2->Draw("PSame");
  
  double RAAxpoint[1] = {0.165};
  double RAAypoint[1] = {2.39};
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
  
  
  TString outputTitle(inputFileName);
  outputTitle.ReplaceAll(".root","");
  c1->SaveAs(outputTitle + ".png");
  c1->SaveAs(outputTitle + ".eps");
  c1->SaveAs(outputTitle + ".pdf");
  c1->SaveAs(outputTitle + ".C");
  
}
