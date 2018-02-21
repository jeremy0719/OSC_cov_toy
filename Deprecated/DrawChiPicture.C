#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <math.h>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"


#include "prospect_style.H"

std::vector<TFile*> UnpackFiles(TString inputList) {
  std::vector<TFile*> outputVector;
  outputVector.clear();
  
  ifstream inputFile(inputList.Data());
  if (!inputFile.is_open()) {
    std::cout << "File " << inputList << " not found!" << std::endl;
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
  
  ifstream inputFile(inputList.Data());
  if (!inputFile.is_open()) {
    std::cout << "File " << inputList << " not found!" << std::endl;
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

void DrawChiPicture(TString inputList, TString legendList){
  
  //Load style
  SetupProspectStyle();
  
  // Make Canvas
  TCanvas *c1 = new TCanvas("c1","c1",0,0,700,600);
  c1->SetLogx();
  c1->SetLogy();
  
  
  // Check that the input file name ends with .list
  // It should be a list of .root files that you want to analyze together
  if (!inputList.EndsWith(".list") || !legendList.EndsWith(".list")) {
    std::cout << "Incorrect argument structure." << std::endl;
    std::cout << " Argument: " << inputList << " must end in '.list' " << std::endl;
    std::cout << " Argument: " << legendList << " must end in '.list' " << std::endl;
    exit(1);
  }
  
  //Get previously made energy-position graphs
  std::vector<TFile*> inputFiles = UnpackFiles(inputList);
  
  //Get Legend Names
  std::vector<TString> legendTitles = UnpackLegend(legendList);
  
  // Check that the inputs have the same size.
  if (inputFiles.size() != legendTitles.size()) {
    std::cout << "File sizes not the same!" << std::endl;
    std::cout << "root files list: " << inputFiles.size() << std::endl;
    std::cout << "legends list   : " << legendTitles.size() << std::endl;
    exit(1);
  }

  // Get comparison regions.
  ifstream infile;
  double Kopp90x[4][2000];  double Kopp90y[4][2000];  
  double Kopp95x[4][2000];  double Kopp95y[4][2000];  

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
      i++;
    }
    Kpoints95[j] = i;
    infile.close();
  }    

  std::cout << "Kopp SBL region loaded." << std::endl;
  
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
  
  std::cout << "Kopp nue region loaded." << std::endl;
  
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
  
  std::cout << "LSN Whitepaper region loaded." << std::endl;
  
  c1->cd();

  // Set contour level
  double contour3Sig[1] = {11.83}; // 3 sigma
  
  // Construct the frame
  TH1F *hr = c1->DrawFrame(0.0045,0.04,1,22);
  hr->GetXaxis()->SetTitle("sin^{2}2#theta_{ee}");
  hr->GetYaxis()->SetTitle("#Deltam^{2}_{14} [eV^{2}]");
  
  std::vector<TH2F*> chiSqGraphs;
  chiSqGraphs.clear();
  
  if (inputFiles.size() > 9) {
    std::cout << "Too many inputs, not enough line styles!" << std::endl;
    exit(1);
  }
  
  int LineStyles[9] = {1, 3, 7, 5, 2, 8, 9, 6, 10};
  
  std::vector<Color_t> colors;
  colors.clear();
  colors.push_back(kBlack);
  colors.push_back(kMagenta+2);
  colors.push_back(kAzure-2);
  colors.push_back(kGreen+2);
  colors.push_back(kOrange+7);
  colors.push_back(kCyan-4);
  colors.push_back(kViolet+1);
  colors.push_back(kBlack);
  colors.push_back(kBlack);
  colors.push_back(kBlack);
  
  for (unsigned int i = 0; i < inputFiles.size(); i++) {
    TH2F* temp = (TH2F*)inputFiles[i]->Get("ChiSquareMap");
    
    temp->SetLineColor(colors[i]);
    temp->SetLineWidth(2);
    temp->SetContour(1, contour3Sig);
    temp->SetLineStyle(LineStyles[i]);

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

  for(int i=0;i<4;i++){
    K95[i]->Draw("FSame");
    K90[i]->Draw("FSame");
  }

  TGraph *RAA95[2];
  for(int i=0;i<2;i++){
    RAA95[i] = new TGraph(RAApoints95[i]-1, RAA95x[i], RAA95y[i]);
    RAA95[i]->SetFillColorAlpha(kGreen+2, 0.20);
    RAA95[i]->SetLineColor(kGreen+2);
    RAA95[i]->SetLineWidth(1);
  }
  
  for(int i=0;i<2;i++){
    RAA95[i]->Draw("FSame");
    RAA95[i]->Draw("CSame");
  }
  
  for (unsigned int i = 0; i < chiSqGraphs.size(); i++) {
    chiSqGraphs[i]->Draw("cont3Same");
    std::cout << legendTitles[i] << std::endl;
    // Calculate the Kopp Best Fit Value
    double chiSquare = chiSqGraphs[i]->GetBinContent(29, 34);
    double sigma = TMath::ErfcInverse(TMath::Prob(chiSquare, 2)) * TMath::Sqrt(2);
    std::cout << " Kopp Best Fit Sigma:     " << sigma << std::endl;
    // Calculate the RAA Best Fit Value
    chiSquare = chiSqGraphs[i]->GetBinContent(35, 36);
    sigma = TMath::ErfcInverse(TMath::Prob(chiSquare, 2)) * TMath::Sqrt(2);
    std::cout << " RAA Best Fit Sigma:     " << sigma << std::endl;
    // Calculate the lowest sin2 2theta
    double minBin = chiSqGraphs[i]->FindFirstBinAbove(11.83, 1);
    double minAngle = chiSqGraphs[i]->GetXaxis()->GetBinCenter(minBin);
    std::cout << " Lowest sin2 2theta: " << minAngle << std::endl;
  }
  
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
 
  for (unsigned int i = 0; i < legendTitles.size(); i++) {
    leg->AddEntry(chiSqGraphs[i], legendTitles[i], "l");
  }
  
  K90[0]->SetMarkerStyle(29);
  K90[0]->SetMarkerSize(1.5);
  K90[0]->SetMarkerColor(kBlack);
  
  RAA95[0]->SetMarkerStyle(29);
  RAA95[0]->SetMarkerSize(1.5);
  RAA95[0]->SetMarkerColor(kGreen+2);
  
  leg->AddEntry(K95[0],"SBL Anomaly (Kopp), 95\% CL","f");
  leg->AddEntry(K90[0],"All #nu_{e} Disappearance Exps (Kopp), 95\% CL","fp");
  leg->AddEntry(RAA95[0],"SBL + Gallium Anomaly (RAA), 95\% CL","fp");
  leg->Draw();
  
  
  TString outputTitle(inputList);
  outputTitle.ReplaceAll(".list","");
  c1->SaveAs(outputTitle + ".png");
  c1->SaveAs(outputTitle + ".eps");
  c1->SaveAs(outputTitle + ".pdf");
  c1->SaveAs(outputTitle + ".C");
  
}
