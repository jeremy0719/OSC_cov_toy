#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <map>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TList.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "PROSPECTStyle.hh"

void usage(int nInputs)
{
  printf("Are you serious!, you think this program could run with %i inputs\n",nInputs);
  printf("Usage: ./PlotFluxDeficit \n");
}

int main(int argc, char** argv)
{
  if(argc !=1)
  {
    usage(argc);
    return -1;
  }
  
  // Remove stupid warnings.
  TString dummy = argv[0];

  //Load style
  SetupProspectStyle();
  gStyle->SetPadTopMargin(0.4);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadBottomMargin(0.2);
  
  gStyle->SetTitleSize(0.07, "xyz");
  gStyle->SetLabelSize(0.07,"xyz");
  
  // Make Canvas
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1200,700);
  
  // Get flux deficit arrays.
  std::ifstream infile;
  std::vector<TString> vExpNames;
  std::vector<double> vRatios;
  std::vector<double> vRatioErrs;
  std::vector<double> vDistances;
  
  infile.open("inputs/FluxDeficit.txt", std::ios::in);

  // Loop through file and store name and variables to vectors
  std::string line;
  while (getline(infile,line)) {
    TString sExpName;
    double dRat, dRatUnc, dRatCor, dDist;
    std::stringstream splitline;
    splitline<<line;
    splitline >> sExpName >> dRat >> dRatUnc >> dRatCor >> dDist;

    if (sExpName == "Experiment") continue;

    vExpNames.push_back(sExpName);
    vRatios.push_back(dRat);
    vRatioErrs.push_back(TMath::Sqrt(TMath::Power(dRatUnc,2) + TMath::Power(dRatCor, 2)));
    vDistances.push_back(dDist);
  }
  infile.close();
  
  std::cout << "Loaded " << vExpNames.size() << " Ratio Values" << std::endl;
  
  // Load the Flux Deficit graphs.
  TFile* fInputSterileFlux = new TFile("./aux/SterileFluxDeficit.root", "READ");
  if (!fInputSterileFlux->IsOpen() || fInputSterileFlux->IsZombie()) {
    std::cout << "File SterileFluxDeficit.root not found!" << std::endl;
    std::cout << "Returning..." << std::endl;
  }
  
  // Get the L/E graph of three neutrino oscillation
  TGraph* gThreeNu = (TGraph*)fInputSterileFlux->Get("gRatio3");
  gThreeNu->SetLineColor(kBlue-3);
  gThreeNu->SetLineWidth(2);
  
  // The the L/E graph of four neutrino oscillation
  TGraph* gFourNu = (TGraph*)fInputSterileFlux->Get("gRatio4");
  gFourNu->SetLineColor(kViolet-3);
  gFourNu->SetLineWidth(2);
  
  std::cout << "Loaded oscillation graphs." << std::endl;
  
  // Construct the frame
  TH1F *hr = c1->DrawFrame(1,0.6,2000,1.3);
  hr->GetXaxis()->SetTitle("L (m)");
  hr->GetXaxis()->SetTitleOffset(1.2);
  hr->GetYaxis()->SetTitle("N_{meas} / N_{HM}");
  
  // Draw the 1.0 line
  TLine* lNullOsc = new TLine(1.0, 1.0, 2000, 1.0);
  lNullOsc->SetLineWidth(2);
  lNullOsc->SetLineStyle(2);
  lNullOsc->Draw();
  
  // Draw the oscillated lines
  gThreeNu->Draw("CSAME");
  gFourNu->Draw("CSAME");
  
  c1->SetLogx();

  // Set-up map between names and graphs
  std::map<TString, TGraphErrors*> mExpGraphs;
  
  // Define the unique names of the experiments.
  // This groups together arrays of experiments.
  std::vector<TString> vUniqNames;
  vUniqNames.push_back("Bugey-3");
  vUniqNames.push_back("Bugey-4");
  vUniqNames.push_back("Rovno88");
  vUniqNames.push_back("Rovno91");
  vUniqNames.push_back("Gosgen");
  vUniqNames.push_back("ILL");
  vUniqNames.push_back("Krasnoyarsk");
  vUniqNames.push_back("SRP");
  vUniqNames.push_back("Nucifer");
  vUniqNames.push_back("LBL");
  
  // set Draw options for experiments
  int iMarkerStyle[10] = {29, 30, 26, 22, 32, 23, 33, 27, 25, 20};
  double dMarkerSize[10] = {3, 3, 2.5, 2.5, 2.5, 2.5, 3, 3, 2, 2};
  double dMarkerColor[10] = {kGreen+2, kGreen+2, kMagenta+2, kMagenta+2, kViolet+2,
                                            kViolet+2, kBlue+2, kBlue+2, kBlue+2, kBlack};
  
  // Build the different graphs into a map.
  for (unsigned int i = 0; i < vUniqNames.size(); i++) {
    TGraphErrors *geTemp = new TGraphErrors();
    
    geTemp->SetName(vUniqNames[i]);
    geTemp->SetMarkerStyle(iMarkerStyle[i]);
    geTemp->SetMarkerSize(dMarkerSize[i]);
    geTemp->SetMarkerColor(dMarkerColor[i]);
    
    mExpGraphs[vUniqNames[i]] = geTemp;
    
  }
  
  
  // Loop through the experiments and add them to the plot
  for (unsigned int i = 0; i < vExpNames.size(); i++) {
    double Baseline = vDistances[i];
    double BaselineErr = 0.0;
    double Ratio = vRatios[i];
    double RatioErr = vRatioErrs[i]/100.;
    
    // Loop through unique names and add point to correct graph
    for (unsigned int j = 0; j < vUniqNames.size(); j++) {
      if (vExpNames[i].Contains(vUniqNames[j])) {
      
        int pointIndex = mExpGraphs[vUniqNames[j]]->GetN();
      
        mExpGraphs[vUniqNames[j]]->SetPoint(pointIndex, Baseline, Ratio);
        mExpGraphs[vUniqNames[j]]->SetPointError(pointIndex, BaselineErr, RatioErr);
      }
      // Special case for LongBaselines
      else if (vExpNames[i].Contains("DayaBay") || vExpNames[i].Contains("PaloVerde") ||
               vExpNames[i].Contains("RENO") || vExpNames[i].Contains("Chooz") ) {
        int pointIndex = mExpGraphs["LBL"]->GetN();
        
        mExpGraphs["LBL"]->SetPoint(pointIndex, Baseline, Ratio);
        mExpGraphs["LBL"]->SetPointError(pointIndex, BaselineErr, RatioErr);
      }
    } // End loop through unique names
  } // Loop through experiments.
  
  // Draw Points
  for (unsigned int j = 0; j < vUniqNames.size(); j++) {
    mExpGraphs[vUniqNames[j]]->Draw("P");
  }
  
  // Create legend Must split in half for style
  
  TLegend* leg1 = new TLegend(0.15,0.6,0.5,0.95);
  leg1->SetTextFont(132);
  leg1->SetTextSize(0.05);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  
  
  // Add to legend
  for (unsigned int j = 0; j < 4; j++) {
    leg1->AddEntry(mExpGraphs[vUniqNames[j]], vUniqNames[j],"p");
  }
  
  TLegend* leg2 = new TLegend(0.4,0.6,0.95,0.95);
  leg2->SetTextFont(132);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  
  
  // Add to legend
  for (unsigned int j = 4; j < 8; j++) {
    leg2->AddEntry(mExpGraphs[vUniqNames[j]], vUniqNames[j],"p");
  }
  
  TLegend* leg3 = new TLegend(0.67,0.6,0.95,0.95);
  leg3->SetTextFont(132);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  
  
  // Add to legend
  for (unsigned int j = 8; j < vUniqNames.size(); j++) {
    leg3->AddEntry(mExpGraphs[vUniqNames[j]], vUniqNames[j],"p");
  }
  leg3->AddEntry(gThreeNu, "3-#nu Osc.", "l");
  leg3->AddEntry(gFourNu, "4-#nu Osc.", "l");
  
  leg1->Draw();
  leg2->Draw();
  leg3->Draw();
  
  TString outputTitle("GlobalFluxDeficitSummary");
  outputTitle.ReplaceAll(".list","");
  c1->SaveAs(outputTitle + ".png");
  c1->SaveAs(outputTitle + ".eps");
  c1->SaveAs(outputTitle + ".pdf");
  c1->SaveAs(outputTitle + ".C");
  
  
  // Make a version without individual Experiment names.
  // Construct the frame
  gStyle->SetPadTopMargin(0.05);

  // Make Canvas
  TCanvas *c2 = new TCanvas("c2","c2",0,0,1200,500);

  TH1F *hr2 = c2->DrawFrame(1,0.6,2000,1.3);
  hr2->GetXaxis()->SetTitle("L (m)");
  hr2->GetXaxis()->SetTitleOffset(1.2);
  hr2->GetYaxis()->SetTitle("N_{meas} / N_{HM}");
  
  // Draw the 1.0 line
  lNullOsc->Draw();
  
  // Draw the oscillated lines
  gThreeNu->Draw("CSAME");
  gFourNu->Draw("CSAME");
  
  c2->SetLogx();
  
  // Loop through experiments and make one TGraph
  TGraphErrors* gAllExps = new TGraphErrors();
  gAllExps->SetName("AllExps");
  gAllExps->SetMarkerStyle(20);
  gAllExps->SetMarkerSize(0.75);
  gAllExps->SetMarkerColor(kBlack);
  
  // Loop through the experiments and add them to the plot
  for (unsigned int i = 0; i < vExpNames.size(); i++) {
    double Baseline = vDistances[i];
    double BaselineErr = 0.0;
    double Ratio = vRatios[i];
    double RatioErr = vRatioErrs[i]/100.;
        
    gAllExps->SetPoint(i, Baseline, Ratio);
    gAllExps->SetPointError(i, BaselineErr, RatioErr);

  } // Loop through experiments.
  
  gAllExps->Draw("PSAME");
  
  
  TLegend* legAll = new TLegend(0.8,0.75,0.95,0.95);
  legAll->SetTextFont(132);
  legAll->SetTextSize(0.07);
  legAll->SetFillColor(0);
  legAll->SetFillStyle(0);
  legAll->SetBorderSize(0);
  
  
  // Add to legend
  legAll->AddEntry(gThreeNu, "3-#nu Osc.", "l");
  legAll->AddEntry(gFourNu, "4-#nu Osc.", "l");
  
  legAll->Draw("same");
  
  outputTitle = "GlobalFluxDeficit";
  outputTitle.ReplaceAll(".list","");
  c2->SaveAs(outputTitle + ".png");
  c2->SaveAs(outputTitle + ".eps");
  c2->SaveAs(outputTitle + ".pdf");
  c2->SaveAs(outputTitle + ".C");
  
  // Make a version with PROSPECT Phase I coverage highlighted

  // Detector position
  double aDetectorX[] = {9.5, 18};
  double aDetectorEX[] = {2.5, 2};
  // Y-Range 1.3 - 0.6
  double aY[] = {(1.3 + 0.6)/2, (1.3 + 0.6)/2};
  double aEY[] = {(1.3 - 0.6)/2, (1.3 - 0.6)/2};
  
  hr2->Draw();
  
  // Draw the 1.0 line
  lNullOsc->Draw();
  
  // Draw the oscillated lines
  gThreeNu->Draw("CSAME");
  gFourNu->Draw("CSAME");
  
  gAllExps->Draw("PSAME");
  
  TGraphErrors* gPhaseI = new TGraphErrors(1, aDetectorX, aY, aDetectorEX, aEY);
  gPhaseI->SetName("PhaseI");
  gPhaseI->SetFillStyle(3005);
  gPhaseI->SetFillColor(kMagenta+3);
  
  gPhaseI->Draw("same2");
  
  TLegend* legI = new TLegend(0.8,0.67,0.95,0.95);
  legI->SetTextFont(132);
  legI->SetTextSize(0.07);
  legI->SetFillColor(0);
  legI->SetFillStyle(0);
  legI->SetBorderSize(0);
  
  
  // Add to legend
  legI->AddEntry(gThreeNu, "3-#nu Osc.", "l");
  legI->AddEntry(gFourNu, "4-#nu Osc.", "l");
  legI->AddEntry(gPhaseI, "Phase I", "f");
  
  legI->Draw("same");
  
  outputTitle = "GlobalFluxDeficitPhaseI";
  outputTitle.ReplaceAll(".list","");
  c2->SaveAs(outputTitle + ".png");
  c2->SaveAs(outputTitle + ".eps");
  c2->SaveAs(outputTitle + ".pdf");
  c2->SaveAs(outputTitle + ".C");
  
  // Make a version with PROSPECT Phase II coverage highlighted
  hr2->Draw();
  
  // Draw the 1.0 line
  lNullOsc->Draw();
  
  // Draw the oscillated lines
  gThreeNu->Draw("CSAME");
  gFourNu->Draw("CSAME");
  
  gAllExps->Draw("PSAME");
  
  // Phase II
  TGraphErrors* gPhaseII = new TGraphErrors(2, aDetectorX, aY, aDetectorEX, aEY);
  gPhaseII->SetName("PhaseII");
  gPhaseII->SetFillStyle(3005);
  gPhaseII->SetFillColor(kMagenta+3);
  
  gPhaseII->Draw("same2");

  TLegend* legII = new TLegend(0.8,0.67,0.95,0.95);
  legII->SetTextFont(132);
  legII->SetTextSize(0.07);
  legII->SetFillColor(0);
  legII->SetFillStyle(0);
  legII->SetBorderSize(0);
  
  
  // Add to legend
  legII->AddEntry(gThreeNu, "3-#nu Osc.", "l");
  legII->AddEntry(gFourNu, "4-#nu Osc.", "l");
  legII->AddEntry(gPhaseII, "Phase II", "f");
  
  legII->Draw("same");

  outputTitle = "GlobalFluxDeficitPhaseII";
  outputTitle.ReplaceAll(".list","");
  c2->SaveAs(outputTitle + ".png");
  c2->SaveAs(outputTitle + ".eps");
  c2->SaveAs(outputTitle + ".pdf");
  c2->SaveAs(outputTitle + ".C");
}
