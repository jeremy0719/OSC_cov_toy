#include "TStyle.h"
#include "TColor.h"

void SetGlobalStyle() {
  gStyle->SetPadTopMargin(0.16);
  gStyle->SetPadRightMargin(0.16);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);
  
  //Set Font
  gStyle->SetTextFont(132);
  gStyle->SetTextSize(0.08);
  gStyle->SetLabelFont(132,"x");
  gStyle->SetLabelFont(132,"y");
  gStyle->SetLabelFont(132,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.06,"z");
  gStyle->SetLabelFont(132,"t");
  gStyle->SetTitleFont(132,"x");
  gStyle->SetTitleFont(132,"y");
  gStyle->SetTitleFont(132,"z");
  gStyle->SetTitleFont(132,"t"); 
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleX(0.50);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(132,"pad");
  gStyle->SetStatFont(132);
  gStyle->SetStatFontSize(0.08);
  
  gStyle->SetHatchesLineWidth(2);
  
  // do not display any of the standard histogram decorations
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);
  
  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Add a greyscale palette for 2D plots
  const Int_t NRGBs = 5;
  const Int_t NCont = 200;

  Double_t stops[NRGBs] = { 0.00, 0.45, 0.50, 0.55, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 1.00, 0.52, 0.52*0.5 };
  Double_t green[NRGBs] = { 0.52*0.5, 0.52, 1.00, 0.00, 0.00 };
  Double_t blue[NRGBs]  = { 0.35*0.5, 0.35, 1.00, 0.52, 0.52*0.5 };

  // Purple . . . 128, 20, 128
  // Blue   . . . 20, 128, 92

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                   NCont);
  gStyle->SetNumberContours(NCont);

}

void SetTwoColorPalette() {
  const Int_t NRGBs = 5;
  const Int_t NCont = 200;
  Double_t stops[NRGBs] = { 0.00, 0.45, 0.50, 0.55, 1.00};
  Double_t red[NRGBs]   = { 0.00, 0.00, 1.00, 0.52, 0.52*0.5};
  Double_t green[NRGBs] = { 0.52*0.5, 0.52, 1.00, 0.00, 0.00};
  Double_t blue[NRGBs]  = { 0.35*0.5, 0.35, 1.00, 0.52, 0.52*0.5};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  
  gStyle->SetNumberContours(NCont);
}

void SetBlackBodyPalette() {
  const Int_t nRGBs = 5;
  const Int_t nCont = 200;
  Double_t stops[nRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00};
  Double_t red[nRGBs]   = { 1.00, 0.90, 0.70, 0.50, 0.00};
  Double_t green[nRGBs] = { 1.00, 0.50, 0.25, 0.05, 0.00};
  Double_t blue[nRGBs]  = { 1.00, 0.90, 0.70, 0.50, 0.00};
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
  
  gStyle->SetNumberContours(nCont);
}

void SetGradientPalette() {
 int ncol=50;
 double dcol = 1./float(ncol+10);
 float hue = 300.0;
 float bright = 0.5;
 TColor **theCols = new TColor*[ncol];
 for (int i=0;i<ncol;i++) theCols[i] = new TColor(9999-i,0.0,0.7,0.7);
 for (int j = 0; j < ncol; j++) {
   float r = 0.0;
   float g = 0.0;
   float b = 0.0;
   float saturation = (ncol-j)*dcol;
   TColor::HLS2RGB(hue,saturation,bright,r,g,b);
   theCols[j]->SetRGB(r,g,b);
   //std::cout << hue << " " << bright <<  " " << saturation << " " << r <<  " " << g << " " << b <<std::endl;
 }
 int ColJul[100];
 for  (int i=0; i<100; i++) ColJul[i]=9999-i;
 gStyle->SetPalette(ncol,ColJul);
}