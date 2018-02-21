#include <TStyle.h>
#include <TROOT.h>
#include <TColor.h>

void SetupProspectStyle(){
  cout << endl << "PROSPECT Experiment plot defaults." << endl;
  
  //
  // KAJ 08-DEC-2006
  // Based on a style from ATLAS,
  // which was based on a style file from BaBar
  // and then lightly seasoned with input from Liu Jianglai
  //
  
  //Creating new style "PROSPECT"
  TStyle *prospectStyle= new TStyle("PROSPECT","Prospect style");
  
  // use plain black on white colors
  Int_t icol=0;
  prospectStyle->SetFrameBorderMode(icol);
  prospectStyle->SetCanvasBorderMode(icol);
  prospectStyle->SetPadBorderMode(icol);
  prospectStyle->SetPadColor(icol);
  prospectStyle->SetCanvasColor(icol);
  prospectStyle->SetStatColor(icol);
  //prospectStyle->SetFillColor(icol);
  
  prospectStyle->SetPalette(1,0); // set a good color palette
  //prospectStyle->SetFillColor(0); // explicitly don't fill
  //Int_t colors[]={4,4,7,7,6,6,6,2,2,2,1};
  //   float levels[]={1,10,14,22,45,60,200,300,400,500,600,700};
  
  // set the paper & margin sizes
  prospectStyle->SetPaperSize(TStyle::kUSLetter);
  
  //Margins, give some additional space on the Left/Bottom where axis labels go
  prospectStyle->SetPadTopMargin(0.07);
  prospectStyle->SetPadLeftMargin(0.111);
  prospectStyle->SetPadRightMargin(0.15);
  prospectStyle->SetPadBottomMargin(0.111);
  
  // set master font
  Int_t font=132;      //times new roman, regular (13), high precision (2)
  //Int_t font=42;     //tnr, bold (4), high precision (2)
  
  //set master font size
  //Double_t tsize=0.11; //should be set between 0.03-0.05, is in units of "% of pad"
  //Double_t tsize=0.09; //should be set between 0.03-0.05, is in units of "% of pad"
  Double_t tsize=0.05; //should be set between 0.03-0.05, is in units of "% of pad"
  //Double_t tsize=0.065; //should be set between 0.03-0.05, is in units of "% of pad"
  
  // Text, Labels, Titles
  // All set the same here
  prospectStyle->SetTextFont(font);
  prospectStyle->SetTextSize(tsize);
  
  prospectStyle->SetLabelFont(font,"xyz");
  prospectStyle->SetLabelSize(tsize,"xyz");
  prospectStyle->SetLabelOffset(0.003,"xyz");
  
  prospectStyle->SetTitleFont(font,"xyz");
  prospectStyle->SetTitleSize(tsize,"xyz");
  //prospectStyle->SetTitleOffset(0.3,"z");
  //prospectStyle->SetTitleOffset(2.0,"y");
  //prospectStyle->SetTitleOffset(0.4,"y");
  prospectStyle->SetTitleOffset(0.9,"y");
  //prospectStyle->SetTitleOffset(1.4,"x");
  prospectStyle->SetTitleOffset(0.9,"x");
  prospectStyle->SetTitleBorderSize(0);
  
  //use bold lines and markers
  //prospectStyle->SetMarkerStyle(20);
  //prospectStyle->SetMarkerSize(1.2);
  prospectStyle->SetHistLineWidth(2.);
  prospectStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  
  //get rid of X error bars and y error bar caps
  //prospectStyle->SetErrorX(0.001);
  
  //do not display any of the standard histogram decorations
  prospectStyle->SetOptTitle(0);     //no title box
  //prospectStyle->SetOptStat(1); //standard stat info
  prospectStyle->SetOptStat(0);    //no stat info
  //prospectStyle->SetOptStat(10);     //number of entries only
  prospectStyle->SetOptFit(11);  //standard fit info
  //prospectStyle->SetOptFit(0);       //no fit info
  
  //legend and stat pad
  prospectStyle->SetStatBorderSize(1);
  prospectStyle->SetStatFont(132);
  prospectStyle->SetStatX(0.95);
  prospectStyle->SetStatY(0.95);
  prospectStyle->SetLegendBorderSize(1);
  
  // put tick marks on top and RHS of plots
  prospectStyle->SetPadTickX(1);
  prospectStyle->SetPadTickY(1);
  // tick marks only on left and bottom
  prospectStyle->SetPadTickX(0);
  prospectStyle->SetPadTickY(0);
  
  // additional series styles (Thanks, Jianglai)
  //   1 = solid
  //   2 = long dash (30 10)
  //   3 = dotted (4 8)
  //   4 = dash-dot (15 12 4 12)
  //   5 = short dash ( 15 15 )
  //   6 = dash-dot-dot   gStyle->SetLineStyleString(1,"[]");
  prospectStyle->SetLineStyleString(2,"[30 10]");
  prospectStyle->SetLineStyleString(3,"[4 8]");
  prospectStyle->SetLineStyleString(4,"[15 12 4 12]");
  prospectStyle->SetLineStyleString(5,"[15 15]");
  prospectStyle->SetLineStyleString(6,"[15 12 4 12 4 12]");
  prospectStyle->SetOptDate(0);
  prospectStyle->SetDateY(.98);
  prospectStyle->SetStripDecimals(kFALSE);
  
  //  prospectStyle->SetEndErrorSize(5.0); //make the end of error bar longer 
  
  
  gROOT->SetStyle("PROSPECT");
  gROOT->ForceStyle();    
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

