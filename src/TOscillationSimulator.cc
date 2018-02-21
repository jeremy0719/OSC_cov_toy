#include "TOscillationSimulator.hh"

TOscillationSimulator::TOscillationSimulator(){}

TOscillationSimulator::~TOscillationSimulator(){}

double TOscillationSimulator::GetDeltam2Term(double deltam2, double energy, double baseline)
{
  return TMath::Power(TMath::Sin(1.27 * deltam2 * baseline / energy), 2);
}

double TOscillationSimulator::ComputeNueDisProb(double deltam2, double sinSq2Theta, double energy, double baseline)
{
  return sinSq2Theta*GetDeltam2Term(deltam2, energy, baseline);
}

double TOscillationSimulator::ComputeNueSurvProb(double deltam2, double sinSq2Theta, double energy, double baseline)
{
  return 1-ComputeNueDisProb(deltam2, sinSq2Theta, energy, baseline);
}

// CAUTION:Something weird going on with this function, use with care
TH2D* TOscillationSimulator::Oscillate2DHistogram(const TH2D* hRef,double deltam2, double sinSq2Theta)
{
  TH2D* hOsc=(TH2D*)hRef->Clone();
  hOsc->Reset();
  for(int i=0;i<hOsc->GetNbinsX();i++)
  {
    double energy = hRef->GetXaxis()->GetBinCenter(i+1);
    for(int j=0;j<hOsc->GetNbinsY();j++)
    {
      double baseline= hRef->GetYaxis()->GetBinCenter(j+1);
    hOsc->Fill(energy,baseline,ComputeNueSurvProb(deltam2,sinSq2Theta,energy,baseline));
    }
  }
  return hOsc;
}
