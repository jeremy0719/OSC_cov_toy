#include "TObject.h"
#include "TMinuit.h"
#include "TH2.h"
#include <iostream>

//#ifndef CHISQUAREDAYABAY_H
//#define CHISQUAREDAYABAY_H 1

class ChiSquareDayaBay : public TObject{
private:
  double m_sinSq2Theta;
  double m_sigmaSource;
  int    m_nEbins;
  int    m_nRbins;
  int    m_nEbinsTot;
  int    m_nRbinsTot;
  int    m_minEbin;
  int    m_maxEbin;
  int    m_minRbin;
  int    m_maxRbin;
  double m_sigmaEbin;
  double m_eErr[100];
  double m_rErr;
  double m_penaltySource;
  double m_penaltyEbins[100];
  TH2*   m_signalNoOsc;
  TH2*   m_signalTestOsc;
  TH2*   m_signalMaxOsc;
  TH2*   m_signalOsc;
  TH2*   m_background;
  TMinuit* m_minuit;
  bool   m_needsUpdate;
  double m_minimum;
  int    m_errorFlag;
  double m_errors[100];
public:
  ChiSquareDayaBay();
  virtual ~ChiSquareDayaBay();
  double Eval();
  double EvalStat();
  int    FindMinimum();
  double Minimum();
  double Nevents(TH2 *hist);
  int    Update();
  void   SetSinSq2Theta(double sinSq2Theta);
  void   SetBounds(TH2 *hist);
  void   SetSigmaSource(double sigma);
  void   SetSigmaEbin();
  void   SetPenaltySource(double alpha);
  void   SetPenaltyEbins(double* alphas);
  void   SetSignalNoOscillation(TH2* hist);
  void   SetSignalTestOscillation(TH2* hist);
  void   SetSignalMaxOscillation(TH2* hist);
  void   SetBackground(TH2* hist);
  //ClassDef(ChiSquareDayaBay,1)
};

//#endif // CHISQUAREDAYABAY_H

////////////////////////////////////////////////////////////////////

void minfcn(Int_t &/*nParameters*/, Double_t */*gin*/, Double_t &value, 
	    Double_t *parameters, Int_t /*iflag*/)
{
  // Silly static function needed for Minuit
  ChiSquareDayaBay* csdyb =
    dynamic_cast<ChiSquareDayaBay*>(gMinuit->GetObjectFit());

  csdyb->SetPenaltySource(parameters[0]);
  csdyb->SetPenaltyEbins(&parameters[1]);

  value = csdyb->Eval();
  return;
}


//////////////////////////////////////////////////////////////////

//ClassImp(ChiSquareDayaBay)

ChiSquareDayaBay::ChiSquareDayaBay()
: TObject(),
  m_sinSq2Theta(0),
  m_sigmaSource(0),
  m_sigmaEbin(0),
  m_nEbinsTot(0),
  m_nRbinsTot(0),
  m_minEbin(0),
  m_maxEbin(0),
  m_minRbin(0),
  m_maxRbin(0),
  m_penaltySource(0),
  m_signalNoOsc(0),
  m_signalTestOsc(0),
  m_signalMaxOsc(0),
  m_signalOsc(0),
  m_background(0),
  m_minuit(0),
  m_needsUpdate(true),
  m_minimum(-1),
  m_errorFlag(0)
{
}


ChiSquareDayaBay::~ChiSquareDayaBay()
{
  if(m_minuit){
    delete m_minuit;
    m_minuit=0;
  }
}

double ChiSquareDayaBay::Nevents(TH2* hist)
{
  std::cout<<"Events: "<<hist->Integral(m_minEbin,m_maxEbin,m_minRbin,m_maxRbin)<<std::endl;
  return hist->Integral(m_minEbin,m_maxEbin,m_minRbin,m_maxRbin);
}

void ChiSquareDayaBay::SetBounds(TH2* hist)
{
  m_minEbin = hist->FindFirstBinAbove(10,1);
  m_maxEbin = hist->FindLastBinAbove(10,1);
  m_minRbin = hist->FindFirstBinAbove(10,2);
  m_maxRbin = hist->FindLastBinAbove(10,2);

  m_nEbinsTot = hist->GetNbinsX();
  m_nRbinsTot = hist->GetNbinsY();
  m_nEbins=m_maxEbin-m_minEbin+1;
  m_nRbins=m_maxRbin-m_minRbin+1;

  std::cout<<"E Threshold: "<<hist->GetXaxis()->GetBinLowEdge(m_minEbin)<<std::endl;
  std::cout<<"MinimumCenter: "<<hist->GetXaxis()->GetBinCenter(m_minEbin)<<std::endl;
  std::cout<<"MaximumCenter: "<<hist->GetXaxis()->GetBinCenter(m_maxEbin)<<std::endl;
  std::cout<<"MinimumCenter: "<<hist->GetYaxis()->GetBinCenter(m_minRbin)<<std::endl;
  std::cout<<"MaximumCenter: "<<hist->GetYaxis()->GetBinCenter(m_maxRbin)<<std::endl;
  std::cout<<"MinE "<<m_minEbin<<" of "<<hist->GetNbinsX()<<std::endl;
  std::cout<<"MaxE "<<m_maxEbin<<" of "<<hist->GetNbinsX()<<std::endl;
  std::cout<<"MinR "<<m_minRbin<<" of "<<hist->GetNbinsY()<<std::endl;
  std::cout<<"MaxR "<<m_maxRbin<<" of "<<hist->GetNbinsY()<<std::endl;
    
}

void ChiSquareDayaBay::SetSinSq2Theta(double sinSq2Theta)
{
  m_sinSq2Theta = sinSq2Theta;
  m_needsUpdate=true;
}

void ChiSquareDayaBay::SetSigmaSource(double sigma)
{
  m_sigmaSource=sigma;
}

void ChiSquareDayaBay::SetSigmaEbin()
{
  
  double errors[m_nEbins];
  
  for (int i = 0; i < m_nEbins; i++) errors[i] = 0.1;
  
  
  for(int i=0;i<m_nEbins;i++){
    m_eErr[i]=errors[i];
  }
}


void ChiSquareDayaBay::SetPenaltySource(double alpha)
{
  m_penaltySource=alpha;
}

void ChiSquareDayaBay::SetPenaltyEbins(double* alphas)
{
  for(int eIdx=0; eIdx<m_nEbins; eIdx++){
    m_penaltyEbins[eIdx]=alphas[eIdx];
  }
}

void ChiSquareDayaBay::SetSignalNoOscillation(TH2* hist)
{
  m_signalNoOsc=hist;
  m_needsUpdate=true;
}

void ChiSquareDayaBay::SetSignalTestOscillation(TH2* hist)
{
  m_signalTestOsc=hist;
  m_needsUpdate=true;
}

void ChiSquareDayaBay::SetSignalMaxOscillation(TH2* hist)
{
  m_signalMaxOsc=hist;
  m_needsUpdate=true;
}

void ChiSquareDayaBay::SetBackground(TH2* hist)
{
  m_background=hist;
}

int ChiSquareDayaBay::Update()
{
  for(int i=0; i<m_nEbins; i++){ m_penaltyEbins[i]=0; }
  for(int i=0; i<(1); i++){ m_errors[i]=0; }

  // Initialize minimizer if needed
  if(!m_minuit){
    m_minuit = new TMinuit(60);//was 60
    m_minuit->SetPrintLevel(-1);
    m_minuit->SetFCN(minfcn);
    m_minuit->SetObjectFit(this);

    m_minuit->mnparm(0,"alphaSource",0,m_sigmaSource,0,0,m_errorFlag);

    for(int eIdx=0; eIdx<m_nEbins; eIdx++){
      m_minuit->mnparm(eIdx+1,Form("alphaE%d",1),
                       0,m_eErr[eIdx],0,0,m_errorFlag);
    }
  }

  if(m_signalOsc){
    delete m_signalOsc;
    m_signalOsc=0;
  }
  m_signalOsc = dynamic_cast<TH2*>(m_signalNoOsc->Clone("signalOsc"));
  m_signalOsc->Add(m_signalMaxOsc,-1.0*m_sinSq2Theta);
  //std::cout<<"EventsOsc: "<<m_signalOsc->Integral(m_minEbin,m_maxEbin,m_minRbin,m_maxRbin)<<std::endl;
  //std::cout<<"EventsUnosc: "<<m_signalNoOsc->Integral(m_minEbin,m_maxEbin,m_minRbin,m_maxRbin)<<std::endl;
  
  m_minimum=-1;
  m_needsUpdate=false;
  return 0;
}

double ChiSquareDayaBay::EvalStat()
{
  // Evaluate the Chi-Square using only the statistical errors.  No minimization needed!
  // Calculate the chi-square, looping over bins in the L vs E event
  // distribution
  double chiSq = 0;
  
  // Make sure that m_signalOsc exists!
  if(m_signalOsc){
    delete m_signalOsc;
    m_signalOsc=0;
  }
  m_signalOsc = dynamic_cast<TH2*>(m_signalNoOsc->Clone("signalOsc"));
  m_signalOsc->Add(m_signalMaxOsc,-1.0*m_sinSq2Theta);
  
  for(int binE=m_minEbin; binE<=m_maxEbin; binE++){
    for(int binR=m_minRbin; binR<=m_maxRbin; binR++){
      
      // Expectation
      double nExp = m_signalTestOsc->GetBinContent(binE,binR);

      // Observation
      double nObs = m_signalOsc->GetBinContent(binE,binR);
      
      // Background
      double nBG = 0;
      if(m_background) nBG = m_background->GetBinContent(binE,binR);
      if (nBG == 0) continue;

      // Total Expectation
      double nTotal = nExp + nBG;
      if (nTotal == 0) continue;
      
      // Difference between Observed and Expected
      double dN = (nObs+nBG) - nTotal;

      // Chi Square Calculation
      chiSq += (dN*dN)/nTotal;

    }
  }
  return chiSq;
}

double ChiSquareDayaBay::Eval()
{
  // Evaluate the Chi-Square using the currently set parameter values

  if(m_needsUpdate) this->Update();

  // Calculate the chi-square, looping over bins in the L vs E event
  // distribution
  double chiSq = 0;

  for(int binE=m_minEbin; binE<=m_maxEbin; binE++){
    for(int binR=m_minRbin; binR<=m_maxRbin; binR++){

      double nExp = m_signalTestOsc->GetBinContent(binE,binR);

      double nObs = m_signalOsc->GetBinContent(binE,binR);
      double nBG = 0;
      if(m_background) nBG = m_background->GetBinContent(binE,binR);
      if (nBG == 0) continue;

      double penaltyEbin=m_penaltyEbins[binE-m_minEbin];
      
      double nTotal = (1+m_penaltySource+penaltyEbin)*nExp + nBG;

      double dN = (nObs+nBG) - nTotal;
      
      if(nTotal<=0){
        if((nObs+nBG)>0){
          std::cout << "Non-zero observed events when 0 expected!" << std::endl;
          std::cout << "Total: " << nTotal << " nExp: " << nExp << " nObs: " << nObs << " nBG: " << nBG << std::endl;
          std::cout <<  " " << m_penaltySource << std::endl;
          std::cout << " sin2 2theta = " << m_sinSq2Theta << std::endl;
          std::cout << " E bin = " << binE << std::endl;
          std::cout << " R bin = " << binR << std::endl;

          m_penaltySource=0.0;
        }
        std::cout << "Negative events expected..." << m_penaltySource << std::endl;
        continue;
      }
      
      chiSq += (dN*dN)/(nExp+nBG);

    }

  }

  // Add Penalty terms
  chiSq += (m_penaltySource*m_penaltySource)/(m_sigmaSource*m_sigmaSource);
  
  for(int eIdx=0;eIdx<m_nEbins;eIdx++){
    chiSq += (m_penaltyEbins[eIdx]*m_penaltyEbins[eIdx])/(m_eErr[eIdx]*m_eErr[eIdx]);
  }
  
  return chiSq;
}

int ChiSquareDayaBay::FindMinimum()
{
  // Find the minimum chi-square, allowing the penalty parameters to vary
  if(m_needsUpdate) this->Update();

  for(int i=0; i<m_nEbins; i++){ m_penaltyEbins[i]=0; }
  for(int i=0; i<(1); i++){ m_errors[i]=0; }
  
  m_minuit->Migrad();

  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  m_minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  m_minimum = amin;
  // Catch parameters and errors
  m_minuit->GetParameter(0,m_penaltySource,m_errors[0]);
  
  for(int eIdx=0;eIdx<m_nEbins;eIdx++){
    m_minuit->GetParameter(1+eIdx,m_penaltyEbins[eIdx],m_errors[eIdx]);
  }
  
  return 0;
}

double ChiSquareDayaBay::Minimum()
{
  // Return the last-found minumum of the chi-square
  return m_minimum;
}
