#include <fstream>

#include <TMath.h>

// The number of generated fake Events
int fNEvents = 1000000;

// The span of time used for detector data taking
double fExposure = 1.0;

//////////////////////////////////////////////////////////
// Reactor Values
//////////////////////////////////////////////////////////

// The HFIR reactor core size
double fReactorRadius = 0.2; // meters
double fReactorHeight = 0.6; // meters

// The HFIR reactor core power
double fReactorPower = 85.0 / 1000.0; // MW turned into GW

// The yearly HFIR reactor Live Time fraction
// Changed the reactor-on from 6 cycles to 7 cycles /yr
double fRxOnFraction = 0.474; // %

// The yearly cosmic data taking fraction with the HFIR reactor off.
// Note that this cannot be larger than 1 - fRxOnFraction
// Also note that we most likely will not be able to take reactor-off data 50% time in an year
double fRxOffFraction = 0.5;

//////////////////////////////////////////////////////////
// Detector Values
//////////////////////////////////////////////////////////

// The size of a segment cell
double fSegWidth = 0.146; // meters

// The distance from the reactor to the front floor edge of the near detector (front position)
// Originally based on reactor center to detector center (14x10) and 14.6 cm segments engineering drawing.
// (5.904-(7*0.146))
double fNearDistX = 4.882; // meters
double fNearDistY = 0.0; // meters
double fNearDistZ = 4.28; // meters

// The number of segments in the near detector;
double fNearSegX = 14.0;
double fNearSegZ = 11.0;

// The length of the near segments.
// Length changed from 1.19 to 1.17
double fNearSegLength = 1.17; // meters

// The distance from the reactor to the far detector front floor edge
// fFarDistZ is originally based on a 28x20 detector.
double fFarDistX = 15.8; // meters
double fFarDistY = 0.0; // meters
double fFarDistZ = 0.9 - 10 * 0.146; // meters

// The number of segments in the far detector
double fFarSegX = 28.0;
double fFarSegZ = 20.0;

// The length of the far detector segments
double fFarSegLength = 1.17; // meters

// A summary of the Detector positons
double fDistX = fNearDistX;
double fDistY = fNearDistY;
double fDistZ = fNearDistZ;

// A summary of the Detector Segments
double fSegX = fNearSegX;
double fSegZ = fNearSegZ;
double fNSeg = fSegX * fSegZ;
double fNFiducialSeg = (fSegX-2) * (fSegZ-2);
double fSegLength = 1.17; // meters

// This function will set the detector type variables
void SetDetectorVariables(int DetectorType) {
  
  // Set the near detector values
  if (DetectorType < 3) {
    fDistX = fNearDistX + DetectorType * 1.5;
    fDistY = fNearDistY;
    fDistZ = fNearDistZ;
    fSegX = fNearSegX;
    fSegZ = fNearSegZ;
  }
  // Set the far detector values
  else {
    fDistX = fFarDistX;
    fDistY = fFarDistY;
    fDistZ = fFarDistZ;
    fSegX = fFarSegX;
    fSegZ = fFarSegZ;
  }
  fNSeg = fSegX * fSegZ;
  fNFiducialSeg = (fSegX-2) * (fSegZ-2);

  return;
}

// Proton Density of the Scintillator (target material)
double fProtonDensity = 5.46e28; // EJ-309. Detector Protons per volume (m^-3)
//double fProtonDensity = 6.39e28; // LAB. Detector Protons per volume (m^-3)

// Predicted Fiducialized Detector Efficiency
double fDetectorEfficiency = 0.436;

// Predicted Detector Position Resolution
// Segment size changed from 14.6 cm  to 14 cm
double fPosResX = 0.146; // meters
double fPosResY = 0.07; // meters
double fPosResZ = 0.146; // meters

// Predicted Detector Energy Resolution
double fEneRes = 0.045; // 4.5%/sqrt(E)

//////////////////////////////////////////////////////////
// Vectors of Energy Spectrum information
//////////////////////////////////////////////////////////

// Create vector to fill with antinu energy
std::vector<double> fAntiNuEnergy;

// Create vector to fill with the energy dependent incoming antineutrino flux
std::vector<double> fAntiNuFlux;

// Create vector to fill with the energy dependent IBD cross section
std::vector<double> fIBDCrossSection;

// This function will load the energy spectrum information into a file.
void ReadEnergySpectrum() {
  std::ifstream EnergySpectrumFile;
  EnergySpectrumFile.open("EnergyTableHeu.txt");
  
  if (!EnergySpectrumFile.is_open()) {
    std::cout << "Energy Spectrum File unavailable, check location and try again!" << std::endl;
    std::exit(1);
  }
  
  // Clear vectors to prepare to fill with energy spectrum information
  fAntiNuEnergy.clear();
  fAntiNuFlux.clear();
  fIBDCrossSection.clear();
  
  // The vector will have a first element of the mean energy
  // of the energy bin in MeV, the second element will be the cross
  // section in the energy bin in cm^2 and the last element will be
  // the number of neutrinos per unit energy (Jx10^?).
  // The values come from Vogel and Engel: Neutrino electromagnetic
  // form factors.  Equation 4 using parameters from table 1 gives
  // the spectrum measurement.
  while(!EnergySpectrumFile.eof()){
    double energy = 0.0;
    double flux = 0.0;
    double xsec = 0.0;
    EnergySpectrumFile >> energy >> xsec >> flux;
    // Save only if Energy != 0
    if (energy != 0.0) {
      fAntiNuEnergy.push_back(energy);
      fAntiNuFlux.push_back(flux / 100.0 / 100.0); // Converting cm^2 into m^2
      fIBDCrossSection.push_back(xsec);
    }
  }
  EnergySpectrumFile.close();
  
  // Total number of energy points
  std::cout << "The total number of Energy Bins used is "
  << fAntiNuEnergy.size() << std::endl;
}

//////////////////////////////////////////////////////////
// Constants Used in weight calculation
//////////////////////////////////////////////////////////

// Conversion of Years to Seconds
double fYearToSeconds = 3.15569e7;

//////////////////////////////////////////////////////////
// Calculate Delta m2
//////////////////////////////////////////////////////////

double fNDeltam2 = 57;

std::vector<double> fDeltam2;
double fDeltam2Bins[57];

void ConstructDeltam2() {
  fDeltam2.clear();
  fDeltam2.push_back(0.0);
  for (int i = 1; i < fNDeltam2; i++) {
    fDeltam2.push_back(TMath::Power(10, -1.4+0.05*(i-1)));
    fDeltam2Bins[i-1] = TMath::Power(10, -1.425 + 0.05*(i-1));
  }
  fDeltam2Bins[(int)fNDeltam2-1] = TMath::Power(10, -1.425 + 0.05*(fNDeltam2-1));
}

//////////////////////////////////////////////////////////
// Calculate sin^2 2theta
//////////////////////////////////////////////////////////

double fNSinSq2Theta = 49;

std::vector<double> fSinSq2Theta;
double fSinSq2ThetaBins[50];

void ConstructSinSq2Theta() {
  fSinSq2Theta.clear();
  for (int i = 0; i < fNSinSq2Theta; i++) {
    fSinSq2Theta.push_back(TMath::Power(10, -2.4+0.05*(i-1)));
    fSinSq2ThetaBins[i] = TMath::Power(10, -2.425 + 0.05*(i-1));
  }
  fSinSq2ThetaBins[(int)fNSinSq2Theta] = TMath::Power(10, -2.425 + fNSinSq2Theta*0.05);
}
