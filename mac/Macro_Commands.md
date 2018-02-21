#####################################################################
#### Author: P T Surukuchi
#### Date: Oct, 2016
#### Example macro with all the possible option so far included 
#### Make sure to use correct datatypes for each variable 
#### Comments starts with a #
#### For cases where there are different inputs for different detectors, the key is succeeded by the detector code or the detector type
#####################################################################

### Add detectors
### Numbering scheme as follows: 
### detector_type*100 +detector_location
### e.g near detector at near location is 1*100+1 =101
### Detector101 = Yes, adds the detector 101
### Detector102 = No, doesn't add the detector 102 
Detector101 = Yes
Detector102 = No
Detector103 = No

### The package prioritizes seconds over days over cycles
### Data-taking time in seconds during 
RxOnRxOnSeconds101 = 2134080
RxOnSeconds102 = 2134080
RxOnSeconds103 = 2134080 

### Data-taking time in days during RxOn
RxOnDays101 = 24.7
RxOnDays102 = 24.7
RxOnDays103 = 24.7 

### Number of RxOn cycles( 1 cycle = 24.7 days) for each detector
RxOnCycles101 = 1
RxOnCycles102 = 1
RxOnCycles103 = 1 

### Data-taking time in seconds during RxOff
#RxOffSeconds101 = 2134080
#RxOffSeconds102 = 2134080
#RxOffSeconds103 = 2134080
### Data-taking time in days during RxOn
RxOffDays101 = 24.7
RxOffDays102 = 24.7
RxOffDays103 = 24.7 
### Number of RxOff cycles( 1 cycle = 24.7 days) for each detector
RxOffCycles101 = 1
RxOffCycles102 = 1
RxOffCycles103 = 1 

#Radius, height in meters
ReactorRadius = 0.2
ReactorHeight = 0.5

#Power in GW
ReactorPower = 0.085

### Energy spectrum file input, the file location has to be wrt to the location of the makefile or location-independent
EnergySpectrumFile = ./inputs/EnergyTableHeu.txt

### Background file
BKGInputFile = ./inputs/SimulatedBackground.txt

### Detector Angular Offsets from the radial axis of the Reactor
AngleOffset101 = 4.9
AngleOffset102 = 4.9
AngleOffset103 = 9.43

### Number of segments
XSegments1 = 14
ZSegments1 = 11

### Segment width and length in meters
SegmentWidth1 = 0.146
SegmentLength1 = 1.17

### Efficiency in fraction
Efficiency1 = 0.436

### Position resolution along to the segment axis in m
YPositionResolution1 =0.07

### Energy resolution in 1/sqrt(E)
EnergyResolution1 =0.045

### Number of outer segments to calculate fiducial segments
XOuterSegmentsLeft1 = 1
XOuterSegmentsRight1 = 1 
ZOuterSegmentsTop1 = 1
ZOuterSegmentsBottom1 = 1

### Position resolution perpendicular to the segment axis in m
#XPositionResolution1 = 2.9
#ZPositionResolution1 = 0.28

### If selected Yes, the reference model for calculating Chi2 will be an oscillated model
IsReferenceOscillated = No
#Input the user-defined oscillation parameters for reference model histograms
ReferenceDeltam2 = 1.78
ReferenceSin22theta = 0.09

### If yes, perform oscillation analysis, if no, perform sensitivity 
UseData = Yes

### If selected yes, the sensitivity/oscillation is generated using toys(instead of default model)
UseToys = Yes

### Number of toys
NToys = 500

### Number of events to start out with
#Events =1000

### Use an input file for deltam2 and sin22theta bins
InputOscBins = No

### The file name for user input for deltam2 and sin22theta bins
### The input file must end with .txt or .root
### If using root file, make sure to include a TH2D histogram with the name OscBins
OscBinsInputFile = ./mac/OscBins.txt
#OscBinsInputFile = ./mac/OscBins.root

### Number of deltam2 sin22theta bins
NDeltam2 = 57
NSinSq2Theta = 49

### Uncertainties to be included in the covariance analysis
### Currently the options are     
##    SSTAT=1, // Signal statisical covariance matrices
##    BSTAT, // Background statistical covariance matrix
##    SNORM, // Signal normalization matrix
##    BNORM, // Background normalization matrix
##    SSHAPE, // Signal shape covariance matrix
##    BSHAPE, // Background shape covariance matrix
##    SESCALE, // Signal energy scale covariance matrix
##    SB2B,// Signal bin to bin covariance matrix
MinimizationUncertainties101 = {SStat}
ToyUncertainties101 = {SStat}

### Files corresponding to a particular signal covariance matrix
### Do not include files for signal and background stat
SignalNormCovFile101 = covInputs/Cov_Bkg_HFIR_102_UC.root

### Files corresponding to a particular toy covariance matrix
##ToyNormCovFile101 = covInputs/Cov_Bkg_HFIR_102_UC.root

### Output setup file name
SetupFileName = Setup.root
### Output oscillation file name
OscillationFileName = Oscillation.root
### Output minimization file name
MinimizationFileName = Minimization.root