### OscSensFitterCC
Package for oscillation fitting and sensitivity estimation for PROSPECT experiment.
#### Author: Pranava Teja Surukuchi
#### Date: Oct, 2017

#### Dependencies:

GNU Make (https://www.gnu.org/software/make/)    
c++11    
ROOT (https://root.cern.ch/) v5.34 and up    
#### Building:    

Define environment variable for the location of the directory
$ export OSCSENSFITTER=locaiton_of_this_directory a.k.a ./
Make all the objects and library needed for oscillation fitting package    
$ make    
Make an example executable(with the same name as the .cc file) and link with the objects and library    
$ make runOscSensChain    
Run the example executable    
$ ./runOscSensChain    
for default macro or,    
$ ./runOscSensChain path_to_macro_file    
for a specific macro    


#### Code Structure:    
To make the code more structured and easy to understand, the following structure has been imposed:    

./ (Home directory: Contains Makefile and all the executables. Executables can be moved to a separate location if needed.)    
./src (Source directory: Contains all the .cc and .hh files that are used to generate objects)    
./obj (Object directory: Contains all .o files)    
./lib (Library directory: Contains library generated from the package, most likely only contains one library so probably not needed)    
./inputs (Input directory: Contains inputs needed for the package)    
./covInputs (Covarince input directory: Contains covariance root file inputs needed for the package)    
./doc (Documentation directory: Doxygen output will be generated here in html and latex format)    
./mac (Macro Directory: Directory to store macro files)    
./aux (Directory with auxiliary files that are used for generating various plots including supplemental plots the oscillation sensitivity plots)    
./untracked (Untracked git Directory: Directory for personal files, not tracked by git)    
./aux/ToyCovExecs Files used to generate toy-based covariance matrices
./tests/ Files used to autorun tests after major changes to the code
#### Documentation:    

$ doxygen Osc_Sens_Doxygen    
The documentation in both latex and html format will be generated in ./doc/html and ./doc/latex respectively    
#### Coding Guidlines:    
The document describes the coding guidelines followed in this package can be found here.    

#### Macros:    

User can provide inputs from a macro file ending with '.mac'    
Each line in the macro file is a esentially a key-value pair delimited by '='    
Find the example macro will all option included in ./mac/sensitivity.mac    
Make sure to use correct datatypes for each variable    
Comments starts with a'#'    
For cases where there are different inputs for different detectors, the key is succeeded by the detector code or the detector type    
To keep track of all changes, any new macro commands have to be included in the ./mac/default.mac file    


#### Generating Covariance Matrices

$ make generateSignalNormCovMatrix
Generate signal covariance matrix

#### Extracting data from calibrated files
Make executable for data extraction
$ make extractData
Extract data from calibrated files to act as inputs
$ ./extractData
for default extraction macro (mac/extract.mac), or
$ ./extractData path_to_macro_file    
for a specific macro    

