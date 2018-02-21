!!!!!!!!!!!!!!!!!!!!!! This code it deprecated as of jun 14, 2017 !!!!!!!!!!!!!!!!!!!!!!

OscSens_CovMatrix

Initial Commit by Karin Gilje 2/8/2016

This code uses a covariance matrix approach to calculating the minimum chi2. The previous sensitivity code used a minimized pulls method. Many details are missing from this readme and a technical note will be produced in order to further explain the usage of the code and the method of the calculation.

— CreateExperiment.C: builds the reactor and detector and calculates random distances for use in the calculation. Run by root -l -b -q ‘CreateExperiment.C+(dettype, exposure)’. Dettype is 0-front near detector, 1-middle near detector, 2-back near detector, 3-Far detector. And exposure is detector running time in years. Creates a file called Setup_HFIR_DetectorSetup_Exposure.root

— SimulateOscillation.C: Simulates the LvsE distribution for 57 values of Delta m^2. Also does an energy smearing of the spectrum. Run by root -l -b -q Setup_HFIR_DetectorSetup_Exposure.root SimulateOscillation.C+. Creates a file called Oscillation_HFIR_DetectorSetup_Exposure.root

— RunMinimizationSingleDetector.C: Runs a chi-square minimization using the covariance matrix approach (no Minuit minimization needed). Run by root -l -b -q Oscillation_HFIR_DetectorSetup_Exposure.root ‘RunMinimization.C+(RxOnFraction)’. Where RxOnFraction (41%) is the live time of the reactor and RxOffFraction is the fraction of time cosmic runs are taken to understand the backgrounds.

— RunPullsMinimization.C: A streamlined version of the minimized pulls approach to compare to the covaraince matrix approach.

Other files/directories:

CovarianceSetup: This directory contains the code to create covariance matrices for input into the chi^2 minimization calculation. inputs: This directory contains the inputs needed for the Kopp and LSN contours.
