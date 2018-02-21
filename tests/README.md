#### This ($OSCSENSFITTER/tests)  contains code to run automatic tests

These are the minimum tests to be performed for any major changes in the code:
1. Generate sensitivity curve and see if there are any variations from the previous versions that are not expected
2. Check if the output matches with inout in case of modeled oscillations (see 'Oscillation Fitting' below)
3. Check if all bins LvsErelative plots are in general close to 1.

Before performing these tests make sure that the environment variable OSCSENSFITTER is set ot the directory OscSens_CovMatrix
$ export OSCSENSFITTER=locaiton_of_this_directory a.k.a ./

Oscillation Fitting:
Run oscillation fitting with a few combinations of oscillation input parameters. To perform this:

$ cd $OSCSENSFITTER/tests
$ python3 generateBFFiles.py
$ bash runBFTests.sh
$ cd $OSCSENSFITTER
$ make testOscillationFits
$ ./testOscillationFits $OSCSENSFITTER/tests
