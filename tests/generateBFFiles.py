#!/usr/bin/python

import os
 
dm2List = [0.2,0.5,1.78,10]
s22List = [0.01, 0.09,0.3,0.8]

print('Creating %s files'%(len(dm2List)*len(s22List)))

for i in range(len(dm2List)):
  for j in range(len(s22List)):
    OSCSENSFITTER = os.environ['OSCSENSFITTER']
    MACROFILENAME='%s/tests/BF_%s_%s.mac' % (OSCSENSFITTER,dm2List[i],s22List[j])
    f = open(MACROFILENAME,'w')
    f.write('Detector101 = Yes\n')
    f.write('RxOnCycles101 =7\n')
    f.write('RxOffCycles101 =7\n')
    f.write('doRelativeMinimization = Yes\n')
    #Make sure this the location of extracted data is right
    f.write('EnergySpectrumFile = %s/inputs/EnergyTableHeu.txt\n' % OSCSENSFITTER)
    f.write('BKGInputFile = %s/inputs/SimulatedBackground.txt\n' % OSCSENSFITTER)
    f.write('DetectorResponseFile101 = %s/inputs/DetectorResponse.root\n' % OSCSENSFITTER)
    f.write('DataRootFile101 = %s/untracked/Truth/ExtractedData.root\n' % OSCSENSFITTER)
    f.write('IsReferenceOscillated = Yes\n')
    REFDM2 = "ReferenceDeltam2 = %s\n" %(dm2List[i])
    REFS22T = "ReferenceSin22theta = %s\n" %(s22List[j])
    f.write(REFDM2)
    f.write(REFS22T)
    MINFILENAME = "MinimizationFileName = %s/tests/BF_%s_%s.root\n" % (OSCSENSFITTER,dm2List[i],s22List[j])
    f.write(MINFILENAME)
    print("Generated %s" % MACROFILENAME)
    f.close()
print("DONE!")

