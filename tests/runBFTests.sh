#!/bin/bash

for FILE in $OSCSENSFITTER/tests/*mac
do
# Run parallely all the fit files
  $OSCSENSFITTER/runOscSensChain $FILE &
  wait
  echo "Process finished"
done
