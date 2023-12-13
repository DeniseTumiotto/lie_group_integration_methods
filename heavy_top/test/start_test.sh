#!/bin/bash

if [[ $# -lt 1 ]]
then
   CONFIG=config.lua
else
   CONFIG=$1
fi

# Test if config file exists
if ! [[ -e $CONFIG ]]
then
   echo "Please create config file '$CONFIG' in the test directory and run again."
   echo "The '$CONFIG' may use 'expandconfig' syntax"
   exit 1
fi

# Maximum parallel processes
MAXPAR=5
echo "Testing with $MAXPAR processes in parallel. Abort if you don't have enough processors"

# Expand config file
rm -f cfg_exp/*
cp $CONFIG cfg_exp
./expandconfig cfg_exp/$CONFIG
rm cfg_exp/$CONFIG

# Loop over all files in directory of expanded config file
#for CONF in cfg_exp/*.exp
for CONF in $(ls cfg_exp/*.exp | shuf)
#in random order!
do
   # Get the number
   NUMBER=${CONF: -11:7}
   # Get current time stamp
   TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

   # Call the integrator
   echo Executing number $NUMBER at $TIMESTAMP
   # DEBUG
   #operf ../heavy_top $CONF out/$(echo ${TIMESTAMP}_${NUMBER}) &> out/$(echo ${TIMESTAMP}_${NUMBER}).err && echo  && echo Success:: number $NUMBER at $TIMESTAMP ; echo Done with number $NUMBER at $TIMESTAMP
   #valkyrie ../heavy_top $CONF out/$(echo ${TIMESTAMP}_${NUMBER}) &> out/$(echo ${TIMESTAMP}_${NUMBER}).err && echo  && echo Success:: number $NUMBER at $TIMESTAMP ; echo Done with number $NUMBER at $TIMESTAMP
   # GUBED

   #sem -j $MAXPAR  "OMP_NUM_THREADS=1 ../heavy_top $CONF out/$(echo ${TIMESTAMP}_${NUMBER}) &> out/$(echo ${TIMESTAMP}_${NUMBER}).err && echo  && echo Success:: number $NUMBER at $TIMESTAMP ; echo Done with number $NUMBER at $TIMESTAMP"

   OMP_NUM_THREADS=1 ../heavy_top $CONF out/$(echo ${TIMESTAMP}_${NUMBER}) &> out/$(echo ${TIMESTAMP}_${NUMBER}).err
   echo
   echo Success:: number $NUMBER at $TIMESTAMP
   echo Done with number $NUMBER at $TIMESTAMP

done
