#!/bin/bash

# Test if config file exists
if ! [[ -e config.lua ]]
then
   echo "Please create config file 'config.lua' in the test directory and run again."
   echo "The 'config.lua' may use 'expandconfig' syntax"
   exit 1
fi

# Expand config file
rm -f cfg_exp/*
cp config.lua cfg_exp
./expandconfig cfg_exp/config.lua
rm cfg_exp/config.lua

# Get current time stamp
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

# Take the first expanded config file and execute the program on it
OMP_NUM_THREADS=1 ../heavy_top cfg_exp/config.lua.0000000.exp out/$(echo ${TIMESTAMP}_try)
