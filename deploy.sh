#!/bin/bash

# Marco Pizzichemi 11.1.2022 marco.pizzichemi@cern.ch

# This script compiles the simulation toolkit and all the support programs
# Run it from main git directory, with 
#
# ./deploy.sh [OPTION]
#
# Argument OPTION is optional. If not given, machine compiler is used.
# The only other possible option is lxplus
#
# ./deploy.sh lxplus
#
# which will source LGC97python3 environment. 
# If you want to add other options, you need to properly modify the SET ENV VARIABLES section

# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "Last command \"${last_command}\" returned with exit code $?."' EXIT


### MAKE BUILD FOLDER (IF NOT THERE ALREADY)
mkdir -p build
CMAKE_ARGS=""

if [ -z "$1" ]
  then
    echo "No argument supplied, so compiling with machine compiler..."
else
  if [ $1 = "lxplus" ]; then
    set --
    echo "Sourcing compiler for lxplus..."
    source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
    CMAKE_ARGS="--verbose -DCMAKE_CXX_COMPILER=`which g++` -DCMAKE_C_COMPILER=`which gcc`"
  else
    echo "Invalid argument $1 - You can either provide no arg (and the script will use machine compiler), or lxplus, and the script will source the lxplus variables"
  fi
fi

### SUPPORT PROGRAMS
echo "Compiling support programs..."
g++ -o build/timeResolution timeResolution.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas
g++ -o build/extractConfiguration extractConfiguration.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore
g++ -o build/calcDOIres DeepLearningDOI/calcDOIres.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas
g++ -o build/filterForDeepLearning_MiniPET DeepLearningDOI/filterForDeepLearning_MiniPET.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas
g++ -o build/filterForDOIscan DeepLearningDOI/filterForDOIscan.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas
g++ -o build/compton_volumes Cluster/compton_volumes.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

### MODULE CALIBRATION
echo "Compiling ModuleCalibration..."
cd build
echo "Running cmake command as: cmake $CMAKE_ARGS ../"
cmake $CMAKE_ARGS ../
make
cd ../
echo "Done."