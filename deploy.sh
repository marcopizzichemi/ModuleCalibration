#!/bin/bash

mkdir -p ../build
if [ -z "$1" ]
  then
    echo "No argument supplied, so compiling with machine compiler..."
    echo "Compiling support programs..."
    g++ -o ../build/timeResolution timeResolution.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas
    g++ -o ../build/extractConfiguration extractConfiguration.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore
    cd DeepLearningDOI
    g++ -o ../../build/calcDOIres calcDOIres.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas
    cd ..
    echo "Compiling ModuleCalibration..."
    cd ../build
    cmake ../code
    make
    cd ../code
    echo "Done."

else
  if [ $1 = "lxplus" ]; then
    set --
    echo "Sourcing compiler for lxplus..."
    source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
    ### old sources, left here for reference
    #source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
    #source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6/setup.sh
    #source /cvmfs/geant4.cern.ch/geant4/10.3/x86_64-slc6-gcc49-opt/bin/geant4.sh
    echo "Compiling support programs..."
    g++ -o ../build/timeResolution timeResolution.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas
    g++ -o ../build/extractConfiguration extractConfiguration.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore
    cd DeepLearningDOI
    g++ -o ../../build/calcDOIres calcDOIres.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

    echo "Compiling ModuleCalibration..."
    cd ../build
    cmake --verbose -DCMAKE_CXX_FLAGS=-std=c++17 -DCMAKE_CXX_COMPILER=`which g++` -DCMAKE_C_COMPILER=`which gcc` ../
    make
    cd ../code
    echo "Done."
  else
    echo "Invalid argument $1 - You can either provide no arg (and the script will use machine compiler), or lxplus, and the script will source the lxplus variables"
  fi
fi
