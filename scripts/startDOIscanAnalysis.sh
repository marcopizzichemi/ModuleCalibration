#!/bin/bash

calibration=$1
zmin=$2
zmax=$3
points=$4

let "points=points-1"

shift 4

for j in "$@"
do
  VAR=""
  for i in `seq 0 ${points}`
  do
    # echo $j $i
    VAR+="output_cry${j}_z${i}.root "
    filterForDOIscan -f Run_scan_doi_ACQ_0.${i}_Angle1_0.0_Angle2_0.0_V1_59_V2_62_t_3600/RootTTrees/ -p TTree_0_0.root -c ${calibration} -o output_cry${j}_z${i}.root -n $j --tagCh 35 --manTagMin 40580.2 --manTagMax 43541.8 --tagCut --cutgCut --photopeakCut --tagPos $i
  done
  # rm -rf crystal${i}.root
  # echo $VAR
  hadd -f crystal${j}.root $VAR
done

for j in "$@"
do
  # echo #j
  calcDOIres -i crystal${j}.root -o plots$j -c ${calibration} -n $j --zmin ${zmin} --zmax ${zmax}
done
