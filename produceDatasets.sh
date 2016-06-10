#!/bin/bash

# run this from the DOI_scan folder

runList='y3'
zlist='z0 z1 z2 z3 z4 z5 z6 z7 z8 z9'

for i in $runList ; do for j in $zlist; do cp $i/$j/Run*/doiData.txt $i/$j/doiData.txt ; done; done

for i in $runList 
do
  cd $i
  doifit `echo $zlist`
  cd ..
done

for i in $runList ; do cp $i/calibration_params.txt ./all/${i}calibration_params.txt ; cp $i/doi.txt ./all/${i}doi.txt; done

#delete the old global files, create the new ones
cd all
rm calibration_params.txt
rm doi.txt

for i in $runList ; do  cat ${i}calibration_params.txt >> calibration_params.txt; done
for i in $runList ; do  cat ${i}doi.txt >> doi.txt; done

cd ..
