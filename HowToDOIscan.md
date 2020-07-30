# **DOI SCAN HOWTO**

## 1. Analyze top dataset with ModuleCalibration


#### Equalize the MPPCs gain

For simplicity, analyze the entire ```TOP``` acq dataset together. Take the example config file in

```
./config/baseConfig.cfg
```

Be sure to check that you are selecting the correct channels for Array and Tag. This is set by the keys

```
digitizer = 16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31
timing    = 16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31
taggingCrystalChannel = 52
taggingCrystalTimingChannel = 52
```

Also set the other parameters properly (saturation etc). The config file should be ready to analyze just quickly the matrix, and dump the position of 511 kev photopeaks. Run

```
ModuleCalibration -c baseConfig.cfg TTree_0.root TTree_1.root TTree_2.root ...
```

where the ```TTree_N.root``` files are the individual top scan files. If you have performed for example the top scan in 16 positions (one on top of each crystal in a 4x4 array), then you can just put 1 file for each of these scan, the statistics should be enough. Otherwise, put N files for each position (try to have the same statistics for each position).

ModuleCalibration will run, write a ```calibration.root``` output file. Take a look at the ```calibration.root``` to see if all the relevant plots make sense. In particular, check

```
Flood Histogram 2D     --> check that all crystals are there
Flood Histogram 3D     --> check that there is nothing strange
TaggingCrystalSpectrum --> check that the photopeak is properly found and fitted
BigSingleSpectra       --> check that all 511 KeV photopeaks are found and fitted
```

You can also check the that individual 3D cuts in the MPPC folders are ok.

ModuleCalibration has also dump a file ```singlePeaksFile.txt``` with the centroid of the positions of 511 KeV photopeaks for single MPPC signals. Now you can put this information in the baseConfig.cfg, in the key

```
gain = 29310,30760,38320,42840,31380,37640,29910,31600,36170,33210,40060,32830,29020,29980,31590,32630
```

where you should list the data written in the ```singlePeaksFile.txt```, but we careful to put them in the same order of the MPPC labels written in the key ```mppc```.

#### Produce calibration file

You can now re-run ModuleCalibration with the same command as before. If needed, you can also run already the timing correction analysis on the TOP dataset, by setting

```
timingCorrection = 1
timingCorrectionForPolished = 1
```

Running the **ModuleCalibration** command exactly as above will overwrite the output file

```
calibration.root
```

(actual file name anyway depends on the value of key ```output``` in config file, so if you want you can change that and avoid overwrite) that will contain the information about the crystal limits (and the timing correction plots if you've run the timing correction analysis).


## 2. Perform lateral scan analysis

Move to one line folder of the lateral scan. Copy the script

```
./scripts/startDOIscanAnalysis.sh
```
in this folder. Run it like this

```
./startDOIscanAnalysis.sh [calibrationFile] [zmin] [zmax] [points] [cryN1 cryN2 ...]

where

[calibrationFile] = full path to calibration.root file you just produced
[zmin]            = min z for the DOI analysis
[zmin]            = max z for the DOI analysis
[points]          = number of vertical points taken
[cryN1 cryN2]     = list of crystal numbers in this line
```

So for example for folder 4 (assuming it's the rightmost line of the array), you would give the command

```
./startDOIscanAnalysis.sh /path/to/calibration.root 4 11 15 12 13 14 15
```

assuming you want to analyze only from z=4 to z=11 and you took 15 z points in the scan, and the crystals in this line are number 12 13 14 and 15.

**Beware**: you need to adapt the startDOIscanAnalysis to your folders, in particular

```
Run_scan_doi_ACQ_0.${i}_Angle1_0.0_Angle2_0.0_V1_59_V2_62_t_3600/RootTTrees/

needs to reflect your subfolders
```

**Beware/2**: adapt the cut in tag photopeak:

```
--tagCh 35 --manTagMin 40580.2 --manTagMax 43541.8

are needed only if the TOP and LATERAL acqs were performed with different tag channels, otherwise you can remove all 3 flags. If they are needed, but the proper values for your case.
```

**Beware/3**: all this works if the lateral scan is in step of 1 mm, otherwise we need a more complicated script.



This will produce

```
plotsN.txt
plotsN.root
```

The ROOT file contains several plots:

```

```

The .txt file has the DOI FWHM results. You are interested in the column ```FIT_calib```. Repeat the operation for all lines.


## 3. Perform timing and energy resolution analysis

If you have run **ModuleCalibration** with timing correction ON, you can now complete the analysis on the TOP dataset. It's more convenient to run on a dataset that includes all crystals, so that it can be done in one go. If you move to the main folder (the one that contains both lateral and top folders), you can merge datasets of different crystals doing

```
hadd -f total.root top_4_sectors/Run_top_scan_ACQ_0.0_Angle1_0.0_Angle2_0.0_V1_59_V2_62_t_43200/RootTTrees/TTree_0_1.root top_4_sectors/Run_top_scan_ACQ_1.0_Angle1_0.0_Angle2_0.0_V1_59_V2_62_t_21600/RootTTrees/TTree_0_1.root top_4_sectors/Run_top_scan_ACQ_2.0_Angle1_0.0_Angle2_0.0_V1_59_V2_62_t_21600/RootTTrees/TTree_0_1.root top_4_sectors/Run_top_scan_ACQ_3.0_Angle1_0.0_Angle2_0.0_V1_59_V2_62_t_21600/RootTTrees/TTree_0_1.root
```

This will produce a complete dataset file in total.root. The command above assumes that your TOP folder is **top_4_sectors**, and that it consisted of 4 acquisitions, but of course you need to adjust it to reflect your specific case.

Now you can analyze the entire dataset with

```
timeResolution --calibration analysis.root --folder ./ --prefix total --output timeTot.root --histoMin -2e-9 --histoMax 2e-9 --histoBins 200 --func 1 --tagFwhm 90e-12
```
