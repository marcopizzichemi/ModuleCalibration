#!/usr/bin/python3
# -*- coding: utf-8 -*-

import math
import os
import stat
import sys
import argparse
import subprocess
from subprocess import Popen, PIPE, STDOUT
import shutil
# import configparser
from io import StringIO
import csv
import pandas as pd

def main(argv):

    parser = argparse.ArgumentParser(description='Python script converting data')
    parser.add_argument('-n','--number'  , help='Chosen crystal number'                ,required=True)
    parser.add_argument('-c','--config'  , help='Config file for ModuleCalibration'    ,required=True)
    parser.add_argument('-p','--produce' , help='Produce data on individual z points'  ,action='store_true')
    parser.add_argument('-v','--view'    , help='View 2D cuts'                         ,action='store_true')
    parser.add_argument('-j','--join'    , help='Join all z point of a given crystal'  ,action='store_true')

    argv = parser.parse_args()

    crystalChoiceString = argv.number
    configFile          = argv.config
    produce             = argv.produce
    view                = argv.view
    join                = argv.join

    if produce == True:
        print ("Producing dataset..." )
    elif view == True:
        print ("Visualizing cuts..." )
    elif join == True:
        print ("Joining dataset..." )
    else:
        print ("You need to choose an action:")
        print ("--produce        Produce data on individual z points ")
        print ("--view           View 2D cuts                        ")
        print ("--join           Join all z point of a given crystal ")
        exit()
    # exit()

    print ("Crystal number       = %s" % crystalChoiceString         )
    print ("Config file          = %s" % configFile         )


    crystalChoice = int(crystalChoiceString)

    # make dataframe
    zPositions = pd.read_csv("/home/marco/cernbox/Universita/NewClearPEM/Programs/ModuleCalibration/code/parametersz_points_jinst_paper.txt/",sep=" ", header=None)
    zPositions.columns = ['i','j','z0','z1','z2','z3','z4','z5','z6','z7','z8','z9']
    crystalSeries = pd.DataFrame({'crystal': [63,62,61,60,59,58,57,56,55,54,53,52,51,50,49,48,47,46,45,44,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0]})

    columnSeries = pd.DataFrame({'column': ["y0","y0","y0","y0","y0","y0","y0","y0","y1","y1","y1","y1","y1","y1","y1","y1","y2","y2","y2","y2","y2","y2","y2","y2","y3","y3","y3","y3","y3","y3","y3","y3","y4","y4","y4","y4","y4","y4","y4","y4","y5","y5","y5","y5","y5","y5","y5","y5","y6","y6","y6","y6","y6","y6","y6","y6","y7","y7","y7","y7","y7","y7","y7","y7"]})

    mppcSeries = pd.DataFrame({'mppc': ["D4","D4","C4","C4","B4","B4","A4","A4","D4","D4","C4","C4","B4","B4","A4","A4","D3","D3","C3","C3","B3","B3","A3","A3","D3","D3","C3","C3","B3","B3","A3","A3","D2","D2","C2","C2","B2","B2","A2","A2","D2","D2","C2","C2","B2","B2","A2","A2","D1","D1","C1","C1","B1","B1","A1","A1","D1","D1","C1","C1","B1","B1","A1","A1"]})

    doiColumnOffsetSeries = pd.DataFrame({'doiColumnOffset': [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0] })

    dt = pd.concat([crystalSeries,columnSeries,mppcSeries,doiColumnOffsetSeries,zPositions], axis = 1)

    row = dt.loc[dt['crystal'] == crystalChoice]


    crystal            = str(row['crystal'].values[0])
    column             = str(row['column'].values[0])
    mppc               = str(row['mppc'].values[0])
    doiColumnOffset    = str(row['doiColumnOffset'].values[0])

    #set digiChannels
    if column == "y0":
        digiChannels="12,13,14,15"
    if column == "y1":
        digiChannels="12,13,14,15"
    if column == "y2":
        digiChannels="8,9,10,11"
    if column == "y3":
        digiChannels="8,9,10,11"
    if column == "y4":
        digiChannels="4,5,6,7"
    if column == "y5":
        digiChannels="4,5,6,7"
    if column == "y6":
        digiChannels="0,1,2,3"
    if column == "y7":
        digiChannels="0,1,2,3"


    for i in `seq 0 $((lenght-1))`
    do
       targetFolder="${yfolder}/z${i}/RootTTrees/"
       targetFile="new_config_${mppc}_${crystal}.cfg"
       cp $baseConfig ${targetFolder}/${targetFile}
       echo "dumpImages         = 1"                                       >> ${targetFolder}/${targetFile}
       echo "digiChannelsForDoi = ${digiChannels}"                         >> ${targetFolder}/${targetFile}
       echo "parallelMPPC       = ${mppc}"                                 >> ${targetFolder}/${targetFile}
       echo "output             = new_parseCalibration_${mppc}_${crystal}" >> ${targetFolder}/${targetFile}
       echo "taggingPosition    = ${zValues[$i]}"                          >> ${targetFolder}/${targetFile}
       echo "doiColumnOffset    = ${doiColumnOffset}"                      >> ${targetFolder}/${targetFile}
       echo ${zValues[$i]}
       cd ${targetFolder}
       time ModuleCalibration -c ${targetFile} TTree_*
       filterForDeepLearning --input TTree_ --calibration new_parseCalibration_${mppc}_${crystal}.root --crystal ${crystal} --output new_out${crystal}.root --photopeakCut --cutgCut
       cd -
    done

    cmd = ['time','./produceDataset.sh',configFile, column,mppc,crystal,doiColumnOffset,str(row['z0'].values[0]), str(row['z1'].values[0]), str(row['z2'].values[0]), str(row['z3'].values[0]), str(row['z4'].values[0]), str(row['z5'].values[0]), str(row['z6'].values[0]), str(row['z7'].values[0]), str(row['z8'].values[0]), str(row['z9'].values[0])]

    print (cmd)
    cmdString = ' '.join(cmd)
    print (cmdString)

    logName = 'log_' + crystal + '.log'
    log = open(logName, 'w')
    subprocess.Popen(cmd,stdout = log,stderr=None).wait()






if __name__ == "__main__":
    main(sys.argv[1:])
