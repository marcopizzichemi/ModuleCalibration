#!/usr/bin/python3
# -*- coding: utf-8 -*-

import math
import os
import stat
import sys
import argparse
import subprocess
from subprocess import Popen, PIPE, STDOUT
import threading
import time
import multiprocessing

def worker(buildPath,name,filesCalib,filesTime,element,histoMin,histoMax,histoBins,fitPercMin,fitPercMax,prefix_name,func,fitCorrection,excludeChannels,excludeChannelsList,inputCalibFolder,inputTimeFolder,tagFwhm):
    """thread worker function"""
    # value = 1
    executable = buildPath + '/ModuleCalibration'
    cmd = [executable,'-c', name]

    # filesMod = ""
    for i in filesCalib.split():
        for file in os.listdir(inputCalibFolder):
            if file.startswith(i):
                cmd.append(inputCalibFolder+ '/' +file)
    # print (filesMod)



    # cmd = ['ModuleCalibration','-c', name, filesMod]
    # cmd = "ModuleCalibration -c " + name + " " + prefix
    print ("Running calibration on Element %s..." %element )
    cmdString = ' '.join(cmd)
    print (cmdString)
    logName = 'log_' + prefix_name + element + '.log'
    log = open(logName, 'w')
    subprocess.Popen(cmd,stdout = log,stderr=None).wait()

    # print (cmd)
    print ("Element %s calibration done" %element )
    print ("Running time analysis on Element %s..." %element )
    #timeResolution --input-folder ./ --input TTree_ --calibration ../calibrationDataset/calibration6.root -o output6.root --histoMin -0.5e-9 --histoMax 1.5e-9 --histoBins 500 --func 0 --fitPercMin 6.0 --fitPercMax 5.0
    executable = buildPath + '/timeResolution'
    cmd = [executable,'-f', inputTimeFolder,'-p', filesTime,'-o', 'time_' + prefix_name + element + '.root', '-c' , prefix_name + element + '.root','--histoMin',str(histoMin) ,'--histoMax',str(histoMax) ,'--histoBins',str(histoBins),'--func',str(func),'--fitPercMin', str(fitPercMin),'--fitPercMax', str(fitPercMax),'--tagFwhm', str(tagFwhm) ]
    if fitCorrection == 1:
        cmd.append('--fitCorrection')
    if excludeChannels == 1:
        cmd.append('--exclude-channels')
        cmd.append(excludeChannelsList)

    cmdString = ""
    cmdString = ' '.join(cmd)
    print (cmdString)
    subprocess.Popen(cmd,stdout = log,stderr=log).wait()
    log.close()
    # print (cmd)
    print ("Element %s time analysis done" %element )

    return



def main(argv):

   #parsing args
   parser = argparse.ArgumentParser(description='Python script to start analysis in parallel')
   parser.add_argument('-c','--config'     , help='Config file'                  ,required=True )
   parser.add_argument('-m','--mppcs'      , help='MPPC analyzed'                ,required=False)
   parser.add_argument('-r','--crystals'   , help='Crystals analyzed'            ,required=False)
   parser.add_argument('-p','--inputCalibFolder'  , help='Folder of input files for Calibration'    ,required=True )
   parser.add_argument('-f','--filesCalib'        , help='File prefix for calibration'  ,required=True )
   parser.add_argument('-w','--inputTimeFolder'   , help='Folder of input files for time resolution'    ,required=True )
   parser.add_argument('-l','--filesTime'         , help='File prefix for time resolution'    ,required=True )
   parser.add_argument('-a','--histoMin'   , help='Min of CTR histograms [s] - default = -15e-9',required=False)
   parser.add_argument('-b','--histoMax'   , help='Max of CTR histograms [s] - default = 15e-9',required=False)
   parser.add_argument('-d','--histoBins'  , help='Number of bins in CTR histograms - default = 300',required=False)
   parser.add_argument('-e','--fitPercMin' , help='Lower bound of CTR fit, in numb of sigmas - default = 6.0',required=False)
   parser.add_argument('-g','--fitPercMax' , help='Upper bound of CTR fit, in numb of sigmas - default = 5.0',required=False)
   parser.add_argument('-o','--output'     , help='Prefix of output file - default = pOutput_',required=False)
   parser.add_argument('-y','--func'       , help='Function for time fit. 0 = crystalball, 1 = gauss+exp   - default = 0',required=False)
   parser.add_argument('-q','--tagFwhm'    , help='FWHM of reference crystal - default = 88.0e-12',required=False)
   parser.add_argument('-k','--fitCorrection'       , help='0 = no fit, 1 = use fit   - default = 0',required=False)
   parser.add_argument('-j','--excludeChannels'       , help='csv list of channels to exclude from time corrections   - default = ""',required=False)
   parser.add_argument('-x','--build'     , help='Path to build folder',required=True )
   args = parser.parse_args()

   # threads = 0
   # if args.threads == None:
   #     args.threads = 8
   prefix_name = "pOutput_"
   inputCalibFolder = "./"
   inputTimeFolder = "./"
   func = 0
   fitCorrection = 0
   excludeChannels = 0

   if args.inputCalibFolder != None:
       inputCalibFolder = args.inputCalibFolder
   if args.inputTimeFolder != None:
       inputTimeFolder = args.inputTimeFolder

   if args.fitCorrection != None:
       fitCorrection = int(args.fitCorrection)
   if args.func != None:
       func = args.func
   if args.output != None:
       prefix_name = args.output
   if args.histoMin == None:
       args.histoMin = -15e-9
   if args.histoMax == None:
       args.histoMax = 15e-9
   if args.histoBins == None:
       args.histoBins = 300
   if args.fitPercMin == None:
       args.fitPercMin = 6.0
   if args.fitPercMax == None:
       args.fitPercMax = 5.0
   if args.excludeChannels != None:
       excludeChannels = 1
   if args.tagFwhm == None:
       args.tagFwhm = 88.0e-12

   elements = ""
   BaseParaStr = ""

   if (args.mppcs == None) and (args.crystals == None):
       print ("ERROR: you need to provide one option between --mppcs or --crystals !!! Aborting ")
       sys.exit()
   if (args.mppcs != None) and (args.crystals != None):
       print ("ERROR: you cannot provide both --mppcs and --crystals options! Aborting ")
       sys.exit()
   if (args.mppcs == None) and (args.crystals != None):
       elements = args.crystals
       elementslist =  args.crystals.split(",")
       BaseParaStr = "parallelCrystal = "

   if (args.mppcs != None) and (args.crystals == None):
       elements = args.mppcs
       elementslist =  args.mppcs.split(",")
       BaseParaStr = "parallelMPPC = "
   #print values
   print ("Config file                 = %s" % args.config )
   print ("Element(s) analyzed         = %s" % elements )
   print ("Calibration files folder    = %s" % inputCalibFolder )
   print ("Calibration files prefix    = %s" % args.filesCalib )
   print ("Time analysis files folder  = %s" % inputTimeFolder )
   print ("Time analysis files prefix  = %s" % args.filesTime )
   print ("Build folder                = %s" % args.build )
   print ("")
   # print ("Number of threads = %s" % args.threads )

   fInName = args.config
   with open(fInName, 'r') as myfile:
     original = myfile.read()

   print (elementslist)



   processList = []
   # create config files and processes
   for i in elementslist:
       # print(i)
       fOutName = fInName[0:len(fInName)-4] + "_" + i + ".cfg"
       fOut = open(fOutName,"w")
       paraStr = BaseParaStr + i
       outputFilePrefix = "parallelOutput = " + prefix_name + i
       fOut.write("### Parallelized analysis ###\n")
       fOut.write(paraStr)
       fOut.write("\n")
       fOut.write(outputFilePrefix)
       fOut.write("\n")
       fOut.write("\n")
       if(excludeChannels):
           fOut.write("### Channels excluded from polished time correction, modification by timeAnalysis ###\n")
           keyString = "excludeChannels = " + args.excludeChannels
           fOut.write(keyString)
           fOut.write("\n")
           fOut.write("\n")
       fOut.write(original)
       fOut.close()
       proc = multiprocessing.Process(target=worker, args=(args.build,fOutName,args.filesCalib,args.filesTime,i,args.histoMin,args.histoMax,args.histoBins,args.fitPercMin,args.fitPercMax,prefix_name,func,fitCorrection,excludeChannels,args.excludeChannels,inputCalibFolder,inputTimeFolder,args.tagFwhm))
       processList.append(proc)

   #start processes
   for i in range(len(processList)):
       processList[i].start()




if __name__ == "__main__":
   main(sys.argv[1:])
