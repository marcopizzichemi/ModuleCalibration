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

def worker(name,files,element,histoMin,histoMax,histoBins,fitPercMin,fitPercMax):
    """thread worker function"""
    # value = 1
    filesMod = files + "*"
    cmd = ['ModuleCalibration','-c', name, filesMod]
    # cmd = "ModuleCalibration -c " + name + " " + prefix
    print ("Running calibration on Element %s..." %element )
    logName = 'log_pOutput_' + element + '.log'
    log = open(logName, 'w')
    subprocess.Popen(cmd,stdout = log,stderr=None).wait()

    # print (cmd)
    print ("Element %s calibration done" %element )
    print ("Running time analysis on Element %s..." %element )
    cmd = ['timeAnalysis','-i', files,'-o', 'time_' + element + '.root', '-c' , 'pOutput_' + element + '.root','--histoMin',str(histoMin) ,'--histoMax',str(histoMax) ,'--histoBins',str(histoBins),'--fitPercMin', str(fitPercMin),'--fitPercMax', str(fitPercMax) ]
    subprocess.Popen(cmd,stdout = log,stderr=log).wait()
    log.close()
    # print (cmd)
    print ("Element %s time analysis done" %element )

    return



def main(argv):

   #parsing args
   parser = argparse.ArgumentParser(description='Python script to start analysis in parallel')
   parser.add_argument('-c','--config'     , help='Config file'               ,required=True )
   parser.add_argument('-m','--mppcs'      , help='MPPC analyzed'             ,required=False)
   parser.add_argument('-r','--crystals'   , help='Crystals analyzed'         ,required=False)
   parser.add_argument('-f','--files'      , help='File prefix'               ,required=True )
   # parser.add_argument('-t','--threads'    , help='Number of parallel threads',required=False)
   parser.add_argument('-a','--histoMin'   , help='Min of CTR histograms [s] - default = -15e-9',required=False)
   parser.add_argument('-b','--histoMax'   , help='Max of CTR histograms [s] - default = 15e-9',required=False)
   parser.add_argument('-d','--histoBins'  , help='Number of bins in CTR histograms - default = 300',required=False)
   parser.add_argument('-e','--fitPercMin' , help='Lower bound of CTR fit, in numb of sigmas - default = 6.0',required=False)
   parser.add_argument('-g','--fitPercMax' , help='Upper bound of CTR fit, in numb of sigmas - default = 5.0',required=False)

   args = parser.parse_args()

   # threads = 0
   # if args.threads == None:
   #     args.threads = 8
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
   print ("Config file        = %s" % args.config )
   print ("Elements analyzed  = %s" % elements )
   print ("File prefix        = %s" % args.files )
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
       outputFilePrefix = "parallelOutput = pOutput_" + i
       fOut.write("### Parallelized analysis ###\n")
       fOut.write(paraStr)
       fOut.write("\n")
       fOut.write(outputFilePrefix)
       fOut.write("\n")
       fOut.write("\n")
       fOut.write(original)
       fOut.close()
       proc = multiprocessing.Process(target=worker, args=(fOutName,args.files,i,args.histoMin,args.histoMax,args.histoBins,args.fitPercMin,args.fitPercMax))
       processList.append(proc)

   #start processes
   for i in range(len(processList)):
       processList[i].start()




if __name__ == "__main__":
   main(sys.argv[1:])
