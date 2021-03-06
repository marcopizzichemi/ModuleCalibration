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

def worker(name,files,mppc,histoMin,histoMax,histoBins,fitPercMin,fitPercMax):
    """thread worker function"""
    # value = 1
    filesMod = files + "*"
    cmd = ['ModuleCalibration','-c', name, filesMod]
    # cmd = "ModuleCalibration -c " + name + " " + prefix
    print ("Running calibration on MPPC %s..." %mppc )
    logName = 'log_pOutput_' + mppc + '.log'
    log = open(logName, 'w')
    subprocess.Popen(cmd,stdout = log,stderr=None).wait()

    # print (cmd)
    print ("MPPC %s calibration done" %mppc )
    print ("Running time analysis on MPPC %s..." %mppc )
    cmd = ['timeAnalysis','-i', files,'-o', 'time_' + mppc + '.root', '-c' , 'pOutput_' + mppc + '.root','--histoMin',str(histoMin) ,'--histoMax',str(histoMax) ,'--histoBins',str(histoBins),'--fitPercMin', str(fitPercMin),'--fitPercMax', str(fitPercMax) ]
    subprocess.Popen(cmd,stdout = log,stderr=log).wait()
    log.close()
    # print (cmd)
    print ("MPPC %s time analysis done" %mppc )

    return



def main(argv):

   #parsing args
   parser = argparse.ArgumentParser(description='Python script to start analysis in parallel')
   parser.add_argument('-c','--config' , help='Config file',required=True)
   parser.add_argument('-p','--path' , help='Full path of folder',required=True)
   parser.add_argument('-f','--files' , help='File prefix',required=True)
   parser.add_argument('-f','--crystal' , help='Crystal analyzed',required=True)
   parser.add_argument('-a','--histoMin' , help='Min in time histograms [s]',required=False)
   parser.add_argument('-b','--histoMax' , help='Max in time histograms [s]',required=False)
   parser.add_argument('-d','--histoBins' , help='Bins in time histograms [s]',required=False)
   parser.add_argument('-e','--fitPercMin' , help='Min in fitting range [perc. of mean]',required=False)
   parser.add_argument('-g','--fitPercMax' , help='Max in fitting range [perc. of mean]',required=False)

   args = parser.parse_args()


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

   #print values
   print ("Config file       = %s" % args.config )
   print ("Folder path       = %s" % args.path )
   print ("File prefix       = %s" % args.files )
   print ("HistoMin          = %s" % args.histoMin)
   print ("histoMax          = %s" % args.histoMax)
   print ("histoBins         = %s" % args.histoBins)
   print ("fitPercMin        = %s" % args.fitPercMin)
   print ("fitPercMax        = %s" % args.fitPercMax)
   print ("")

   # read base config file
   fInName = args.config
   with open(fInName, 'r') as myfile:
     original = myfile.read()

   # create output filename
   fOutName = fInName[0:len(fInName)-4] + "_" + args.crystal + ".cfg"

   mppclist = args.mppcs.split(",")
   processList = []
   # create config files and processes
   for i in mppclist:
       # print(i)
       fOutName = fInName[0:len(fInName)-4] + "_" + i + ".cfg"
       fOut = open(fOutName,"w")
       paraStr = "parallelMPPC = " + i
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
