// compile with
// g++ -o ../build/extractConfiguration extractConfiguration.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

// small program to extrac config file and pc info from the results of ModuleCalibration.
// run extractConfiguration without arguments for usage info

#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"
#include "TEllipse.h"
#include "TFormula.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "TCutG.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TGraphDelaunay.h"
#include "TVector.h"

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>

void usage()
{
  std::cout << "\t\t" << "[-i | --input]   <moduleCalibration file> " << std::endl
            << "\t\t" << "[--all|--hostname|--pwd|--files|--config]   - choice of info to be printed. At least one has to be chosen!"   << std::endl
            << "\t\t" << "[-o | --output]  <output file>             - OPTIONAL: if no output file given, results will be print on screen> ] " << std::endl
            << "\t\t" << std::endl;
}

int main (int argc, char** argv)
{
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }


  std::string inputFileName;
  std::string outputFileName;
  bool outToTextFile   = false;
  bool includeAll      = false;
  bool includeConfig   = false;
  bool includePWD      = false;
  bool includeHostname = false;
  bool includeFiles    = false;
  bool choiceMade      = false;

  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "all", no_argument, 0, 0 },
      { "hostname", no_argument, 0, 0 },
      { "pwd", no_argument, 0, 0 },
      { "files", no_argument, 0, 0 },
      { "config", no_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			inputFileName = (char *)optarg;
    }
		else if (c == 'o'){
      outputFileName = (char *)optarg;
      outToTextFile = true;
    }
		else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      outputFileName = (char *)optarg;
      outToTextFile = true;
    }
    else if (c == 0 && optionIndex == 2){
      includeAll = true;
      choiceMade = true;
    }
    else if (c == 0 && optionIndex == 3){
      includeHostname = true;
      choiceMade = true;
    }
    else if (c == 0 && optionIndex == 4){
      includePWD = true;
      choiceMade = true;
    }
    else if (c == 0 && optionIndex == 5){
      includeFiles = true;
      choiceMade = true;
    }
    else if (c == 0 && optionIndex == 6){
      includeConfig = true;
      choiceMade = true;
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}

  if(!choiceMade)
  {
    std::cout	<< "Usage: " << argv[0] << std::endl;
		usage();
		return 1;
  }
  TFile *file = TFile::Open(inputFileName.c_str());
  file->cd("Configuration");
  std::ofstream outFile;
  if(outToTextFile)
  {
    outFile.open(outputFileName.c_str(),std::ifstream::out);
    if(includeAll | includeHostname)
    {
      outFile << "Hostname                 = ";
      outFile << ((TNamed*) gDirectory->Get("Hostname"))->GetTitle() << std::endl;
    }
    if(includeAll | includePWD)
    {
      outFile << "Working Directory        = ";
      outFile << ((TNamed*) gDirectory->Get("PWD"))->GetTitle() << std::endl;
    }
    if(includeAll | includeFiles)
    {
      outFile << "Input File(s):" << std::endl;
      outFile << ((TNamed*) gDirectory->Get("InputFiles"))->GetTitle() << std::endl;
    }
    if(includeAll | includeConfig)
    {
      outFile << std::endl;
      outFile << std::endl;
      outFile << "------------------------------------------- Configuration File -------------------------------------------------- " << std::endl;
      outFile << ((TNamed*) gDirectory->Get("ConfigFile"))->GetTitle() << std::endl;
      outFile << "---------------------------------------  END OF Configuration File -------------------------------------------------- " << std::endl;
    }
    outFile.close();
  }
  else
  {
    std::cout << std::endl;
    if(includeAll | includeHostname)
    {
      std::cout << "Hostname                 = ";
      std::cout << ((TNamed*) gDirectory->Get("Hostname"))->GetTitle() << std::endl;
    }
    if(includeAll | includePWD)
    {
      std::cout << "Working Directory        = ";
      std::cout << ((TNamed*) gDirectory->Get("PWD"))->GetTitle() << std::endl;
    }
    if(includeAll | includeFiles)
    {
      std::cout << "Input File(s):" << std::endl;
      std::cout << ((TNamed*) gDirectory->Get("InputFiles"))->GetTitle() << std::endl;
    }
    if(includeAll | includeConfig)
    {
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "------------------------------------------- Configuration File -------------------------------------------------- " << std::endl;
      std::cout << ((TNamed*) gDirectory->Get("ConfigFile"))->GetTitle() << std::endl;
      std::cout << "---------------------------------------  END OF Configuration File -------------------------------------------------- " << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
  return 0;
}
