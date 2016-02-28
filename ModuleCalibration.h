#include "TROOT.h"
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
#include "TFile.h"
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
#include "TCutG.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TPaveText.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdlib.h> 
#include <stdio.h> 
#include <unistd.h>
#include <cmath> 
#include <assert.h>
#include <cstddef>

#include "ConfigFile.h"
#include "InputFile.h"
#include "Element.h"
#include "Crystal.h"
#include "Module.h"
#include "Mppc.h"

#define ENERGY_RESOLUTION 0.12
#define ENERGY_RESOLUTION_SATURATION_CORRECTION 0.25

/*
  int ncrystalsx                = config.read<int>("ncrystalsx",2);                 // number of crystals in x direction per mppc - default to 2 if the key is not found in the config file
  int ncrystalsy                = config.read<int>("ncrystalsy",2);                 // number of crystals in y direction per mppc - default to 2 if the key is not found in the config file
  int nmppcx                    = config.read<int>("nmppcx",2);                     // number of mppc in x direction per mppc - default to 2 if the key is not found in the config file
  int nmppcy                    = config.read<int>("nmppcy",2);                     // number of mppc in y direction per mppc - default to 2 if the key is not found in the config file
  int nmodulex                  = config.read<int>("nmodulex",1);                   // number of modules in x direction per mppc - default to 1 if the key is not found in the config file
  int nmoduley                  = config.read<int>("nmoduley",1);                   // number of modules in y direction per mppc - default to 1 if the key is not found in the config file
  int histo1Dmax                = config.read<int>("histo1Dmax");                   // max of the 1D charge histograms (in ADC channels)
  int histo1Dbins               = config.read<int>("histo1Dbins");                  // number of bins of the 1D charge histograms
  int histo2DchannelBin         = config.read<int>("histo2DchannelBin");            // number of bins of the 2D flood histograms, for single channels
  int histo2DglobalBins         = config.read<int>("histo2DglobalBins");            // number of bins of the 2D flood histograms, for entire module
  int histo3DchannelBin         = config.read<int>("histo3DchannelBin");            // number of bins of the 3D flood histograms, for single channels
  int histo3DglobalBins         = config.read<int>("histo3DglobalBins");            // number of bins of the 3D flood histograms, for entire module
  int taggingPeakMin            = config.read<int>("taggingPeakMin",8000);          // min range of tagging crystal photopeak, in ADC channels - to help TSpectrum
  int taggingPeakMax            = config.read<int>("taggingPeakMax",12000);         // max range of tagging crystal photopeak, in ADC channels - to help TSpectrum
  bool saveAnalysisTree         = config.read<bool>("saveAnalysisTree");            // choice to save or not the analysis TTree, in a file temp.root
  float taggingPosition         = config.read<float>("taggingPosition");            // position of the tagging bench in mm 
  bool usingTaggingBench        = config.read<bool>("usingTaggingBench");           // true if the input is using tagging bench, false if not
  int taggingCrystalChannel     = config.read<int>("taggingCrystalChannel");        // input channel where the tagging crystal information is stored
  bool correctingSaturation     = config.read<bool>("correctingSaturation");        // true if saturation correction is applied, false if it's not
  bool correctingForDOI         = config.read<bool>("correctingForDOI",0);            // true if the energy correction using DOI info is computed    
  float energyResolution        = config.read<float>("expectedEnergyResolution",0); // energy resolution input by the user, if any, otherwise 0
  bool usingRealSimData         = config.read<bool>("usingRealSimData",0);
  // --- paramenters for roto-translations to separate the nXn peaks
  // lateral, not corners
  double base_lateralQ1         = config.read<double>("lateralQ1",0.905);           // right and left
  double base_lateralQ2         = config.read<double>("lateralQ2",1.1);             // top and bottom
  double base_lateralDeltaU     = config.read<double>("lateralDeltaU",1);           // used for right and left
  double base_lateralDeltaV     = config.read<double>("lateralDeltaV",1);           // used for top and bottom
  double base_lateralRescaleRL  = config.read<double>("lateralRescaleRL",1.5);      // used for right and left 
  double base_lateralRescaleTB  = config.read<double>("lateralRescaleTB",2);        // used for top and bottom
  // corners                                                                    
  double base_cornerQ1          = config.read<double>("cornerQ1",0.675);            // rotation around Z
  double base_cornerQ2          = config.read<double>("cornerQ2",1.41);             // rotation around X
  double base_cornerDeltaU      = config.read<double>("cornerDeltaU",3);            // translations
  double base_cornerDeltaV      = config.read<double>("cornerDeltaV",2.1);          // translations
  double base_cornerRescale     = config.read<double>("cornerRescale",4);           // rescale factor
  bool   onlyuserinput          = config.read<double>("onlyuserinput",0);           // ignore 2d automatic fitting  
  // set output file name                                                   
  std::string outputFileName = config.read<std::string>("output");
*/