//------------------------------------------------------------//
//                                                            //
//  PROGRAM FOR matrix calibration USING DT5740               //
//                                                            //
//------------------------------------------------------------//

// compile with cmake and make

// this program is meant to read a "short" acquisition data and find the limits that will be used by the reconstruction algorithm to
// do the images. Limits will be searched on u,v,w and e, i.e

// u = weighted average in the x direction
// v = weighted average in the y direction
// w = ratio highest charge / total charge collected -> converted to z
// e = energy of the incident particle, to select photopeaks

// CAREFUL: this program is meant to run on a single module, so if the acquisition is done in parallel on the two modules, it is
// supposed to be run twice, with different configuration files that will select the proper input channels, and with, as an output,
// different text files

// CAREFUL: at the moment this program is used to characterize different modules. So instead of really looking for the limits, it takes
// limits (for the crystal positions in u,v) and computes some relevat paramenters, especially the peak position, energy resolutions and
// width (FWHM) of the w histograms for each crystal.
// The user can input the crystal limits in the config file, only crystals that are specified in the config file will be analized (see newBoardConfig.cfg)

// PROGRAM STRUCTURE
// This file performes the relevant analysis, while the module structure is defined by the Element, Module, Mppc and Crystal classes. Module, Mppc and Crystal
// inherit from Element, but also have some specific method
// Config files are read in various parts of the program, via the ConfigFile class, to input the user parameters
// So ModuleCalibration does:
// - some initial check
// - read the input, via the InputFile class (see those files). This will read the Tchain of input root files and creat a specific TTree that will be used for the analysis
// - read config parameters
// - Create the modules, mppcs, crystals, according to the user requests
// - Run on all the elements and produce the desided plots. Mainly 2D and 3D histograms, crystal and mppc spectra, w histograms and so on (see in the code below)
// - Creates some Canvases for nice plot storing
// - Creates an output file with a directory structure and saves the plots
// - Optionally, saves the analysis TTree in a root file






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
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "TCutG.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TGraphDelaunay.h"

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
// #include <thread>
#include <getopt.h>

#include "ConfigFile.h"
#include "InputFile.h"
#include "Element.h"
#include "Crystal.h"
#include "Module.h"
#include "Mppc.h"
#include "ModuleCalibration.h"
// #include "Detector.h"

// #include <omp.h>
//test

#define ENERGY_RESOLUTION 0.07
#define ENERGY_RESOLUTION_SATURATION_CORRECTION 0.12


void extractFromCrystalBall(TH1F* histo,double fitPercMin,double fitPercMax, int divs, double* res)
{
  //first, dummy gaussian fit
  TCanvas *cTemp  = new TCanvas("temp","temp");
  TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  histo->Fit(gaussDummy,"QN");

  double f1min = histo->GetXaxis()->GetXmin();
  double f1max = histo->GetXaxis()->GetXmax();
  // std::cout << f1min << " " << f1max << std::endl;
  TF1* f1  = new TF1("f1","crystalball");
  f1->SetLineColor(kBlack);
  f1->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  double fitMin = gaussDummy->GetParameter(1) - fitPercMin*(gaussDummy->GetParameter(2));
  double fitMax = gaussDummy->GetParameter(1) + fitPercMax*(gaussDummy->GetParameter(2));
  if(fitMin < f1min)
  {
    fitMin = f1min;
  }
  if(fitMax > f1max)
  {
    fitMax = f1max;
  }
  histo->Fit(f1,"Q","",fitMin,fitMax);
  double min,max,min10,max10;
  // int divs = 3000;
  double step = (f1max-f1min)/divs;
  // is [0] the max of the function???
  for(int i = 0 ; i < divs ; i++)
  {
    if( (f1->Eval(f1min + i*step) < f1->GetParameter(0)/2.0) && (f1->Eval(f1min + (i+1)*step) > f1->GetParameter(0)/2.0) )
    {
      min = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) > f1->GetParameter(0)/2.0) && (f1->Eval(f1min + (i+1)*step) < f1->GetParameter(0)/2.0) )
    {
      max = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) < f1->GetParameter(0)/10.0) && (f1->Eval(f1min + (i+1)*step) > f1->GetParameter(0)/10.0) )
    {
      min10 = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) > f1->GetParameter(0)/10.0) && (f1->Eval(f1min + (i+1)*step) < f1->GetParameter(0)/10.0) )
    {
      max10 = f1min + (i+0.5)*step;
    }
  }
  res[0] = f1->GetParameter(1);  // res[0] is mean
  res[1] = max-min;              // res[1] is FWHM
  delete cTemp;
}


// Double_t thetaFunction(Double_t *x, Double_t *par);

int main (int argc, char** argv)
{

  // gStyle->SetOptStat(0);

  //----------------------------------------------------------//
  //  Check input                                             //
  //----------------------------------------------------------//
  if(argc<2)
  {
    std::cout << " Usage: " 										<< std::endl;
    std::cout << " ModuleCalibration [-c config file ] <input-files> " 					<< std::endl;
    std::cout << "   note: -c option is optional, but if you use it, it has to be the first argument"	<< std::endl;
    std::cout << "         without it, a configFile.cfg will be assumed" 				<< std::endl;
    return 0;
  }




  //----------------------------------------------------------//
  //  Say hello to the user                                   //
  //----------------------------------------------------------//
  std::cout<<"\n"<<std::endl;
  std::cout<<"###########################################################"<<std::endl;
  std::cout<<"#                                                         #"<<std::endl;
  std::cout<<"#                   Module calibration                    #"<<std::endl;
  std::cout<<"#                                                         #"<<std::endl;
  std::cout<<"###########################################################"<<std::endl;
  std::cout<<"\n\n"<<std::endl;
  std::cout<<"=====>   C O N F I G U R A T I O N   <====\n"<<std::endl;


  //----------------------------------------------------------//
  //  Enable Implicit MultiThreading of ROOT                  //
  //----------------------------------------------------------//
  // unsigned int nthreads = std::thread::hardware_concurrency();
  // std::cout << nthreads << " concurrent threads are supported " << std::endl;
  // ROOT::EnableImplicitMT(nthreads-1);
  // std::cout << nthreads-1 << " concurrent threads will be used " << std::endl;

  //----------------------------------------------------------//
  //  Import input and parse the config file                  //
  //----------------------------------------------------------//
  // Set a default config file name
  std::string ConfigFileName = "config.cfg";
  // and assume as default that there is no config file name from command line
  // then check
  if(std::string(argv[1]) == std::string("-c")) // first argument is -c, then the config file name is passed by command line
  {
    ConfigFileName = argv[2];
    std::cout << "Configuration file: '" << argv[2] << "'"<< std::endl;
  }
  else // the config file was indeed the default one
  {
    std::cout << "Configuration file set to default: config.cfg "<< std::endl;
  }

  //PUT THE ENTIRE CONFIG FILE IN A STRINGSTREAM, TO SAVE IT LATER IN THE OUTPUT FILE
  std::stringstream streamConfigFile;
  std::string strConfig;
  std::ifstream inConfig;
  inConfig.open (ConfigFileName.c_str(), std::ifstream::in);
  while (std::getline(inConfig, strConfig))
  {
    streamConfigFile << strConfig << std::endl;
  }
  inConfig.close();
  ConfigFile config(ConfigFileName); // create a ConfigFile object
  InputFile input(config); // read the input chain of root files, passing the inputs and the config object


  /////////////////////////////////////////////
  // READ CONFIG FILE                        //
  ////////////////////////////////////////////7
  // std::string chainName         = config.read<std::string>("chainName","adc");
  // int digitizerTotalCh                 = config.read<int>("digitizerTotalCh");
  // int digitizerType               = config.read<int>("digitizerType",0);
  std::string fname                       = config.read<std::string>("chainName","adc");
  int ncrystalsx                = config.read<int>("ncrystalsx",2);                 // number of crystals in x direction per mppc - default to 2 if the key is not found in the config file
  int ncrystalsy                = config.read<int>("ncrystalsy",2);                 // number of crystals in y direction per mppc - default to 2 if the key is not found in the config file
  int nmppcx                    = config.read<int>("nmppcx",2);                     // number of mppc in x direction per mppc - default to 2 if the key is not found in the config file
  int nmppcy                    = config.read<int>("nmppcy",2);                     // number of mppc in y direction per mppc - default to 2 if the key is not found in the config file
  int nmodulex                  = config.read<int>("nmodulex",1);                   // number of modules in x direction per mppc - default to 1 if the key is not found in the config file
  int nmoduley                  = config.read<int>("nmoduley",1);                   // number of modules in y direction per mppc - default to 1 if the key is not found in the config file
  int histo1Dmax                = config.read<int>("histo1Dmax",25000);                   // max of the 1D charge histograms (in ADC channels)
  int histo1Dbins               = config.read<int>("histo1Dbins",250);                  // number of bins of the 1D charge histograms
  int histoLYmax                = config.read<int>("histoLYmax",20000);                   // max of the 1D charge histograms (in ADC channels)
  int histoLYbins               = config.read<int>("histoLYbins",500);                  // number of bins of the 1D charge histograms
  double qe                     = config.read<double>("qe",0.35);                   // mppc QE
  double gainMPPC               = config.read<double>("gainMPPC",1.25e6);           // mppc gain
  double chargeBinningADC       = config.read<double>("chargeBinningADC",156e-15);  // adc charge binning
  double sourceMeV              = config.read<double>("sourceMeV",0.511);           // gamma peak in MeV
  int histo2DchannelBin         = config.read<int>("histo2DchannelBin",250);            // number of bins of the 2D flood histograms, for single channels
  int histo2DglobalBins         = config.read<int>("histo2DglobalBins",1000);            // number of bins of the 2D flood histograms, for entire module
  int histo3DglobalBins         = config.read<int>("histo3DglobalBins",100);            // number of bins of the 3D flood histograms, for entire module
  int taggingPeakMin            = config.read<int>("taggingPeakMin",0);          // min range of tagging crystal photopeak, in ADC channels - to help TSpectrum
  int taggingPeakMax            = config.read<int>("taggingPeakMax",0);         // max range of tagging crystal photopeak, in ADC channels - to help TSpectrum
  // int clusterLevelPrecision     = config.read<int>("clusterLevelPrecision",10);     // precision of the level search when separating the cluster of 3D points
  float taggingPosition         = config.read<float>("taggingPosition",0);            // position of the tagging bench in mm
  bool usingTaggingBench        = config.read<bool>("usingTaggingBench",0);           // true if the input is using tagging bench, false if not
  int taggingCrystalChannel     = config.read<int>("taggingCrystalChannel",16);        // input channel where the tagging crystal information is stored
  bool calcDoiResWithDelta      = config.read<bool>("calcDoiResWithDelta",0);         // alternative calcolation of doi res, based on deltas. only if it's usingTaggingBench
  std::string calcDoiFileName   = config.read<std::string>("calcDoiFileName","");      // name (and path if not in this folder) of calibration_params.txt
  int pointsFromDoi             = config.read<int>("pointsFromDoi",0);                 // points measured in file calcDoiFileName. this is MANDATORY if calcDoiFileName is specified!
  bool correctingSaturation     = config.read<bool>("correctingSaturation");        // true if saturation correction is applied, false if it's not
  bool correctingForDOI         = config.read<bool>("correctingForDOI",0);              // true if the energy correction using DOI info is computed
  float energyResolution        = config.read<float>("expectedEnergyResolution",0);     // energy resolution input by the user, if any, otherwise 0
  bool usingRealSimData         = config.read<bool>("usingRealSimData",0);
  float moduleLateralSideX      = config.read<float>("moduleLateralSideX",7.0);         //
  float moduleLateralSideY      = config.read<float>("moduleLateralSideY",7.0);         //
  bool backgroundRun            = config.read<bool>("backgroundRun",0);                 // whether this is a background run or not
  bool lateralRun               = config.read<bool>("lateralRun",0);                    // whether this is a lateral irradiation run or not
  float userBroadCut            = config.read<float>("userBroadCut",0.0);            // if in backgroundRun, cut to get rid of low energy events is not done on photopeak search but by user input (default 0ch)
  float thresholdKev            = config.read<float>("thresholdKev",0.0);
  float wThreshold              = config.read<float>("wThreshold",0.1);                 // Threshold for w plots limits
  double crystalz               = config.read<double>("crystalz",15);
  // double DoiResolutionVsIJmax   = config.read<double>("DoiResolutionVsIJmax",10);       // max of the 2d DoiResolution values plot (starts from 0) - it's mm
  double EnergyResolutionVsIJmax= config.read<double>("EnergyResolutionVsIJmax",0.3);   // max of the 2d EnergyResolution values plot (starts from 0)
  double LYvsIJmax              = config.read<double>("LYvsIJmax",40000);
  double PeakPositionVsIJmax    = config.read<double>("PeakPositionVsIJmax",12000);     // max of the 2d PeakPosition values plot (starts from 0)  - it's ADC channels
  int wHistogramsBins           = config.read<int>("wHistogramsBins",250);
  int doiColumnOffset           = config.read<int>("doiColumnOffset",0);                // for DOI output, fix the column i by adding this quantity. if not stated, 0 by default
  double energyCorrectionMin    = config.read<double>("energyCorrectionMin",0.25);      // once the wmin and wmax are found for each w histo, choose at which point to start and to stop the linear fitting
  double energyCorrectionMax    = config.read<double>("energyCorrectionMax",0.75);      // (as percentage from min to max)
  double lambda511              = config.read<double>("lambda511",12.195); //everything in mm
  // bool wAllChannels             = config.read<bool>("wAllChannels",0);                  // whether we use the sum of all channels to compute w (true = 1) of just the neighbours (false = 0). Deafult to false.
  bool comptonAnalysis          = config.read<bool>("comptonAnalysis",0);               //whether to perform or not the compton recovery analysis part. Default to false
  bool lightYieldComputation    = config.read<bool>("lightYieldComputation",0);         //whether to perform or not the light yield calculation. Default to false
  float peakSearchRangeMin      = config.read<int>("peakSearchRangeMin",0);             //lower limit for search of 511KeV peak - wide limitation to help peak search
  float peakSearchRangeMax      = config.read<int>("peakSearchRangeMax",histo1Dmax);    //upper limit for search of 511KeV peak - wide limitation to help peak search
  bool saturationRun            = config.read<bool>("saturationRun",0);                 //whether this is a saturation run or not
  float histoSingleChargeMax    = config.read<float>("histoSingleChargeMax",0);         // max in the histograms of charge for saturationRun
  float histoSingleChargeBin    = config.read<float>("histoSingleChargeBin",0);         // max in the histograms of charge for saturationRun
  float histoSumChargeMax    = config.read<float>("histoSumChargeMax",0);         // max in the histograms of charge for saturationRun
  float histoSumChargeBin    = config.read<float>("histoSumChargeBin",0);         // max in the histograms of charge for saturationRun
  float saturationPeakFractionLow    = config.read<float>("saturationPeakFractionLow",0.06);         // lower limit for saturation peak fitting, expressed in fraction of peak position. limit will be = (peak - peak*saturationPeakFractionLow)
  float saturationPeakFractionHigh    = config.read<float>("saturationPeakFractionHigh",0.06);         // upper limit for saturation peak fitting, expressed in fraction of peak position. limit will be = (peak + peak*saturationPeakFractionHigh)
  bool performSaturationPeakSearch          = config.read<bool>("performSaturationPeakSearch",1);    //perform of nor satuartion peak search
  bool backgroundSaturationRun              = config.read<bool>("backgroundSaturationRun",0);
  // set output file name
  std::string outputFileName = config.read<std::string>("output");
  std::string digitizer_s   = config.read<std::string>("digitizer");
  std::string saturationPeakEnergy_s = config.read<std::string>("saturationPeakEnergy","");
  std::string saturationPeakMin_s = config.read<std::string>("saturationPeakMin","");
  std::string saturationPeakMax_s = config.read<std::string>("saturationPeakMax","");
  std::string loadAnalysisTreeName         = config.read<std::string>("loadAnalysisTreeName","0");       // look for input analysis ttree in config file. if no input, it will produced by this program
  std::string saveAnalysisTreeName         = config.read<std::string>("saveAnalysisTreeName","0");       // look for filename to save analysis ttree in config file. if no filename, analysis tree won't be saved
  bool saveAnalysisTree          = config.read<bool>("saveAnalysisTree",0);            // choice to save or not the analysis TTree, in a file (name chosen above)
  bool calcDoiResWithCalibration = config.read<bool>("calcDoiResWithCalibration",0);
  std::string calibrationFileName       = config.read<std::string>("calibrationFileName","0");
  bool usingAllChannelsForEnergySpectra = config.read<bool>("usingAllChannelsForEnergySpectra",0);      // whether to use all digitizer channels for energy spectra (1) or just neighbours (0)
  bool usingAllChannels            = config.read<bool>("usingAllChannels",0);
  bool wAllChannels                = config.read<bool>("wAllChannels",0);                  // whether we use the sum of all channels to compute w (true = 1) of just the neighbours (false = 0). Deafult to false.
  int   taggingCrystalBins = config.read<int>("taggingCrystalBins",1200);
  float taggingSpectrumMin = config.read<float>("taggingSpectrumMin",0.0);
  float taggingSpectrumMax = config.read<float>("taggingSpectrumMax",12000.0);
  bool TagEdgeCalculation  = config.read<bool>("tagEdgeCalculation",0);
  int digitizerType               = config.read<int>("digitizerType",0);       // type of digitizer. 0 = desktop, 1 = vme, 2 = petiroc
  bool cuttingOnTagPhotopeak = config.read<bool>("cuttingOnTagPhotopeak",1);
  int CTRbins = config.read<int>("CTRbins",500);
  float CTRmin = config.read<float>("CTRmin",-5e-9);
  float CTRmax = config.read<float>("CTRmax",5e-9);
  float histo3Dmin = config.read<float>("histo3Dmin",0);
  float histo3Dmax = config.read<float>("histo3Dmax",1);
  int DeltaTimeBins  = config.read<int>("DeltaTimeBins",500);
  float DeltaTimeMin = config.read<float>("DeltaTimeMin",-5e-9);
  float DeltaTimeMax = config.read<float>("DeltaTimeMax",5e-9);
  bool taggingForTiming = config.read<bool>("taggingForTiming",0);
  float photopeakSigmasMin = config.read<float>("photopeakSigmasMin",2.0);      // how many sigmas far from mean is the lower bound of cut on photopeak  - default = 2.0
  float photopeakSigmasMax = config.read<float>("photopeakSigmasMax",4.0);     // how many sigmas far from mean is the upper bound of cut on photopeak  - default = 4.0
  float TaggingPhotopeakSigmasMin = config.read<float>("TaggingPhotopeakSigmasMin",1.5);      // how many sigmas far from mean is the lower bound of cut on Tagging photopeak  - default = 1.5
  float TaggingPhotopeakSigmasMax = config.read<float>("TaggingPhotopeakSigmasMin",2.0);     // how many sigmas far from mean is the upper bound of cut on Tagging photopeak  - default = 2.0
  int WrangeBinsForTiming = config.read<int>("WrangeBinsForTiming",10);
  bool smearTaggingTime = config.read<bool>("smearTaggingTime",0);// whether to smear the time stamp of external tagging. Needed for simulations, where the tagging time stamp is always 0 (i.e. the gamma emission time) - default = 0
  float tagFitLowerFraction = config.read<float>("tagFitLowerFraction",0.06);  // enRes = 2.355*sigma/peak --> sigma = enRes*Peak/2.355   EnRes about 0.15 -->  sigma = 0.06* peak  -> limits -1sigma +2 sigma
  float tagFitUpperFraction = config.read<float>("tagFitUpperFraction",0.12);
  bool timingCorrection = config.read<bool>("timingCorrection",1);// perform or not the timing correction (for example, it could make no sense to perform it fully in polished arrays) - deafult = 1
  bool timingCorrectionForPolished = config.read<bool>("timingCorrectionForPolished",1);// produce the plots for the simple combination of timestamps that can be done for polished crystals (i.e. when DOI info is not available)
  int adcChannels = config.read<int>("digitizerTotalCh");// number of output channels in the adc in use (32 for the CAEN DT5740)
  float noiseSigmas =  config.read<float>("noiseSigmas",1.0); // how many sigmas of noise far from "0" to cut the dataset
  int taggingCrystalTimingChannel  = config.read<int>("taggingCrystalTimingChannel",16);     // input timing channel where the tagging crystal information is stored                - default = 16
  float marginWZgraph = config.read<float>("marginWZgraph",0.1); // distance from crystal limit for calculation of beginW and endW
  float tagCrystalPeakResolutionFWHM = config.read<float>("tagCrystalPeakResolutionFWHM",0.07);
  float photopeakFitRangeMin = config.read<float>("photopeakFitRangeMin",1.2);
  float photopeakFitRangeMax = config.read<float>("photopeakFitRangeMax",1.2);

  //----------------------------------------------------------//
  //  Load and save TTree                                     //
  //----------------------------------------------------------//
  //decide whether to load the ttree from file or produce it, and whether to save it or not to file
  // bool loadAnalysisTree;
  // if(loadAnalysisTreeName.compare("0") == 0) // no file provided to load ttree
  //   loadAnalysisTree = false;
  // else  // file provided
  //   loadAnalysisTree = true;
  // if(!saveAnalysisTree) // user didn't set to save ttree
  // {
  //   if(saveAnalysisTreeName.compare("0") != 0) // but gave a name to save the ttree
  //   {
  //     saveAnalysisTree = true;  // then set the ttree to be saved
  //   }
  // }
  // else //user set the flag to save ttree
  // {
  //   if(saveAnalysisTreeName.compare("0") == 0) // but didn't give a name to save the ttree
  //   {
  //     saveAnalysisTreeName = "temp.root"; // then set the ttree to be saved into temp.root file
  //   }
  // }
  // if(loadAnalysisTree) // anyway there's no point saving the analysis tree if it's not produced by this analysis but just loaded from file
  // {
  //   saveAnalysisTree = false;
  // }
  // std::cout<<"\n"<<std::endl;
  // std::cout<<"###########################################################"<<std::endl;
  // std::cout<<"#                                                         #"<<std::endl;
  // if(loadAnalysisTree)
  // std::cout << "# Analysis TTree loaded from file " << loadAnalysisTreeName.c_str() << std::endl;
  // else
  // std::cout << "# Analysis TTree produced by this program " << std::endl;
  // if(saveAnalysisTree)
  // std::cout << "# Analysis TTree will be saved to file " << saveAnalysisTreeName.c_str() << std::endl;
  // else
  // std::cout << "# No analysis TTree file will be saved " << std::endl;
  // std::cout<<"#                                                         #"<<std::endl;
  // std::cout<<"###########################################################"<<std::endl;
  // std::cout<<"\n"<<std::endl;
  // if(!loadAnalysisTree)
  // {
  //   if(digitizerType == 0 || digitizerType == 1)
  //   {
  //     input.ImportTChain(argc,argv);
  //   }
  //   // input.PrepareTTree();
  //   // input.FillTree(argc,argv);                        // create the TTree that will be used in analysis
  // }
  // else    // otherwise load it from the indicated file
  // {
  //   TFile *fTemp = new TFile(loadAnalysisTreeName.c_str());
  //   TTree *TempTree =  (TTree*) fTemp->Get("adc");
  //   input.SetTree(TempTree);
  // }

  // GET INFO ON PC AND FILES
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  std::string HostNameString(hostname);
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::string PWDstring(cwd);
  // std::cout << cwd << std::endl;

  ROOT::v5::TFormula::SetMaxima(10000,10000,10000);

  // list of all input files
  std::stringstream InputFiles;
  if(std::string(argv[1]) == std::string("-c")) // first argument is -c, then the config file name is passed by command line
  {
    for (int i = 3; i < argc ; i++) // run on the remaining arguments to add all the input files
    {
      // std::cout << "Adding file " << argv[i] << std::endl;
      InputFiles << argv[i] << std::endl;
    }
  }
  else // the config file was indeed the default one
  {
    for (int i = 1; i < argc ; i++) // run on the remaining arguments to add all the input files
    {
      // std::cout << "Adding file " << argv[i] << std::endl;
      InputFiles << argv[i] << std::endl;
    }
  }
  // std::cout << InputFiles.str() << std::endl;


  //set outputFileName
  // change name if parallel keys are given
  std::string parallelOutput = config.read<std::string>("parallelOutput","0");
  if (parallelOutput.compare("0") != 0)
  {
    outputFileName = parallelOutput;
  }
  outputFileName += ".root";


  TFile *calibrationFile = NULL;
  if(calcDoiResWithCalibration)
  {
    calibrationFile = new TFile(calibrationFileName.c_str());
  }

  //----------------------------------------------------------//
  //  Creating the module and elements                        //
  //----------------------------------------------------------//
  // create the elements of the module
  Module***   module;
  Mppc***     mppc;
  Crystal***  crystal;
  // make an array of module pointers
  module = new Module**[nmodulex];
  for(int j = 0; j < nmodulex ; j++) module[j] = new Module* [nmoduley];
  // make an array of mppc pointers
  mppc = new Mppc**[nmodulex*nmppcx];
  for(int j = 0; j < nmodulex*nmppcx ; j++) mppc[j] = new Mppc* [nmoduley*nmppcy];
  // make an array of crystal pointers
  crystal = new Crystal**[nmodulex*nmppcx*ncrystalsx];
  for(int j = 0; j < nmodulex*nmppcx*ncrystalsx ; j++) crystal[j] = new Crystal* [nmoduley*nmppcy*ncrystalsy];
  // fill the elements
  input.FillElements(module,mppc,crystal);
  //----------------------------------------------------------//



  //----------------------------------------------------------//
  //  Plots and spectra                                       //
  //----------------------------------------------------------//

  //properly set timing correction flags
  // the problem is back compatibility
  // timingCorrection is used to enable the complete time correciton analysis
  // timingCorrectionForPolished enables only the correction for polished crystals
  // but the second is part of the first, so the second controls and if statement that
  // would not allow the first to happen. Therefore, if the first is set on but the second is off
  // the second needs to be enabled too
  if(timingCorrection == true && timingCorrectionForPolished == false)
  {
    timingCorrectionForPolished = true;
  }

  // doi bench specific part
  std::ofstream doiFile;
  std::ofstream AltDoiFile;
  std::ifstream altDoiCalcFile;

  // saturation dataset part
  std::ofstream saturationFile;
  if(usingTaggingBench && calcDoiResWithDelta)
  {
    if( (calcDoiFileName.compare("") == 0)  )
    {
      std::cout << "No file with doi Calibration provided!" << std::endl;
      calcDoiResWithDelta = false;
    }
    if( pointsFromDoi == 0)
    {
      std::cout << "Need to specify the number of points in doi calibration file!" << std::endl;
      calcDoiResWithDelta = false;
    }
  }
  //saturation run part
  std::vector<SaturationPeak_t> saturationPeak;
  if(saturationRun)
  {
    saturationFile.open("saturationData.txt", std::ofstream::out);
    if(histoSingleChargeMax == 0)// calc histoSingleChargeMax from the histo1Dmax,
    {
      histoSingleChargeMax = histo1Dmax * chargeBinningADC / 2.0 ; // the /2.0 it's because the histo1D is sum spectrum
    }
    if(histoSingleChargeBin == 0)// set the histoSingleChargeBin to histo1Dbins
    {
      histoSingleChargeBin = histo1Dbins;
    }
    if(histoSumChargeMax == 0)// set it to the same as single
    {
      histoSumChargeMax = histoSingleChargeMax ; // the /2.0 it's because the histo1D is sum spectrum
    }
    if(histoSumChargeBin == 0)// set it to the same as single
    {
      histoSumChargeBin = histoSingleChargeBin;
    }
    //read the config file to set the search areas for peaks
    // std::string                 saturationPeakEnergy_s,saturationPeakMin_s,saturationPeakMax_s;
    std::vector <std::string>   saturationPeakEnergy_f,saturationPeakMin_f,saturationPeakMax_f;
    std::vector <float>         saturationPeakEnergy,saturationPeakMin,saturationPeakMax;
    config.split( saturationPeakEnergy_f, saturationPeakEnergy_s, "," );
    config.split( saturationPeakMin_f, saturationPeakMin_s, "," );
    config.split( saturationPeakMax_f, saturationPeakMax_s, "," );
    for(unsigned int i = 0 ; i < saturationPeakEnergy_f.size() ; i++)
    {
      config.trim(saturationPeakEnergy_f[i]);
      saturationPeakEnergy.push_back(atof(saturationPeakEnergy_f[i].c_str()));
    }
    for(unsigned int i = 0 ; i < saturationPeakMin_f.size() ; i++)
    {
      config.trim(saturationPeakMin_f[i]);
      saturationPeakMin.push_back(atof(saturationPeakMin_f[i].c_str()));
    }
    for(unsigned int i = 0 ; i < saturationPeakMax_f.size() ; i++)
    {
      config.trim(saturationPeakMax_f[i]);
      saturationPeakMax.push_back(atof(saturationPeakMax_f[i].c_str()));
    }
    // store data in saturationPeak vector of struct
    for(unsigned int i = 0 ; i < saturationPeakEnergy.size() ; i++)
    {
      SaturationPeak_t tempPeak;
      tempPeak.energy  = saturationPeakEnergy[i];
      tempPeak.peakMin = saturationPeakMin[i];
      tempPeak.peakMax = saturationPeakMax[i];
      saturationPeak.push_back(tempPeak);
    }
  }

  //points form DOI input file
  std::vector<inputDoi_t> inputDoi;
  inputDoi_t tempInputDoi(pointsFromDoi);
  if(usingTaggingBench)
  {
    doiFile.open("doiData.txt", std::ofstream::out);
    if(calcDoiResWithDelta)
    {
      altDoiCalcFile.open(calcDoiFileName.c_str(), std::ios::in);
      if(!altDoiCalcFile)
      {
        std::cout << "Doi calibration file not found!" << std::endl;
        calcDoiResWithDelta = false;
      }
      else //finally, everything is fine...
      {
        AltDoiFile.open("altDoiRes.txt", std::ofstream::out);
        while(altDoiCalcFile >> tempInputDoi)
        {
          inputDoi.push_back(tempInputDoi);
          tempInputDoi.clear();
        }
      }
    }
  }

  // simulation dataset spectific part
  TCut SingleCrystalInteraction;
  TCut SingleEnergyDeposition;
  if(usingRealSimData) // only if this is a sim dataset
  {
    SingleCrystalInteraction = "CrystalsHit == 1";
    SingleEnergyDeposition = "NumbOfInteractions == 1";
  }

  //variable that will be used only if the user requires the TagEdgeCalculation
  double tagPeakHgEntries = 0.0;

  // MAIN LOOP v2
  // get the TTree, to plot the spectra
  // TTree* tree = input.GetTree();
  // fill the tchain with input files

  //----------------------------------------------------------//
  //  Get TChain of input files                               //
  //----------------------------------------------------------//
  TChain* tree = new TChain(fname.c_str());  // create the input tchain and the analysis ttree
  // first, create the adc channels variables and branches
  ChainAdcChannel        = new Int_t [adcChannels];
  ChainDesktopAdcChannel = new Short_t [adcChannels]; // input from ADC desktop
  ChainVMEadcChannel     = new UShort_t [adcChannels]; // input from VME
  ChainTimeStamp         = new Float_t[adcChannels];
  TDCBinning             = new Float_t[adcChannels];
  // DigitizerChannelOn     = new bool[adcChannels];
  bChainAdcChannel       = new TBranch* [adcChannels];
  bChainTimeStamp        = new TBranch* [adcChannels];
  if(std::string(argv[1]) == std::string("-c")) // first argument is -c, then the config file name is passed by command line
  {
    for (int i = 3; i < argc ; i++) // run on the remaining arguments to add all the input files
    {
      std::cout << "Adding file " << argv[i] << std::endl;
      tree->Add(argv[i]);
    }
  }
  else // the config file was indeed the default one
  {
    for (int i = 1; i < argc ; i++) // run on the remaining arguments to add all the input files
    {
      std::cout << "Adding file " << argv[i] << std::endl;
      tree->Add(argv[i]);
    }
  }
  // set branches for reading the input files
  tree->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag, &bChainExtendedTimeTag);
  tree->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag, &bChainDeltaTimeTag);
  if(usingRealSimData)
  {
    tree->SetBranchAddress("RealX", &RealX, &bRealX);
    tree->SetBranchAddress("RealY", &RealY, &bRealY);
    tree->SetBranchAddress("RealZ", &RealZ, &bRealZ);
    // fchain->SetBranchAddress("Tagging", &simTaggingCharge, &bsimTaggingCharge);
    // fchain->SetBranchAddress("TaggingTimeStamp", &simTaggingTime, &bsimTaggingTime);
    tree->SetBranchAddress("CrystalsHit",&CrystalsHit, &bCrystalsHit);
    tree->SetBranchAddress("NumbOfInteractions",&NumbOfInteractions, &bNumbOfInteractions);
    // fchain->SetBranchAddress("TotalCryEnergy",&TotalCryEnergy, &bTotalCryEnergy);
  }
  for(int i=0; i<adcChannels; i++)
  {
    if(digitizerType == 0)
    {
      std::stringstream sname;
      sname << "ch" << i;
      tree->SetBranchAddress(sname.str().c_str(), &ChainDesktopAdcChannel[i], &bChainAdcChannel[i]);
    }
    if(digitizerType == 1)
    {
      std::stringstream sname;
      sname << "ch" << i;
      tree->SetBranchAddress(sname.str().c_str(), &ChainVMEadcChannel[i], &bChainAdcChannel[i]);
      sname.str("");
      sname << "t" << i;
      tree->SetBranchAddress(sname.str().c_str(), &ChainTimeStamp[i],&bChainTimeStamp[i]);
    }
  }

  TH1F* TaggingCrystalSpectrum = NULL;
  TH1F *TriggerSpectrumHighlight = NULL;



  // Loop on modules, mppcs and crystal
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      // different variables logic, trying to use directly the original files:
      // almost everything need to be produced channel by channel (otherwise the
      // definition of neighbours is impossible)
      // so the variables in the original files are just:
      //
      // - BoardTimeTags ... ch0 ... chN , t0 ... tN
      //
      // and with the pre-processing we were genetating some additional variables:
      //
      // - TriggerChannel   = the channel with highest charge among all channels in the event
      // - FloodX,Y         = weighted average on charge of detector positions, with multiple coices (i.e. summing on all detectors or only on trigger channel + 8 neighbours)
      // - FloodZ           = ratio (charge in trigger channel)/(sum charge) where sum chage can be the sum of all channels or the sum of trigger channel + 8 neighbours
      // - BadEvent         = counter of events with charge higher of saturation parameter for that specific channel. these events are "bad" because given the logarithmic formula of saturation correction, it would be mathematically impossible to correct them. It's just a counter that provides a feedback to user (typically if all makes sense, badEvents need to be just a small fraction of all events)
      // - Tagging          = charge measured in the external tagging crystal, if any
      // - TaggingTimeStamp = time stamp of the external tagging crystal, if any
      // - ZPosition        = z position of the tagging bench, for DOI measurements, for this event
      //
      // These variables can be computed on the fly directly on the original data (with some limitations, see later)
      //



      // write down the channel collection

      // std::stringstream TriggerChannel;

      // TriggerChannel


      //useful strings etc
      // First cut is a general, global cut
      // it will be composed of all cuts that will apply generally to all the events
      // For example
      // 1. A global XYZ cut that allows us to consider only the events with "reasonable" uvw coordinates
      // 2. A cut on photopeak of the tagging crystal, if there is any such crystal involved
      std::cout << "Generating global spectra..." << std::endl;


      std::stringstream var;
      std::stringstream sname;
      std::stringstream cut;
      TCut BasicCut;

      //FIXME what do we do with the "old" CutXYZ?   <--------------------
      // std::stringstream CutXYZstream;
      // CutXYZstream << "FloodX > " << -moduleLateralSideX << " && FloodX < " << moduleLateralSideX  << "&& FloodY > " << -moduleLateralSideY <<   " && FloodY < " << moduleLateralSideY << "&& FloodZ > " << histo3Dmin << " && FloodZ <  " << histo3Dmax;
      // TCut CutXYZ = CutXYZstream.str().c_str();

      // get a std::vector of all module channels
      // std::vector<int> moduleChannels = module[iModule][jModule]->GetChannels();
      //get the detector data
      std::vector<detector_t> detector;
      detector = module[iModule][jModule]->GetDetector();

      // // //DEBUG OUTPUT
      // // output the detector[] data for check
      // // struct detector_t
      // // // {
      // // //   int digitizerChannel;
      // // //   int timingChannel;
      // // //   std::string label;
      // // //   int i;
      // // //   int j;
      // // //   float saturation;
      // // //   int plotPosition;
      // // //   float xPosition;
      // // //   float yPosition;
      // // //   float pedestal;
      // // //   int OnForDOI;
      // // //   std::vector<int> neighbourChannels;
      // // //   bool OnForModular;
      // // //   //     bool operator<(const masks_t& rhs) const { meanx < rhs.meanx; }
      // // // };
      // std::cout << "////////////////////////////////////////////////" << std::endl;
      // std::cout << "//             Detector Structure             //" << std::endl;
      // std::cout << "////////////////////////////////////////////////" << std::endl;
      // for(unsigned int iDet = 0; iDet < detector.size() ; iDet++)
      // {
      //   std::cout << "------------------------------------------------" << std::endl;
      //   std::cout << "Index              = " << iDet << std::endl;
      //   std::cout << "digitizerChannel   = " << detector[iDet].digitizerChannel << std::endl;
      //   std::cout << "timingChannel      = " << detector[iDet].timingChannel << std::endl;
      //   std::cout << "label              = " << detector[iDet].label << std::endl;
      //   std::cout << "i                  = " << detector[iDet].i << std::endl;
      //   std::cout << "j                  = " << detector[iDet].j << std::endl;
      //   std::cout << "saturation         = " << detector[iDet].saturation << std::endl;
      //   std::cout << "plotPosition       = " << detector[iDet].plotPosition << std::endl;
      //   std::cout << "xPosition          = " << detector[iDet].xPosition << std::endl;
      //   std::cout << "yPosition          = " << detector[iDet].yPosition << std::endl;
      //   std::cout << "pedestal           = " << detector[iDet].pedestal << std::endl;
      //   std::cout << "OnForDOI           = " << detector[iDet].OnForDOI << std::endl;
      //   std::cout << "Neighbours         = ";
      //   for(unsigned int iNeigh = 0; iNeigh < detector[iDet].neighbourChannels.size(); iNeigh++)
      //   {
      //     std::cout << detector[iDet].neighbourChannels[iNeigh] << " " ;
      //   }
      //   std::cout << std::endl;
      //   std::cout << "OnForModular       = " << detector[iDet].OnForModular << std::endl;
      //   std::cout << std::endl;
      // }
      // std::cout << "////////////////////////////////////////////////" << std::endl;
      // std::cout << "//         End of Detector Structure          //" << std::endl;
      // std::cout << "////////////////////////////////////////////////" << std::endl;
      // // //END of DEBUG OUTPUT

      std::vector<int> allModuleChID;
      std::vector<float> saturationData;
      std::vector<float> pedestalData;
      std::vector<int> detChData;
      for(unsigned int iDet = 0; iDet < detector.size() ; iDet++)
      {
        allModuleChID.push_back(iDet); //dummy... we need to get all the indexes of detector[] that identify all the channels in the array, and of course these are all the indexes of the detector[] vector.
        detChData.push_back(detector[iDet].digitizerChannel);
        saturationData.push_back(detector[iDet].saturation);
        pedestalData.push_back(detector[iDet].pedestal);
      }

      module[iModule][jModule]->SetDetChData(detChData);
      module[iModule][jModule]->SetSaturationData(saturationData);
      module[iModule][jModule]->SetPedestalData(pedestalData);


      // //DEBUG OUTPUT
      // std::cout << std::endl;
      // for(unsigned int iDet = 0; iDet < moduleChannels.size() ; iDet++)
      // {
      //   std::cout << moduleChannels[iDet] << " ";
      // }
      // std::cout << std::endl;


      // BasicCut += CutXYZ;  // add the module physical constrains the to global cut of accepted events
      //
      TCut taggingPhotopeakCut = "" ;

      // int Tagging =

      // prepare GLOBAL SPECTRA
      // Flood histogram
      // TString nameModule;
      // std::stringstream varModule;
      sname << "Flood Histogram 2D - " << module[iModule][jModule]->GetName();
      // varModule << "FloodY:FloodX >> " << nameModule;
      TH2F *spectrum2dModule = new TH2F(sname.str().c_str(),sname.str().c_str(),histo2DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo2DglobalBins,-moduleLateralSideY,moduleLateralSideY);
      // tree->Draw(varModule.str().c_str(),"","COLZ");
      spectrum2dModule->SetName(sname.str().c_str());
      spectrum2dModule->SetTitle(sname.str().c_str());
      spectrum2dModule->GetXaxis()->SetTitle("U");
      spectrum2dModule->GetYaxis()->SetTitle("V");
      // module[iModule][jModule]->SetFloodMap2D(spectrum2dModule);
      sname.str("");





      if(usingTaggingBench || calcDoiResWithCalibration || taggingForTiming)//trigger spectrum
      {
        if(taggingPeakMax == 0)
        {
          taggingPeakMax = taggingSpectrumMax;
        }
        if(taggingPeakMin == 0)
        {
          taggingPeakMin = taggingSpectrumMin;
        }

        std::stringstream tagStream;
        tagStream << "ch" << taggingCrystalChannel; // if using original files

        TaggingCrystalSpectrum =  new TH1F("TaggingCrystalSpectrum","TaggingCrystalSpectrum",taggingCrystalBins,taggingSpectrumMin,taggingSpectrumMax);
        var << tagStream.str().c_str() << ">> TaggingCrystalSpectrum";
        // 	std::cout << nameModule << " ... ";
        TaggingCrystalSpectrum->SetName("TaggingCrystalSpectrum");
        TaggingCrystalSpectrum->GetXaxis()->SetTitle("ADC");
        TaggingCrystalSpectrum->GetYaxis()->SetTitle("Counts");
        TaggingCrystalSpectrum->SetTitle("Spectrum of Tagging Crystal");
        tree->Draw(var.str().c_str(),"");
        var.str("");


        //restrict the region where to look for peaks. Fix for tspectrum...
        TaggingCrystalSpectrum->GetXaxis()->SetRangeUser(taggingPeakMin,taggingPeakMax);
        //find peak in the tagging crystal
        TSpectrum *sTagCrystal;
        sTagCrystal = new TSpectrum(5);
        Int_t TagCrystalPeaksN = sTagCrystal->Search(TaggingCrystalSpectrum,1,"",0.5);
        Double_t *TagCrystalPeaks  = sTagCrystal->GetPositionX();
        Double_t *TagCrystalPeaksY = sTagCrystal->GetPositionY();
        // float saturationPeakFraction
        //delete s;
        // float distPeak = INFINITY;
        float maxPeak = 0.0;
        int peakID = 0;
        for (int peakCounter = 0 ; peakCounter < TagCrystalPeaksN ; peakCounter++ )
        {
          // //DEBUG
          // std::cout << peakCounter << "\t" << TagCrystalPeaks[peakCounter] << "\t" << TagCrystalPeaksY[peakCounter] << std::endl;
          // if( fabs(CrystalPeaks[peakCounter] - 0.5*(saturationPeak[iSaturation].peakMin+saturationPeak[iSaturation].peakMax)) < distPeak)//take closest peak to the center of search range selected by the user
          if(TagCrystalPeaksY[peakCounter] > maxPeak)
          {
            // distPeak = CrystalPeaks[peakCounter];
            maxPeak = TagCrystalPeaksY[peakCounter];
            peakID = peakCounter;
          }
        }
        // //DEBUG
        // std::cout << "chosen " << peakID << "\t" << TagCrystalPeaks[peakID] << "\t" << TagCrystalPeaksY[peakID] << std::endl;
        TF1 *gaussTag = new TF1("gaussTag",
                                "gaus",
                                TagCrystalPeaks[peakID] - tagFitLowerFraction*TagCrystalPeaks[peakID],
                                TagCrystalPeaks[peakID] + tagFitUpperFraction*TagCrystalPeaks[peakID]);


        // TF1 *gaussTag = new TF1("gaussTag",
        //                         "gaus",
        //                         taggingPeakMin,
        //                         taggingPeakMax);

        gaussTag->SetParameter(0,TagCrystalPeaksY[peakID]); // heigth
        gaussTag->SetParameter(1,TagCrystalPeaks[peakID]); //mean
        gaussTag->SetParameter(2,(tagCrystalPeakResolutionFWHM*TagCrystalPeaks[peakID])/2.355); // sigma
        //

        // TaggingCrystalSpectrum->Fit("gaussTag",
        //                             "Q",
        //                             "",
        //                             taggingPeakMin,
        //                             taggingPeakMax);

        TaggingCrystalSpectrum->Fit("gaussTag",
                                    "Q",
                                    "",
                                    TagCrystalPeaks[peakID] - tagFitLowerFraction*TagCrystalPeaks[peakID],
                                    TagCrystalPeaks[peakID] + tagFitUpperFraction*TagCrystalPeaks[peakID]);
        // // DEBUG
        // std::cout << "fit pars " << gaussTag->GetParameter(1) << "\t" << gaussTag->GetParameter(2) << std::endl;
        TaggingCrystalSpectrum->GetXaxis()->SetRangeUser(taggingSpectrumMin,taggingSpectrumMax);
        //define a TCut for this peak
        double tagPhotopeakMin = gaussTag->GetParameter(1) - TaggingPhotopeakSigmasMin*gaussTag->GetParameter(2);
        double tagPhotopeakMax = gaussTag->GetParameter(1) + TaggingPhotopeakSigmasMax*gaussTag->GetParameter(2);
        std::stringstream tagString;
        tagString << tagStream.str().c_str() << " > " << tagPhotopeakMin << "&& " << tagStream.str().c_str() <<" < " << tagPhotopeakMax;
        taggingPhotopeakCut = tagString.str().c_str();
        //highlighted spectrum
        TriggerSpectrumHighlight = new TH1F("TriggerSpectrumHighlight","",taggingCrystalBins,taggingSpectrumMin,taggingSpectrumMax);
        var.str("");
        var << tagStream.str().c_str() << " >> TriggerSpectrumHighlight";
        TriggerSpectrumHighlight->SetLineColor(3);
        TriggerSpectrumHighlight->SetFillColor(3);
        TriggerSpectrumHighlight->SetFillStyle(3001);
        tree->Draw(var.str().c_str(),taggingPhotopeakCut);
        var.str("");
        if(TagEdgeCalculation)
        {
          tagPeakHgEntries = TriggerSpectrumHighlight->GetEntries();
        }
        if(cuttingOnTagPhotopeak) // only if tagging on photopeak is chosen
        {
          BasicCut += taggingPhotopeakCut;  // add the cut on photopeak of the tagging crystal.
        }
        //but always save the tagging cut in the output file
        taggingPhotopeakCut.SetName("taggingPhotopeakCut");
        module[iModule][jModule]->SetTaggingPhotopeakCut(taggingPhotopeakCut);        // energy cut, events in the photopeak
        // 	delete sTagCrystal;
        delete gaussTag;
        // 	std::cout << " done" << std::endl;


        //save tagging crystal timing channel
        module[iModule][jModule]->SetTaggingTimingChannel(taggingCrystalTimingChannel);
      }

      sname << "Flood Histogram 3D - Module " << module[iModule][jModule]->GetName();
      TH3I* spectrum3dModule = new TH3I(sname.str().c_str(),sname.str().c_str(),histo3DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo3DglobalBins,-moduleLateralSideY,moduleLateralSideY,histo3DglobalBins,histo3Dmin,histo3Dmax);
      spectrum3dModule->GetXaxis()->SetTitle("U");
      spectrum3dModule->GetYaxis()->SetTitle("V");
      spectrum3dModule->GetZaxis()->SetTitle("W");
      sname.str("");




      // if(usingRealSimData)
      // {
      //   // GLOBAL SPECTRA but with events confined in one crystal (from sim data)
      //   // Flood histogram
      //   sname << "SIM - Crystal Hit = 1 - Flood Histogram 2D - " << module[iModule][jModule]->GetName();
      //   varModule << "FloodY:FloodX >> " << nameModule;
      //   //       std::cout << nameModule << " ... ";
      //   TH2F *spectrum2dModuleSingleCrystalHit = new TH2F(nameModule,nameModule,histo2DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo2DglobalBins,-moduleLateralSideY,moduleLateralSideY);
      //   tree->Draw(varModule.str().c_str(),SingleCrystalInteraction,"COLZ");
      //   spectrum2dModuleSingleCrystalHit->SetName(nameModule);
      //   spectrum2dModuleSingleCrystalHit->SetTitle(nameModule);
      //   spectrum2dModuleSingleCrystalHit->GetXaxis()->SetTitle("U");
      //   spectrum2dModuleSingleCrystalHit->GetYaxis()->SetTitle("V");
      //   module[iModule][jModule]->SetFloodMap2DSingleCrystalHit(spectrum2dModuleSingleCrystalHit);
      //   sname.str("");
      //   var.str("");
      //
      //
      // }



      // run on MPPCs
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
        for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
        {
          //but proceed only if the MPPC is "on" for modular analysis
          if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForModular())
          {

            //----------------------------------------------------------//
            //  Get channels, variables etc.                            //
            //----------------------------------------------------------//

            //get ch number for this MPPC
            int channel = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
            std::cout << "Generating spectra for MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " ..." << std::endl;


            //------------------------------------------//
            // CREATION OF STRINGS FOR PLOTTING         //
            //------------------------------------------//
            // in order to compute all th plots on the fly, it is necessary to put together the proper strings to pass to the draw commands
            // this has to include also the saturation correction, the pedestal position and the noise
            // furthermore, we need to prepare some collection of elements for the channels that are relevant for computation of energy, u-v coordinates and w coordinates
            // it all starts from the detector array, produced and filled by the InputFile.cc part of the program. Declaration is in ./include/Detector.h
            // this vector of struct holds an entry for each mppc. each entry has data on the mppc. In particular, the digitizer (aka charge) channel number and
            // the timing channel number, a vector with the digitizer channel numbers of all the neighbour channels, and the info on saturation, pedestal and noise
            // The steps are the following (more explainations written just before each step below):

            // 1. Get the index of this MPPC in the detector array                                                                 ---> thisChannelID
            // 2. Get the indexes of the neighbour MPPCs in the detector array                                                     ---> neighbourChID[]
            // 3. Translate indexes into strings, taking also into account saturation (if requested), pedestals and noise          ---> thisChannel
            //                                                                                                                     ---> neighbourChannels
            //                                                                                                                     ---> allChannels
            // 4. Generate the TriggerChannel condition                                                                            ---> TriggerChannel
            // 5. Get the relevant channels for this MPPC                                                                          ---> relevantForUV
            //                                                                                                                     ---> relevantForW
            //                                                                                                                     ---> relevantForEnRes
            // 6. Compose the FloodX and FloodY strings                                                                            ---> FloodX
            //                                                                                                                     ---> FloodY
            // 7. Compose the FloodZ string                                                                                        ---> FloodZ
            // 8. Compose an additional trigger condition, cutting the noise                                                       ---> relevantForNoiseCut
            //                                                                                                                     ---> noiseCut
            //                                                                                                                     ---> summed to CutTrigger


            // 1. Get the index of this MPPC in the detector array
            //    Simply, run on all the entries of detector array and get the one where detector[i].digitizerChannel == mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel()
            //    ---> thisChannelID
            int thisChannelID;
            // get mppc channel ID in detector[]
            for(unsigned int iDet = 0; iDet < detector.size(); iDet++) // run on all detector[] entries
            {
              if(detector[iDet].digitizerChannel == channel) // if the digitizer channel in detector[iDet] is the same as this mppc digitize channel
              {
                thisChannelID = iDet; // get the index of detector[] that identifies it
              }
            }
            // // DEBUG crosscheck
            // std::cout << "channel from GetDigitizerChannel           = " << channel  << std::endl
            //           << "detector[thisChannelID].digitizerChannel   = " << detector[thisChannelID].digitizerChannel  << std::endl
            //           << "thisChannelID                              = " << thisChannelID << std::endl;
            // // end of debug


            // 2. Get the indexes of the neighbour MPPCs in the detector array
            //    Same logic as in 1., run on all entries on detector[] and store in vector neighbourChID
            //    ---> neighbourChID[]
            std::vector<int> neighbourChID;
            for(unsigned int iNeigh = 0; iNeigh < detector[thisChannelID].neighbourChannels.size(); iNeigh++)
            {
              // for each entry, look for the index in detector[] that corresponds to the entry value
              for(unsigned int iDet = 0; iDet < detector.size(); iDet++)
              {
                if(detector[iDet].digitizerChannel == detector[thisChannelID].neighbourChannels[iNeigh])
                {
                  neighbourChID.push_back(iDet); // get the index of detector[] that identifies it
                }
              }
            }

            // //DEBUG
            // std::cout << "neighbourChID[iNeigh]\t      detector[thisChannelID].neighbourChannels[iNeigh]"<< std::endl;
            // for(unsigned int iNeigh = 0; iNeigh < neighbourChID.size(); iNeigh++)
            // {
            //   std::cout << neighbourChID[iNeigh] << "\t\t"
            //             << detector[thisChannelID].neighbourChannels[iNeigh] << std::endl;
            // }
            // std::cout << std::endl;
            // // end of DEBUG


            // 3. Translate indexes into strings, taking also into account saturation (if requested), pedestals and noise
            //    Two info need to be stored together:
            //           a) the channel ID (i.e. the digitizer channel)
            //           b) the complete string that will be passed to the Draw command
            //    So it is necessary to create a multi_channel_t struct, with these two entries (see definition above)
            //    Then, 3 arrays of this multi_channel_t elements are created:
            //    ---> thisChannel        = holds the string and digitizerChannel for this MPPC, so it's not a vector but one struct
            //    ---> neighbourChannels  = vector holding the string and digitizerChannel for each of the neighbour channels to this MPPC
            //    ---> allChannels        = vector holding the string and digitizerChannel for all channels in this MPPC array
            //    The 3 arrays will be used in step 4. and 5.
            multi_channel_t thisChannel;
            std::vector<multi_channel_t> neighbourChannels;
            std::vector<multi_channel_t> allChannels;
            // compose the strings for variables, corrected or not by saturation
            if(correctingSaturation)
            {
              // this channel
              var.str("");
              var <<  "("
              << -detector[thisChannelID].saturation
              << " * TMath::Log(1.0 - ( ((ch"
              << detector[thisChannelID].digitizerChannel
              << " - "
              << detector[thisChannelID].pedestal
              << "))/("
              << detector[thisChannelID].saturation
              << "))))";
              thisChannel.string = var.str();
              thisChannel.detectorIndex = thisChannelID;
              var.str("");

              // neighbour channels
              for(unsigned int iNeigh = 0; iNeigh < neighbourChID.size(); iNeigh++)
              {
                var <<  "("
                << -detector[neighbourChID[iNeigh]].saturation
                << " * TMath::Log(1.0 - ( ((ch"
                << detector[neighbourChID[iNeigh]].digitizerChannel
                << " - "
                << detector[neighbourChID[iNeigh]].pedestal
                << "))/("
                << detector[neighbourChID[iNeigh]].saturation
                << "))))";
                multi_channel_t tempMultiChannel;
                tempMultiChannel.string = var.str();
                tempMultiChannel.detectorIndex = neighbourChID[iNeigh];
                neighbourChannels.push_back(tempMultiChannel);
                var.str("");
              }

              //all channels
              for(unsigned int iAll = 0; iAll < allModuleChID.size(); iAll++)
              {
                var <<  "("
                << -detector[allModuleChID[iAll]].saturation
                << " * TMath::Log(1.0 - ( ((ch"
                << detector[allModuleChID[iAll]].digitizerChannel
                << " - "
                << detector[allModuleChID[iAll]].pedestal
                << "))/("
                << detector[allModuleChID[iAll]].saturation
                << "))))";
                multi_channel_t tempMultiChannel;
                tempMultiChannel.string = var.str();
                tempMultiChannel.detectorIndex = allModuleChID[iAll];
                allChannels.push_back(tempMultiChannel);
                var.str("");
              }

              // std::cout << thisChannel << std::endl;
              // std::cout << std::endl;
              //
              // for(unsigned int iNeigh = 0; iNeigh < neighbourChannels.size(); iNeigh++)
              // {
              //   std::cout << neighbourChannels[iNeigh] << std::endl;
              // }
              // std::cout << std::endl;
              //
              // for(unsigned int iAll = 0; iAll < allChannels.size(); iAll++)
              // {
              //   std::cout << allChannels[iAll] << std::endl;
              // }
              // std::cout << std::endl;
            }
            else
            {
              // this channel
              var.str("");
              var <<  "(ch" <<  detector[thisChannelID].digitizerChannel << "-" << detector[thisChannelID].pedestal << ")";
              thisChannel.string = var.str();
              thisChannel.detectorIndex = thisChannelID;
              var.str("");

              // neighbour channels
              for(unsigned int iNeigh = 0; iNeigh < neighbourChID.size(); iNeigh++)
              {
                var <<  "(ch" <<  detector[neighbourChID[iNeigh]].digitizerChannel << "-" << detector[neighbourChID[iNeigh]].pedestal << ")";
                multi_channel_t tempMultiChannel;
                tempMultiChannel.string = var.str();
                tempMultiChannel.detectorIndex = neighbourChID[iNeigh];
                neighbourChannels.push_back(tempMultiChannel);
                var.str("");
              }

              //all channels
              for(unsigned int iAll = 0; iAll < allModuleChID.size(); iAll++)
              {
                var <<  "(ch" <<  detector[allModuleChID[iAll]].digitizerChannel << "-" << detector[allModuleChID[iAll]].pedestal << ")";
                multi_channel_t tempMultiChannel;
                tempMultiChannel.string = var.str();
                tempMultiChannel.detectorIndex = allModuleChID[iAll];
                allChannels.push_back(tempMultiChannel);
                var.str("");
              }

              // std::cout << thisChannel << std::endl;
              //
              // for(unsigned int iNeigh = 0; iNeigh < neighbourChannels.size(); iNeigh++)
              // {
              //   std::cout << neighbourChannels[iNeigh] << " ";
              // }
              // std::cout << std::endl;
              //
              // for(unsigned int iAll = 0; iAll < allChannels.size(); iAll++)
              // {
              //   std::cout << allChannels[iAll] << " ";
              // }
              // std::cout << std::endl;
            }




            // 4. Generate the TriggerChannel condition
            //    the idea is:        TriggerChannel     = "max(ch0,max(ch1,...max(chN))) == chI", for I == mppc channel and 0...N are all the channels in this module
            //    in this case the MPPC channels involved are all the channels, so the allChannels array
            std::stringstream TriggerChannel;
            int parenthesisCounter = 0;
            for(unsigned int iDet = 0; iDet < allChannels.size()-1; iDet++)
            {
              TriggerChannel << "max(" << allChannels[iDet].string << ",";
              parenthesisCounter++;
            }
            TriggerChannel << allChannels[allChannels.size()-1].string;
            for(unsigned int iDet = 0; iDet < parenthesisCounter; iDet++)
            {
              TriggerChannel << ")";
            }

            //DEBUG
            // std::cout << TriggerChannel.str().c_str() << std::endl;

            cut.str("");
            cut << TriggerChannel.str().c_str() << " == " << thisChannel.string.c_str();
            TCut CutTrigger = cut.str().c_str();



            //----------------------------------------------------------//
            // Draw plots                                               //
            //----------------------------------------------------------//

            // simple Raw Spectrum for this MPPC
            // it will be stored in the MPPC class and then retreived later for saving in the output file
            sname << "Raw Spectrum - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " - Module " << module[iModule][jModule]->GetName();
            var << thisChannel.string.c_str() << " >> " << sname.str();
            TH1F* spectrumRaw = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,0,histo1Dmax);
            tree->Draw(var.str().c_str(),BasicCut);
            spectrumRaw->GetXaxis()->SetTitle("ADC Channels");
            spectrumRaw->GetYaxis()->SetTitle("N");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetRawSpectrum(spectrumRaw);
            var.str("");
            sname.str("");

            //trigger selected spectrum
            sname << "Trigger Spectrum - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            var << thisChannel.string.c_str() << " >> " << sname.str().c_str();
            TH1F* spectrumTrigger = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,0,histo1Dmax);
            tree->Draw(var.str().c_str(),CutTrigger+BasicCut);
            spectrumTrigger->GetXaxis()->SetTitle("ADC Channels");
            spectrumTrigger->GetYaxis()->SetTitle("N");
            // mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTriggerSpectrum(spectrumTrigger);
            var.str("");
            sname.str("");


            //set a very broad cut on the trigger spectrum (single channel, not sum) to get rid of low energy events
            std::stringstream broadCutstream;
            if(!backgroundRun)
            {
              if(thresholdKev > 0.0)
              {
                TSpectrum *sTrigger;
                sTrigger = new TSpectrum(20);
                Int_t TriggerCrystalPeaksN    = sTrigger->Search(spectrumTrigger,2,"goff",0.2);
                Double_t *TriggerCrystalPeaks  = sTrigger->GetPositionX();
                // Float_t *TriggerCrystalPeaksY = sTrigger->GetPositionY();
                //delete s;
                float TriggermaxPeak = 0.0;
                int TriggerpeakID    = 0;
                for (int TriggerpeakCounter = 0 ; TriggerpeakCounter < TriggerCrystalPeaksN ; TriggerpeakCounter++ )
                {
                  if(TriggerCrystalPeaks[TriggerpeakCounter] > TriggermaxPeak)
                  {
                    TriggermaxPeak = TriggerCrystalPeaks[TriggerpeakCounter];
                    TriggerpeakID = TriggerpeakCounter;
                  }
                }
                // now if 511KeV or 662KeV correspond to TriggerCrystalPeaks[TriggerpeakID], it means that a broad cut, energy > 200-250KeV is approximately that divided by 2.5 (assuming 0 is 0 and scale is linear)

                broadCutstream << thisChannel.string.c_str() << ">" << (TriggerCrystalPeaks[TriggerpeakID] / ((1000*sourceMeV)/thresholdKev) );
              }
              else
              {
                broadCutstream << "";
              }

            }
            else //otherwise, for backgroundRun set the broadcut to the value chosen by the user (or default to 0 adc channels)
            {
              if(userBroadCut > 0)
              {
                broadCutstream << thisChannel.string.c_str() << ">" << userBroadCut;
              }
              else
              {
                broadCutstream << "";
              }

            }
            TCut broadCut = broadCutstream.str().c_str();
            // to make things easier, we add this directly to the CutTrigger. FIXME Not a clean solution, but ehi...
            CutTrigger += broadCut;
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTriggerSpectrum(spectrumTrigger);
            var.str("");

            if( (thresholdKev > 0.0) || ( (backgroundRun) && (userBroadCut > 0))  )
            {
              //prepare an highlighted plot to show the broad cut
              sname.str("");
              sname << "Trigger Spectrum Hg - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
              var << thisChannel.string.c_str() << " >> " << sname.str().c_str();
              TH1F* spectrumTriggerHighlighted = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,0,histo1Dmax);
              tree->Draw(var.str().c_str(),CutTrigger+BasicCut);
              spectrumTriggerHighlighted->GetXaxis()->SetTitle("ADC Channels");
              spectrumTriggerHighlighted->GetYaxis()->SetTitle("N");
              mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTriggerSpectrumHighlighted(spectrumTriggerHighlighted);
              var.str("");
              sname.str("");
            }





            // //trigger spectrum of this mppc, expressed in charge, WITHOUT any consideration on where the crystal of interaction.
            // //useful for getting the peaks of lu-176 emitted from OTHER crystals, that interact in crystal(s) coupled to this detector.
            // //it's the same as the trigger spectrum done before, but expressed in charge
            // sname << "Trigger Charge Spectrum - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            // var << thisChannel.string.c_str() << "*( "<< chargeBinningADC <<  ") >> " << sname.str().c_str();
            // TH1F* spectrumTriggerCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histoSingleChargeBin,0,histoSingleChargeMax);
            // tree->Draw(var.str().c_str(),BasicCut+CutTrigger);
            // spectrumTriggerCharge->GetXaxis()->SetTitle("Charge [C]");
            // spectrumTriggerCharge->GetYaxis()->SetTitle("N");
            // sname.str("");
            // var.str("");
            // mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTriggerChargeSpectrum(spectrumTriggerCharge);
            //
            // //even more basic, the raw spectrum of this channel, expressed in terms of charge
            // sname << "Raw Charge Spectrum - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            // var << thisChannel.string.c_str() << "*( "<< chargeBinningADC <<  ") >> " << sname.str().c_str();
            // TH1F* spectrumRawCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histoSingleChargeBin,0,histoSingleChargeMax);
            // tree->Draw(var.str().c_str(),BasicCut);
            // spectrumRawCharge->GetXaxis()->SetTitle("Charge [C]");
            // spectrumRawCharge->GetYaxis()->SetTitle("N");
            // sname.str("");
            // var.str("");
            // mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetRawChargeSpectrum(spectrumRawCharge);


            // 5. Get the relevant channels for this MPPC
            //    Relevant channels are split in 3
            //    5.1. relevant for u-v coordinates
            //    5.2. relevant for w coordinate
            //    5.3. relevant for energy resolution (or better for sum spectrum)
            //    The 3 conditions are controlled by the corresponding config file keys
            //    5.1. usingAllChannels  -> 0 = just neighbour channels
            //                         -> 1 = all channels in the module
            //    5.2. wAllChannels      -> 0 = just neighbour channels
            //                         -> 1 = all channels in the module
            //    5.3. usingAllChannelsForEnergySpectra  -> 0 = just neighbour channels
            //                                         -> 1 = all channels in the module
            //    Remember that the string for the channels, already including the saturation correction if any, have been created in step 3.
            //    multi_channel_t thisChannel;
            //    std::vector<multi_channel_t> neighbourChannels;
            //    std::vector<multi_channel_t> allChannels;


            // 5.1. relevant for u-v coordinates
            std::vector<multi_channel_t> relevantForUV;
            if(usingAllChannels) // all channels used for u-v
            {
              for(unsigned int iRel = 0 ; iRel < allChannels.size(); iRel++)
              {
                relevantForUV.push_back(allChannels[iRel]);
              }
            }
            else // trigger channel (i.e. this channel) plus neighbour channels used for u-v
            {
              relevantForUV.push_back(thisChannel);
              for(unsigned int iRel = 0 ; iRel < neighbourChannels.size(); iRel++)
              {
                relevantForUV.push_back(neighbourChannels[iRel]);
              }
            }


            // 5.2. relevant for w coordinate
            std::vector<multi_channel_t> relevantForW;
            if(wAllChannels)// all channels used for w
            {
              for(unsigned int iRel = 0 ; iRel < allChannels.size(); iRel++)
              {
                relevantForW.push_back(allChannels[iRel]);
              }
            }
            else // trigger channel (i.e. this channel) plus neighbour channels used for w
            {
              relevantForW.push_back(thisChannel);
              for(unsigned int iRel = 0 ; iRel < neighbourChannels.size(); iRel++)
              {
                relevantForW.push_back(neighbourChannels[iRel]);
              }
            }


            // 5.3. relevant for energy resolution
            std::vector<multi_channel_t> relevantForEnRes;
            if(usingAllChannelsForEnergySpectra) // all channels used for energy spectra
            {
              for(unsigned int iRel = 0 ; iRel < allChannels.size(); iRel++)
              {
                relevantForEnRes.push_back(allChannels[iRel]);
              }
            }
            else // trigger channel (i.e. this channel) plus neighbour channels used for energy spectra
            {
              relevantForEnRes.push_back(thisChannel);
              for(unsigned int iRel = 0 ; iRel < neighbourChannels.size(); iRel++)
              {
                relevantForEnRes.push_back(neighbourChannels[iRel]);
              }
            }


            // 6. Compose the FloodX and FloodY strings
            //    Standard Anger logic of charges (corrected by saturation and pedestal) with mppc position, for the channels defined as relevant for u-v
            std::stringstream FloodX;
            std::stringstream FloodY;
            FloodX << "(" ;
            FloodY << "(" ;
            for(unsigned int iRel = 0 ; iRel < relevantForUV.size() - 1; iRel++)
            {
              FloodX << "("
                     << relevantForUV[iRel].string << "*(" << detector[relevantForUV[iRel].detectorIndex].xPosition
                     << ")) + ";
              FloodY << "("
                     << relevantForUV[iRel].string << "*(" << detector[relevantForUV[iRel].detectorIndex].yPosition
                     << ")) + ";
            }
            FloodX << "("
                   << relevantForUV[relevantForUV.size() - 1].string << "*(" << detector[relevantForUV[relevantForUV.size() - 1].detectorIndex].xPosition
                   << "))";
            FloodX << ") / (" ;
            FloodY << "("
                   << relevantForUV[relevantForUV.size() - 1].string << "*(" << detector[relevantForUV[relevantForUV.size() - 1].detectorIndex].yPosition
                   << "))";
            FloodY << ") / (" ;
            for(unsigned int iRel = 0 ; iRel < relevantForUV.size() - 1; iRel++)
            {
              FloodX << "("
                     << relevantForUV[iRel].string
                     << ") + ";
              FloodY << "("
                     << relevantForUV[iRel].string
                     << ") + ";
            }
            FloodX << "("
                   << relevantForUV[relevantForUV.size() - 1].string
                   << ")";
            FloodX << ")";
            FloodY << "("
                   << relevantForUV[relevantForUV.size() - 1].string
                   << ")";
            FloodY << ")";


            // 7. Compose the FloodZ string
            //    Ratio of charge of this mppc versus charges of the relevant mppcs for w computation. all charges are corrected by saturation and pedestal
            std::stringstream FloodZ;
            FloodZ << "( "<< thisChannel.string << ") / (";
            for(unsigned int iRel = 0 ; iRel < relevantForW.size() - 1; iRel++)
            {
              FloodZ << "("
                     << relevantForW[iRel].string
                     << ") + ";
            }
            FloodZ << "("
                   << relevantForUV[relevantForW.size() - 1].string
                   << ")";
            FloodZ << ")";


            // 8. Compose an additional trigger condition, cutting the noise
            //    the idea is to ignore events when the charge recorded by one of the detectors involved is == 0
            //    ((ch12 > 0)+ (ch11>0) + (ch0>0)+ (ch4>0) + (ch7>0) +(ch10>0) +(ch13>0) +(ch14>0)+ (ch15>0)) > 8
            //    the channels involved are defined as the biggest collection between relevantForUV and relevantForW
            std::vector<multi_channel_t> relevantForNoiseCut;
            if(relevantForUV.size() > relevantForW.size())
            {
              relevantForNoiseCut = relevantForUV;
            }
            else
            {
              relevantForNoiseCut = relevantForW;
            }

            std::stringstream noiseCut;
            noiseCut << "(";
            for(unsigned int iRel = 0 ; iRel < relevantForNoiseCut.size() -1; iRel++)
            {
              noiseCut << "("
                       << relevantForNoiseCut[iRel].string << " > "
                       << noiseSigmas * detector[relevantForNoiseCut[iRel].detectorIndex].noise
                       << ") + ";
            }
            noiseCut << "("
                       << relevantForNoiseCut[relevantForNoiseCut.size() -1].string << " > "
                       << noiseSigmas * detector[relevantForNoiseCut[relevantForNoiseCut.size() -1].detectorIndex].noise
                       << ")";
            noiseCut << ") > " << relevantForNoiseCut.size() -1;
            // std::cout << noiseCut.str() << std::endl;
            // std::cout << std::endl;

            TCut CutNoise = noiseCut.str().c_str();
            // to make things easier, we add this directly to the CutTrigger. FIXME Not a clean solution, but ehi...
            CutTrigger += CutNoise;


            //standard 2d plot
            sname << "Flood Histogram 2D - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            TH2F* spectrum2dMPPC = new TH2F(sname.str().c_str(),sname.str().c_str(),histo2DchannelBin,-moduleLateralSideX,moduleLateralSideX,histo2DchannelBin,-moduleLateralSideY,moduleLateralSideY);
            var << FloodY.str() << " : " << FloodX.str() << " >> " << sname.str().c_str();
            tree->Draw(var.str().c_str(),BasicCut+CutTrigger,"COLZ");
            spectrum2dMPPC->GetXaxis()->SetTitle("U");
            spectrum2dMPPC->GetYaxis()->SetTitle("V");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetFloodMap2D(spectrum2dMPPC);
            sname.str("");
            var.str("");

            spectrum2dModule->Add(spectrum2dMPPC);

            // 3D spectrum for this mppc
            // little trick to try and use less bins
            // take the mean x and y and their sigma from previous 2dplot, define the limit of this 3dplot in x and y accordingly
            double minX3Dplot = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetMean(1) - 3.0*mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetRMS(1);
            double maxX3Dplot = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetMean(1) + 3.0*mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetRMS(1);
            double minY3Dplot = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetMean(2) - 3.0*mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetRMS(2);
            double maxY3Dplot = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetMean(2) + 3.0*mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetRMS(2);

            //DEBUG
            // 	  std::cout << "###### Main Program " << std::endl;
            // 	  std::cout << minX3Dplot << " " << maxX3Dplot << " " << minY3Dplot << " "<< maxY3Dplot << std::endl;
            // 	  std::cout << "----------------- " << std::endl;
            sname << "Flood Histogram 3D - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            var << FloodZ.str() << " : "
                << FloodY.str() << " : "
                << FloodX.str() << " >> " << sname.str().c_str();

            int histo3DchannelBin =  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetHisto3DchannelBin();
            TH3I* spectrum3dMPPC = new TH3I(sname.str().c_str(),sname.str().c_str(),histo3DchannelBin,minX3Dplot,maxX3Dplot,histo3DchannelBin,minY3Dplot,maxY3Dplot,histo3DchannelBin,histo3Dmin,histo3Dmax);//FIXME temp
            // if(usingTaggingBench | calcDoiResWithCalibration)
            // {
            tree->Draw(var.str().c_str(),BasicCut+CutTrigger);
            // }
            // else
            // {
            // tree->Draw(var.str().c_str(),BasicCut+CutTrigger);
            // }
            spectrum3dMPPC->GetXaxis()->SetTitle("U");
            spectrum3dMPPC->GetYaxis()->SetTitle("V");
            spectrum3dMPPC->GetZaxis()->SetTitle("W");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetFloodMap3D(spectrum3dMPPC);
            sname.str("");
            var.str("");

            // another 3d plot same data as before, but with histo limits same as global plot,
            // to make some plot possible
            sname << "spectrum3dMPPC_large - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            var << FloodZ.str() << " : "
                << FloodY.str() << " : "
                << FloodX.str() << " >> " << sname.str().c_str();
            TH3I* spectrum3dMPPC_large = new TH3I(sname.str().c_str(),sname.str().c_str(),histo3DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo3DglobalBins,-moduleLateralSideY,moduleLateralSideY,histo3DglobalBins,histo3Dmin,histo3Dmax);
            tree->Draw(var.str().c_str(),BasicCut+CutTrigger);
            spectrum3dMPPC_large->GetXaxis()->SetTitle("U");
            spectrum3dMPPC_large->GetYaxis()->SetTitle("V");
            spectrum3dMPPC_large->GetZaxis()->SetTitle("W");
            sname.str("");
            spectrum3dModule->Add(spectrum3dMPPC_large);
            var.str("");



            //----------------------------------------------------------//
            // Separate Crystals                                        //
            //----------------------------------------------------------//
            // automatic crystal finder on the mppc
            // it runs always, unless the user has set onlyuserinput //FIXME onlyuserinput is not used anymore!
            std::vector<TString> FloodXYZ;
            FloodXYZ.push_back(FloodX.str());
            FloodXYZ.push_back(FloodY.str());
            FloodXYZ.push_back(FloodZ.str());

            TCutG**** cutg; // prepare the graphical cuts
            //const int numbOfCrystals = 4;
            cutg = new TCutG***[2]; // two planes of cuts, their intersection will create a 3d cut
            int right_ncrystalsx;
            if(usingTaggingBench)
            { //with tagging crystal setup, use only one row of crystals
              right_ncrystalsx =1;
            }
            else
            {
              right_ncrystalsx =ncrystalsx;
            }
            for(int iCut =0 ; iCut < 2 ; iCut++)
            {
              cutg[iCut] = new TCutG**[right_ncrystalsx];
              for(int iCry = 0; iCry < right_ncrystalsx ; iCry++)
              {
                cutg[iCut][iCry] = new TCutG*[ncrystalsy];
              }
            }
            bool found = false;
            if(usingTaggingBench) //if it's a tagging bench run, check first if the mppc is on for DOI bench measurement
            {
              if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForDoi())
              {
                found = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->FindCrystalCuts(cutg,1,ncrystalsy,FloodXYZ);
              }
            }
            else
            {
              found = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->FindCrystalCuts(cutg,ncrystalsx,ncrystalsy,FloodXYZ);
            }

            // run on all the possible crystals (i.e. all the crystals coupled to this mppc)
            for(int iCry = 0; iCry < right_ncrystalsx ; iCry++)
            {
              for(int jCry = 0; jCry < ncrystalsy ; jCry++)
              {
                if(crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)]->GetIsOnForModular())
                {
                  if(found)
                  {
                    // get a pointer to this crystal
                    Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
                    CurrentCrystal->SetCrystalOn(true);
                    //store the cutg in the crystal
                    CurrentCrystal->SetZXCut(cutg[0][iCry][jCry]);
                    CurrentCrystal->SetZYCut(cutg[1][iCry][jCry]);

                    //get charge ch numbers used for w
                    std::vector<int> channelsNumRelevantForW;
                    for(unsigned int iRel = 0; iRel < relevantForW.size(); iRel++)
                    {
                      channelsNumRelevantForW.push_back(relevantForW[iRel].detectorIndex);
                    }

                    CurrentCrystal->SetRelevantForW(channelsNumRelevantForW);

                    // create "sum channels" string, a string to have the variable "sum of all channels" to be used later
                    // it can be either all channels selected as input
                    // or just trigger channel plus neighbourChannels
                    // already found this strings, so just need to copy

                    // first get the input digitizer channels
                    // std::vector<std::string> digitizer_f;
                    // config.split( digitizer_f, digitizer_s, "," );
                    // std::vector<int> digitizer;
                    // for(unsigned int i = 0 ; i < digitizer_f.size() ; i++)
                    // {
                    //   config.trim(digitizer_f[i]);
                    //   digitizer.push_back(atoi(digitizer_f[i].c_str()));
                    // }
                    // create the string
                    std::stringstream sSumChannels;
                    std::string SumChannels;
                    std::vector<int> neighbours = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetNeighbours(); //get neighbours anyway, they will be useful later

                    sSumChannels << "(";
                    for(unsigned int iRel = 0 ; iRel < relevantForEnRes.size() - 1 ; iRel++)
                    {
                      sSumChannels << "(" <<relevantForEnRes[iRel].string << ") + ";
                    }
                    sSumChannels << "(" <<relevantForEnRes[relevantForEnRes.size() - 1].string << ") ) ";

                    SumChannels = sSumChannels.str();

                    // if(usingAllChannelsForEnergySpectra) // add all channels in the digitizer list
                    // {
                    //   sSumChannels << "ch" <<  digitizer[0];
                    //   for(unsigned int iDigi = 1 ; iDigi < digitizer.size() ; iDigi++)
                    //   {
                    //     sSumChannels << "+ch" << digitizer[iDigi];
                    //   }
                    //   SumChannels = sSumChannels.str();
                    //   sSumChannels.str("");
                    // }
                    // else // add only this channel and the neighbourChannels
                    // {
                    //   sSumChannels << "ch" << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
                    //
                    //   for(unsigned int iNeighbour = 0 ; iNeighbour < neighbours.size() ; iNeighbour++)
                    //   {
                    //     sSumChannels << "+ch" << neighbours[iNeighbour];
                    //   }
                    //   SumChannels = sSumChannels.str();
                    //   sSumChannels.str("");
                    // }


                    // std::cout << SumChannels << std::endl;
                    //----------------------------------------------------------//





                    std::cout << "Generating spectra for crystal " << CurrentCrystal->GetID() << " ..." << std::endl;

                    // a cut for events falling in this crystal. This basically means
                    // 1. "Physical" constrains are satisfied (u and v are not beyond the broad physical x-y limits of the array, w is between 0 and 1). This is set at module level - FIXME not implemented in this version! (is it useful?)
                    // 2. If a tagging crystal is used, events have charge in the channel of tagging crystal in the photopeak of the tagging crystal. This is set at module level
                    // 3. The MPPC where the crystal is has the maximum signal in the event ("trigger" channel).
                    // 4. The "broad" cut on energy is applied, which means that the user decided to discard events that in the "trigger" spectrum of this mppc are "too low". This is set at mppc level
                    // 5. None of the channels involved is below noiseSigmas*noise for that channel
                    // 6. Events are in the density area of this crystal, defined by the two tcutg
                    //
                    // BasicCut = 1. and 2.
                    // CutTrigger = 3. and 4. and 5.
                    // CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName() = 6.

                    TCut CrystalCut = BasicCut + CutTrigger + CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName();

                    // a 3d historgram for this crystal, mainly to check the 3d cut
                    sname << "Flood Histogram 3D - Crystal " << CurrentCrystal->GetID();
                    var << FloodZ.str() << " : "
                        << FloodY.str() << " : "
                        << FloodX.str() << " >> " << sname.str().c_str();
                    // std::cout << var.str() << std::endl;
                    // NB histo3DchannelBin is taken from the MPPC element earlier in the code
                    TH3I* spectrum3dCrystal = new TH3I(sname.str().c_str(),sname.str().c_str(),histo3DchannelBin,minX3Dplot,maxX3Dplot,histo3DchannelBin,minY3Dplot,maxY3Dplot,histo3DchannelBin,histo3Dmin,histo3Dmax);
                    tree->Draw(var.str().c_str(),CrystalCut);
                    spectrum3dCrystal->GetXaxis()->SetTitle("U");
                    spectrum3dCrystal->GetYaxis()->SetTitle("V");
                    spectrum3dCrystal->GetZaxis()->SetTitle("W");
                    CurrentCrystal->SetFloodMap3D(spectrum3dCrystal);
                    sname.str("");
                    var.str("");

                    // a 2d historgram for this crystal, mainly to check the 3d cut
                    sname << "Flood Histogram 2D - Crystal " << CurrentCrystal->GetID();
                    var << FloodY.str() << " : "
                        << FloodX.str() << " >> " << sname.str().c_str();
                    TH2F* spectrum2dCrystal = new TH2F(sname.str().c_str(),sname.str().c_str(),histo2DchannelBin,-moduleLateralSideX,moduleLateralSideX,histo2DchannelBin,-moduleLateralSideY,moduleLateralSideY);
                    tree->Draw(var.str().c_str(),CrystalCut);
                    spectrum2dCrystal->GetXaxis()->SetTitle("U");
                    spectrum2dCrystal->GetYaxis()->SetTitle("V");
                    CurrentCrystal->SetFloodMap2D(spectrum2dCrystal);
                    sname.str("");
                    var.str("");


                    if(saturationRun) //only if this is a saturation analysis run
                    {
                      //sum spectrum for events localized in this crystal, but expressed in Q
                      sname << "Sum Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      // var << SumChannels << " >> " << sname.str();
                      var << "(" << SumChannels << ")*( "<< chargeBinningADC <<  ") >> " << sname.str();
                      TH1F* spectrumSumCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histoSumChargeBin,0,histoSumChargeMax);
                      tree->Draw(var.str().c_str(),CrystalCut);
                      spectrumSumCharge->GetXaxis()->SetTitle("Charge [C]");
                      spectrumSumCharge->GetYaxis()->SetTitle("N");
                      sname.str("");
                      var.str("");
                      CurrentCrystal->SetSumChargeSpectrum(spectrumSumCharge);
                      // CurrentCrystal->SetSaturationFits(gaussFitSaturation);


                      //sum of all channels minus the trigger channel, for events confined in this crystal, expressed in Q
                      //this should be almost not saturated. useful for sources, to compare with background spectrum and get an energy measurement of relevant parts (to use in the
                      // identification of relevant part of bg spectrum in single channels spectrum of bg)
                      sname << "All-Minus-Trigger Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      // var << SumChannels << " >> " << sname.str();
                      var << "(" << SumChannels << " - "
                          << thisChannel.string
                          << ")*( "
                          << chargeBinningADC <<  ") >> " << sname.str();
                      TH1F* spectrumAMTCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histoSumChargeBin,0,histoSumChargeMax);
                      tree->Draw(var.str().c_str(),CrystalCut);
                      spectrumAMTCharge->GetXaxis()->SetTitle("Charge [C]");
                      spectrumAMTCharge->GetYaxis()->SetTitle("N");
                      sname.str("");
                      var.str("");
                      CurrentCrystal->SetAMTChargeSpectrum(spectrumAMTCharge);
                      // CurrentCrystal->SetSaturationFits(gaussFitSaturation);


                      //single MPPC spectrum for this crystal. for the saturation Run the important quantity is the q_max, i.e. th maximum charge that can be
                      //"seen" by the SiPM, corresponding to a situation where all the pixels are firing. So these plots are generate not in ADCch but in charge, using the chargeBinningADC
                      sname << "Single Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      var << thisChannel.string
                          << "*( "
                          << chargeBinningADC
                          <<  ") >> " << sname.str();
                      TH1F* spectrumSingleCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histoSingleChargeBin,0,histoSingleChargeMax);
                      tree->Draw(var.str().c_str(),CrystalCut);
                      spectrumSingleCharge->GetXaxis()->SetTitle("Charge [C]");
                      spectrumSingleCharge->GetYaxis()->SetTitle("N");
                      sname.str("");
                      var.str("");


                      //spectrum of charge in this mppc, but when the trigger is NOT itself nor a neighbouring channel
                      sname << "Not-Neighbour Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      // std::vector<int> neighbours = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetNeighbours();
                      TCut NotNeighbours;
                      std::stringstream sneigh;
                      sneigh << TriggerChannel.str() <<" != " << thisChannel.string ;

                      NotNeighbours += sneigh.str().c_str();
                      for(unsigned int iRel = 0 ; iRel < neighbourChannels.size(); iRel++)
                      {
                        std::stringstream sneigh2;
                        sneigh2.str("");
                        sneigh2 << TriggerChannel.str() <<" != " <<   neighbourChannels[iRel].string;
                        NotNeighbours += sneigh2.str().c_str();
                      }
                      var << thisChannel.string
                          << "*( "
                          << chargeBinningADC
                          <<  ") >> " << sname.str();
                      TH1F* NotNeighboursSingleCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histoSingleChargeBin,0,histoSingleChargeMax);
                      // std::cout << NotNeighbours << std::endl;
                      tree->Draw(var.str().c_str(),BasicCut+NotNeighbours);
                      NotNeighboursSingleCharge->GetXaxis()->SetTitle("Charge [C]");
                      NotNeighboursSingleCharge->GetYaxis()->SetTitle("N");

                      sname.str("");
                      var.str("");
                      CurrentCrystal->SetNotNeighboursSingleCharge(NotNeighboursSingleCharge);




                      //look for peaks in the areas selected by the user
                      //here there could be more peaks the user is interested in (for example, in na22)
                      //and they could be anywhere, so it's pointless to implement a very complicated algorithm
                      //to find all possible peaks...
                      std::vector<TF1*> gaussFitSaturation;
                      TH1F* spectrumWhereToSearch;  //the spectrum where to search peaks is different if the run is backgroun saturation or not.

                      if(backgroundSaturationRun)
                      {
                        spectrumWhereToSearch = (TH1F*) NotNeighboursSingleCharge->Clone();
                      }
                      else
                      {
                        spectrumWhereToSearch = (TH1F*) spectrumSingleCharge->Clone();
                      }
                      // in background saturation runs, in fact, the spectrum where to search is the trigger spectrum without any crystal consideration, where the 202 Kev and 307 Kev peaks
                      // are much more clear. otherwise the single channel spectrum cut on the events in the crystal is used
                      if(performSaturationPeakSearch)
                      {
                        spectrumWhereToSearch->SetName("Investigated Spectrum");
                        float maxPeakHight = 0;
                        for(unsigned int iSaturation = 0 ; iSaturation < saturationPeak.size(); iSaturation++)
                        {
                          sname << "Peak " << saturationPeak[iSaturation].energy << " KeV - Crystal " << CurrentCrystal->GetID() << " - MPPC " <<  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();

                          spectrumWhereToSearch->GetXaxis()->SetRangeUser(saturationPeak[iSaturation].peakMin,saturationPeak[iSaturation].peakMax);

                          TSpectrum *s;
                          s = new TSpectrum(20);
                          // 		Input[i].SumSpectraCanvas->cd(j+1);
                          Int_t CrystalPeaksN = s->Search(spectrumWhereToSearch,1,"nobackground goff",0.3);
                          Double_t *CrystalPeaks = s->GetPositionX();
                          Double_t *CrystalPeaksY = s->GetPositionY();
                          // float saturationPeakFraction
                          //delete s;
                          // float distPeak = INFINITY;
                          float maxPeak = 0.0;
                          int peakID = 0;
                          for (int peakCounter = 0 ; peakCounter < CrystalPeaksN ; peakCounter++ )
                          {
                            // if( fabs(CrystalPeaks[peakCounter] - 0.5*(saturationPeak[iSaturation].peakMin+saturationPeak[iSaturation].peakMax)) < distPeak)//take closest peak to the center of search range selected by the user
                            if(CrystalPeaksY[peakCounter] > maxPeak)
                            {
                              // distPeak = CrystalPeaks[peakCounter];
                              maxPeak = CrystalPeaksY[peakCounter];
                              peakID = peakCounter;
                            }
                          }


                          TF1 *satGauss = new TF1(sname.str().c_str(),  "gaus",CrystalPeaks[peakID]-CrystalPeaks[peakID]*saturationPeakFractionLow,CrystalPeaks[peakID]+CrystalPeaks[peakID]*saturationPeakFractionHigh);
                          satGauss->SetParameter(1,CrystalPeaks[peakID]);
                          satGauss->SetParameter(0,CrystalPeaksY[peakID]);
                          spectrumWhereToSearch->Fit(sname.str().c_str(),"QN","",CrystalPeaks[peakID]-CrystalPeaks[peakID]*saturationPeakFractionLow,CrystalPeaks[peakID]+CrystalPeaks[peakID]*saturationPeakFractionHigh);
                          gaussFitSaturation.push_back(satGauss);
                          saturationFile << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << "\t"
                                         << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetI() << "\t"
                                         << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetJ() << "\t"
                                         << CurrentCrystal->GetID() << "\t"
                                         << CurrentCrystal->GetI() << "\t"
                                         << CurrentCrystal->GetJ() << "\t"
                                         << saturationPeak[iSaturation].energy << "\t"
                                         << satGauss->GetParameter(1) << "\t"
                                         << satGauss->GetParError(1)
                                         << std::endl;
                          spectrumWhereToSearch->GetXaxis()->SetRangeUser(0,histoSingleChargeMax);
                          if(CrystalPeaksY[peakID] > maxPeakHight )
                          {
                            maxPeakHight = CrystalPeaksY[peakID];
                          }

                          sname.str("");
                        }
                        spectrumWhereToSearch->GetYaxis()->SetRangeUser(0,maxPeakHight*2.0);
                        CurrentCrystal->SetSaturationFits(gaussFitSaturation);
                      }
                      CurrentCrystal->SetInvestigatedSpectrum(spectrumWhereToSearch);
                      CurrentCrystal->SetSingleChargeSpectrum(spectrumSingleCharge);

                    }
                    else
                    {

                      //-------------------------------------------------------------------------
                      //standard sum spectrum with cut on crystal events, xyz and trigger channel
                      //------------------------------------------------------------------------
                      //draw charge spectrum
                      sname << "Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      var << SumChannels << " >> " << sname.str();
                      TH1F* spectrumCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut);
                      spectrumCharge->GetXaxis()->SetTitle("ADC Channels");
                      spectrumCharge->GetYaxis()->SetTitle("N");
                      sname.str("");
                      var.str("");
                      if(lightYieldComputation)
                      {
                        //SumSpectrum in Ph/MeV -- CAREFUL this is not as accurate as measuring LY on PMTs
                        sname << "Light Yield @ "<< sourceMeV <<  " MeV - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                        var << "(" << SumChannels << ")* " << chargeBinningADC << "/(" << gainMPPC <<"*1.6e-19* " << qe << "*"<<  sourceMeV << ") >> " << sname.str();
                        TH1F* spectrumLY = new TH1F(sname.str().c_str(),sname.str().c_str(),histoLYbins,1,histoLYmax);
                        tree->Draw(var.str().c_str(),CrystalCut);
                        spectrumLY->GetXaxis()->SetTitle("[Ph/MeV]");
                        spectrumLY->GetYaxis()->SetTitle("N");
                        CurrentCrystal->SetLYSpectrum(spectrumLY);
                        sname.str("");
                        var.str("");

                        TSpectrum *s;
                        s = new TSpectrum(20);
                        // 		Input[i].SumSpectraCanvas->cd(j+1);
                        Int_t CrystalPeaksN = s->Search(spectrumLY,1,"",0.5);
                        Double_t *CrystalPeaks = s->GetPositionX();
                        Double_t *CrystalPeaksY = s->GetPositionY();
                        //delete s;
                        float maxPeak = 0.0;
                        int peakID = 0;
                        for (int peakCounter = 0 ; peakCounter < CrystalPeaksN ; peakCounter++ )
                        {

                          if(CrystalPeaks[peakCounter] > maxPeak)
                          {
                            maxPeak = CrystalPeaks[peakCounter];
                            peakID = peakCounter;
                          }
                        }
                        //fit the spectra - TODO use the gaussian plus fermi?
                        if (energyResolution == 0)
                        {
                          if (correctingSaturation)
                            energyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
                          else
                            energyResolution = ENERGY_RESOLUTION;
                        }
                        float par0 = CrystalPeaksY[peakID];
                        float par1 = CrystalPeaks[peakID];
                        float par2 = (CrystalPeaks[peakID]*energyResolution)/2.35;
                        float fitmin = par1-1.0*par2;
                        float fitmax = par1+1.0*par2;
                        sname << "gaussLY - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                        TF1 *gauss = new TF1(sname.str().c_str(),  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitmin,fitmax);
                        gauss->SetParameter(0,par0);
                        gauss->SetParameter(1,par1);
                        gauss->SetParameter(2,par2);
                        spectrumLY->Fit(gauss,"Q","",fitmin,fitmax);
                        //store the mean and sigma in the crystal
                        if(gauss->GetParameter(1) > 0) // otherwise the fit was very wrong..)
                        CurrentCrystal->SetLY(gauss->GetParameter(1),std::abs(gauss->GetParameter(2)));
                        CurrentCrystal->SetLYFit(gauss);


                      }
                      //photons*1.25e6*1.6e-19/156e-15

                      //prepare the photopeak cuts and stuff
                      TCut PhotopeakEnergyCutCorrected = "";
                      TCut PhotopeakEnergyCut  = "";
                      double EnergyCutMin;
                      double EnergyCutMax;
                      int bin3,bin4;
                      double wbin3,wbin4,meanW20;
                      //automatically look for the 511Kev peak to find the photopeak energy cut
                      //find peaks in each crystal spectrum, with TSpectrum
                      if(!backgroundRun)// do it only if this is NOT a background run
                      {
                        TSpectrum *s;
                        s = new TSpectrum(20);
                        // 		Input[i].SumSpectraCanvas->cd(j+1);
                        Int_t CrystalPeaksN = s->Search(spectrumCharge,2,"",0.5);
                        Double_t *CrystalPeaks =  s->GetPositionX();
                        Double_t *CrystalPeaksY = s->GetPositionY();
                        //delete s;
                        float maxPeak = 0.0;
                        int peakID = 0;
                        for (int peakCounter = 0 ; peakCounter < CrystalPeaksN ; peakCounter++ ) //highest peak?
                        {
                          // std::cout << peakCounter << "\t" << CrystalPeaks[peakCounter] << "\t" << CrystalPeaksY[peakCounter] << std::endl;
                          if(CrystalPeaksY[peakCounter] > maxPeak)
                          {
                            maxPeak = CrystalPeaksY[peakCounter];
                            peakID = peakCounter;
                          }
                        }
                        // std::cout << "chosen " << peakID << "\t" << CrystalPeaks[peakID] << "\t" << CrystalPeaksY[peakID] << std::endl;
                        //fit the spectra - TODO use the gaussian plus fermi?
                        if (energyResolution == 0)
                        {
                          if (correctingSaturation)
                            energyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
                          else
                            energyResolution = ENERGY_RESOLUTION;
                        }
                        Float_t par0 = CrystalPeaksY[peakID];
                        Float_t par1 = CrystalPeaks[peakID];
                        Float_t par2 = (CrystalPeaks[peakID]*energyResolution)/2.35;
                        Float_t fitmin = par1-photopeakFitRangeMin*par2;
                        Float_t fitmax = par1+photopeakFitRangeMax*par2;
                        // std::cout << par0   << "\t"
                        //           << par1   << "\t"
                        //           << par2   << "\t"
                        //           << fitmin << "\t"
                        //           << fitmax << std::endl;

                        sname << "gaussCharge - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                        TF1 *gauss = new TF1(sname.str().c_str(),  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitmin,fitmax);
                        gauss->SetParameter(0,par0);
                        gauss->SetParameter(1,par1);
                        gauss->SetParameter(2,par2);

                        spectrumCharge->Fit(gauss,"Q","",fitmin,fitmax);
                        // std::cout << "pars " << par0 << "\t" << par1 << "\t" << par2 << std::endl;

                        if(usingTaggingBench && TagEdgeCalculation)
                        {
                          std::cout << "Crystal "
                                    << CurrentCrystal->GetID()
                                    << " - MPPC "
                                    << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel()
                                    << ", Ratio = "
                                    << fabs(gauss->GetParameter(0)*
                                            gauss->GetParameter(2)*
                                            TMath::Sqrt(2.0*TMath::Pi())) /
                                            tagPeakHgEntries
                                    << std::endl;
                        }
                        //store the mean and sigma in the crystal
                        if(gauss->GetParameter(1) > 0) // otherwise the fit was very wrong..)
                          CurrentCrystal->SetPhotopeak(gauss->GetParameter(1),std::abs(gauss->GetParameter(2)));
                        CurrentCrystal->SetFit(gauss);
                        // 		std::cout << "Photopeak Mean for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakPosition() << std::endl;
                        // 		std::cout << "Photopeak Sigma for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakSigma() << std::endl;
                        // 		std::cout << "Photopeak Energy Resolution FWHM for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakEnergyResolution() << std::endl;
                        //Compute the energy Tcut
                        std::stringstream streamEnergyCut;
                        EnergyCutMin = gauss->GetParameter(1) - photopeakSigmasMin*std::abs(gauss->GetParameter(2));
                        EnergyCutMax = gauss->GetParameter(1) + photopeakSigmasMax*std::abs(gauss->GetParameter(2));
                        streamEnergyCut << SumChannels << " > " << EnergyCutMin << " && " << SumChannels << " < " << EnergyCutMax;
                        PhotopeakEnergyCut  = streamEnergyCut.str().c_str();
                        sname.str("");

                        // then prepare the highlighted spectrum and store it in the crystal
                        sname << "Hg Charge Spectrum - Crystal " << CurrentCrystal->GetID();
                        var << SumChannels << " >> " << sname.str();
                        TH1F* spectrumChargeHighlighted = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                        tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                        spectrumChargeHighlighted->GetXaxis()->SetTitle("ADC Channels");
                        spectrumChargeHighlighted->GetYaxis()->SetTitle("N");
                        CurrentCrystal->SetHighlightedSpectrum(spectrumChargeHighlighted);
                        var.str("");
                        sname.str("");
                      }
                      //-----------------------------------------------------------------------
                      CurrentCrystal->SetSpectrum(spectrumCharge);
                      //

                      // Light sharing plot
                      // draw a NxN canvas with light collected in all channles involved for energy calculation, when trigger channel is in photopeak
                      sname << "Light collected on MPPC "
                            << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel()
                            << " for photoelectric event in crystal "
                            << CurrentCrystal->GetID();
                      var << thisChannel.string.c_str() << " >> " << sname.str();
                      TH1F* LScentralSpectrum = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                      LScentralSpectrum->GetXaxis()->SetTitle("ADC Channels");
                      LScentralSpectrum->GetYaxis()->SetTitle("N");
                      CurrentCrystal->SetLScentralSpectrum(LScentralSpectrum);
                      sname.str("");
                      var.str("");

                      // now run on the neighbour channels that are involved in the energy res calculation, and plot
                      // the charge deposited when central mmpc is in photopeak
                      for(unsigned int iNeig = 0; iNeig < neighbours.size(); iNeig++)
                      {
                        //find the ID of this neighbour in the detector index
                        int iNeighDet;
                        for(unsigned int iDet = 0; iDet < detector.size(); iDet++)
                        {
                          if(detector[iDet].digitizerChannel == neighbours[iNeig])
                          {
                            iNeighDet = iDet;
                          }
                        }
                        int  neighID = 0;
                        for(unsigned int iNeigh = 0; iNeigh < neighbourChannels.size(); iNeigh++)
                        {
                          if( detector[iNeighDet].digitizerChannel == neighbourChannels[iNeigh].detectorIndex)
                          {
                            neighID = iNeigh;
                          }
                        }

                        // Light sharing plot
                        sname << "Light collected on MPPC "
                              << detector[iNeighDet].label
                              << " for photoelectric event in crystal "
                              << CurrentCrystal->GetID();
                        var << neighbourChannels[neighID].string.c_str() << " >> " << sname.str();
                        TH1F* LSSpectrum = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                        tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                        LSSpectrum->GetXaxis()->SetTitle("ADC Channels");
                        LSSpectrum->GetYaxis()->SetTitle("N");
                        sname.str("");
                        var.str("");

                        int plotPos = -1;
                        for(int iTimeMppc = 0 ; iTimeMppc < nmppcx ; iTimeMppc++)
                        {
                          for(int jTimeMppc = 0 ; jTimeMppc < nmppcy ; jTimeMppc++)
                          {
                            if(mppc[(iModule*nmppcx)+iTimeMppc][(jModule*nmppcy)+jTimeMppc]->GetDigitizerChannel() == neighbours[iNeig])
                            {
                              plotPos = mppc[(iModule*nmppcx)+iTimeMppc][(jModule*nmppcy)+jTimeMppc]->GetCanvasPosition();
                            }
                          }
                        }
                        if(plotPos != -1)
                        {
                          CurrentCrystal->AddLSSpectrum(LSSpectrum,plotPos);
                        }
                      }




                      //w histogram with cut on crystal events, xyz and trigger channel and cut on photopeak
                      //if it's a background run, it won't cut on photopeak (since there will be none and the string will be empty)
                      sname << "W histogram - Crystal " << CurrentCrystal->GetID();
                      //var << "(ch" << channel << "/(" << SumChannels << ")) >> " << sname.str();
                      var  << FloodZ.str() << " >> " << sname.str();
                      TH1F* spectrumHistoW = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                      spectrumHistoW->GetXaxis()->SetTitle("W");
                      spectrumHistoW->GetYaxis()->SetTitle("N");
                      // 		int bin1 = spectrumHistoW->FindFirstBinAbove(spectrumHistoW->GetMaximum()/2.0);
                      // 		int bin2 = spectrumHistoW->FindLastBinAbove(spectrumHistoW->GetMaximum()/2.0);
                      bin3 = spectrumHistoW->FindFirstBinAbove(wThreshold*spectrumHistoW->GetMaximum());
                      bin4 = spectrumHistoW->FindLastBinAbove(wThreshold*spectrumHistoW->GetMaximum());
                      // for the linear fit, let's try to stay in the really linear part of the crystal (in terms of ADCch vs w)
                      // so the fit rang will not be between these two points but closer to the center,
                      // for the moment let's say the central 50% of the w width
                      double wmin = spectrumHistoW->GetBinCenter(bin3);
                      double wmax = spectrumHistoW->GetBinCenter(bin4);
                      //DEBUG
                      meanW20 = (wmax + wmin) / 2.0;
                      wbin3 = wmin + (wmax-wmin)*energyCorrectionMin;
                      wbin4 = wmin + (wmax-wmin)*energyCorrectionMax;
                      std::stringstream ssCut20w;
                      ssCut20w << "(" << FloodZ.str() << ") > " << spectrumHistoW->GetBinCenter(bin3) << " && " << "(" << FloodZ.str() << ") < "<<  spectrumHistoW->GetBinCenter(bin4);
                      TCut w20percCut = ssCut20w.str().c_str();  //cut for w to get only the "relevant" part - TODO find a reasonable way to define this
                      CurrentCrystal->SetW20percCut(w20percCut);
                      CurrentCrystal->SetHistoW(spectrumHistoW);
                      var.str("");
                      sname.str("");
                      //
                      if(usingTaggingBench || calcDoiResWithCalibration) //done on the standard W because there's no ADCvsW correction in tagging bench
                      {
                        sname << "gaussW histogram - Crystal " << CurrentCrystal->GetID();
                        TF1 *gaussW = new TF1(sname.str().c_str(),  "gaus",histo3Dmin,histo3Dmax);
                        int binmax = spectrumHistoW->GetMaximumBin();
                        double maximum = spectrumHistoW->GetXaxis()->GetBinCenter(binmax);
                        int nentries= spectrumHistoW->GetEntries();
                        gaussW->SetParameter(0,maximum);
                        gaussW->SetParameter(1,spectrumHistoW->GetMean());
                        gaussW->SetParameter(2,spectrumHistoW->GetRMS());
                        Int_t fitStatus = spectrumHistoW->Fit(sname.str().c_str(),"QR");
                        if(!fitStatus)
                        doiFile << CurrentCrystal->GetI() + doiColumnOffset << "\t" << CurrentCrystal->GetJ() << "\t" << gaussW->GetParameter(1) << "\t" << taggingPosition <<"\t" << gaussW->GetParameter(2)/TMath::Sqrt(nentries) <<"\t"<<TMath::Sqrt(nentries)<< std::endl;
                        CurrentCrystal->SetHistoWfit(gaussW);
                        sname.str("");

                        // if(calcDoiResWithDelta) // FIXME for now with this approach it is not possible to have different ZPositions in the same dataset, so this calcDoiResWithDelta should not be used!
                        // {
                        //   for(unsigned int kAltDoi = 0; kAltDoi < inputDoi.size(); kAltDoi++)
                        //   {
                        //     int ik = inputDoi[kAltDoi].i;
                        //     int jk = inputDoi[kAltDoi].j;
                        //     if( ((CurrentCrystal->GetI() + doiColumnOffset ) == ik) && CurrentCrystal->GetJ() == jk)
                        //     {
                        //       sname << "Alternate DOI res - Crystal " << CurrentCrystal->GetID();
                        //       var << "(FloodZ * " << inputDoi[kAltDoi].m << "+" << inputDoi[kAltDoi].q << ") - (ZPosition  ) >> " << sname.str();
                        //       TH1F* spectrumAltDoiRes = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,-10,10); //FIXME boundaries hardcoded
                        //       tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                        //       spectrumAltDoiRes->GetXaxis()->SetTitle("Doi_calculated - DoiTagPosition [mm]");
                        //       spectrumAltDoiRes->GetYaxis()->SetTitle("N");
                        //       var.str("");
                        //       sname.str("");
                        //
                        //       sname << "gaussFit Alternate - Crystal " << CurrentCrystal->GetID();
                        //       TF1 *gaussFitAlt = new TF1(sname.str().c_str(),  "gaus",-10,10); //FIXME boundaries hardcoded
                        //       Int_t fitStatusAlt = spectrumAltDoiRes->Fit(sname.str().c_str(),"QR");
                        //       if(!fitStatusAlt)
                        //       {
                        //         AltDoiFile << ik << " " << jk << " " << taggingPosition << " " << gaussFitAlt->GetParameter(2)*2.355 << std::endl;
                        //       }
                        //
                        //
                        //       CurrentCrystal->SetHistoAltDoiRes(spectrumAltDoiRes);
                        //       CurrentCrystal->SetHistoAltDoiFit(gaussFitAlt);
                        //       sname.str("");
                        //       // for(int iik =0 ; iik < pointsFromDoi ; iik++)
                        //       // {
                        //       //   sigmaWdoiCentral->Fill(inputDoi[kAltDoi].sw[iik] * inputDoi[kAltDoi].sqrt_nentries[iik]);
                        //       // }
                        //     }
                        //   }
                        // }
                      }



                      if(!backgroundRun)
                      {
                        //histogram of w versus adc channels
                        //it will be useful for doi correction
                        //long long int nPoints;
                        sname << "ADC channels vs. W - Crystal " << CurrentCrystal->GetID();
                        var << SumChannels << " : " <<FloodZ.str() << " >> " << sname.str();
                        TH2F* spectrum2dADCversusW = new TH2F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax,histo1Dbins,0,histo1Dmax);
                        tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+w20percCut,"COLZ");
                        spectrum2dADCversusW->GetXaxis()->SetTitle("W");
                        spectrum2dADCversusW->GetYaxis()->SetTitle("ADC channels");
                        double parM;
                        double parQ;
                        var.str("");
                        sname.str("");


                        if(correctingForDOI)
                        {
                          //mayhem
                          //i.e. use FitSlicesY to get the gaussian fit of each slices of the TH2F. Sliced in x (bin by bin)
                          // first define the gaussian function
                          //TF1 *gaussFitSlice = new TF1("gaussFitSlice","[0]*exp(-0.5*((x-[1])/[2])**2)",EnergyCutMin,EnergyCutMax);
                          //gaussFitSlice->SetParameter(1,(EnergyCutMin+EnergyCutMax)/2.0);
                          //gaussFitSlice->SetParameter(2,0.15*(EnergyCutMin+EnergyCutMax)/2.0);
                          // 		    gaussFitSlice->SetRange();
                          spectrum2dADCversusW->FitSlicesY(0, bin3, bin4, 0, "QNRG5S");
                          sname << spectrum2dADCversusW->GetName() << "_1";
                          TH1D *spectrum2d_1 = (TH1D*)gDirectory->Get(sname.str().c_str()); // _1 is the TH1D automatically created by ROOT when FitSlicesX is called, holding the TH1F of the mean values
                          sname.str("");
                          sname << "linearCrystal - Crystal " << CurrentCrystal->GetID();
                          TF1 *linearCrystal = new TF1(sname.str().c_str(),  "[0]*x + [1]",wbin3,wbin4);
                          spectrum2d_1->Fit(sname.str().c_str(),"QR");

                          parM = linearCrystal->GetParameter(0); // m parameter for the linear fit to correct energy res for DOI
                          parQ = linearCrystal->GetParameter(1); // q parameter for the linear fit to correct energy res for DOI
                          CurrentCrystal->SetSlicesMean(spectrum2d_1);
                          CurrentCrystal->SetSlicesMeanFit(linearCrystal);
                          sname.str("");


                        }

                        // 		spectrum2dADCversusW->SetName(sname.str().c_str());
                        CurrentCrystal->SetADCversusW(spectrum2dADCversusW);

                        if(correctingForDOI)
                        {
                          //spectrum corrected for DOI
                          sname << "Charge Spectrum Corrected - Crystal " << CurrentCrystal->GetID();
                          std::stringstream baseVar;

                          // 		    baseVar << "(("  <<  SumChannels<< " ) - ( ( FloodZ - " <<  meanW20 << " ) * ( " << parM << ") ))";
                          //calculate ADC(w0), it's always the same. w0=meanw20
                          // 		    double adcw0 = parM * meanW20 + parQ;

                          baseVar << "(("  <<  SumChannels<< " ) * ( (" << parM * meanW20 + parQ << ") / (" << parM << "* " << FloodZ.str() << " + "<< parQ << " )))";
                          // 		    std::cout << std::endl;
                          // 		    std::cout << bin3 << " " << bin4 << " "  << baseVar.str() << std::endl;
                          var << baseVar.str() << " >> " << sname.str();
                          TH1F* spectrumChargeCorrected = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                          tree->Draw(var.str().c_str(),CrystalCut,"COLZ");
                          spectrumChargeCorrected->GetXaxis()->SetTitle("ADC channels");
                          spectrumChargeCorrected->GetYaxis()->SetTitle("Counts");
                          CurrentCrystal->SetCorrectedSpectrum(spectrumChargeCorrected);
                          var.str("");
                          sname.str("");


                          spectrumChargeCorrected->GetXaxis()->SetRangeUser(peakSearchRangeMin,peakSearchRangeMax);
                          TSpectrum *s_corr;
                          s_corr = new TSpectrum(5);
                          // 		Input[i].SumSpectraCanvas->cd(j+1);


                          Int_t CrystalPeaksN_corr = s_corr->Search(spectrumChargeCorrected,1,"",0.5);
                          Double_t *CrystalPeaks_corr =  s_corr->GetPositionX();
                          Double_t *CrystalPeaksY_corr = s_corr->GetPositionY();
                          // 		  delete s_corr;
                          float maxPeak_corr = 0.0;
                          int peakID_corr = 0;

                          // std::cout << "-----------------------------------------------" << std::endl;
                          for (int peakCounter = 0 ; peakCounter < CrystalPeaksN_corr ; peakCounter++ )
                          {
                            // //DEBUG
                            // std::cout << peakCounter << " " << CrystalPeaks_corr[peakID_corr] << " " << CrystalPeaksY_corr[peakID_corr] << std::endl;
                            // if( (CrystalPeaks_corr[peakCounter] >= peakSearchRangeMin) && (CrystalPeaks_corr[peakCounter] <= peakSearchRangeMax) ) //look for the 511 peak only in the selected range (which could be all histogram if nothing is set)
                            // {
                            if(CrystalPeaksY_corr[peakCounter] > maxPeak_corr) // then look for tallest peak in the range
                            {
                              maxPeak_corr = CrystalPeaksY_corr[peakCounter];
                              peakID_corr = peakCounter;
                            }
                            // }
                          }
                          // //DEBUG
                          // std::cout << "chosen - " << peakID_corr << " " << CrystalPeaks_corr[peakID_corr] << " " << CrystalPeaksY_corr[peakID_corr] << std::endl;

                          //std::cout << CrystalPeaks[0] << std::endl;
                          //std::cout << CrystalPeaksY[0] << std::endl;
                          //fit the spectra - TODO use the gaussian plus fermi?
                          //float energyResolution;
                          if (energyResolution == 0)
                          {
                            if (correctingSaturation)
                              energyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
                            else
                              energyResolution = ENERGY_RESOLUTION;
                          }
                          float par0_corr = CrystalPeaksY_corr[peakID_corr];
                          float par1_corr = CrystalPeaks_corr[peakID_corr];
                          float par2_corr = (CrystalPeaks_corr[peakID_corr]*energyResolution)/2.35;
                          float fitmin_corr = par1_corr-photopeakFitRangeMin*par2_corr;
                          float fitmax_corr = par1_corr+photopeakFitRangeMax*par2_corr;

                          //DEBUG
                          // std::cout << fitmin_corr << " " << fitmax_corr << std::endl;

                          sname << "gauss_corr - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                          TF1 *gauss_corr = new TF1(sname.str().c_str(),  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitmin_corr,fitmax_corr);
                          gauss_corr->SetParameter(0,par0_corr);
                          gauss_corr->SetParameter(1,par1_corr);
                          gauss_corr->SetParameter(2,par2_corr);
                          spectrumChargeCorrected->Fit(gauss_corr,"Q","",fitmin_corr,fitmax_corr);

                          spectrumChargeCorrected->GetXaxis()->SetRangeUser(1,histo1Dmax);

                          //store the mean and sigma in the crystal
                          if(gauss_corr->GetParameter(1) > 0) // otherwise the fit was very wrong..)
                            CurrentCrystal->SetPhotopeakCorrected(gauss_corr->GetParameter(1),std::abs(gauss_corr->GetParameter(2)));
                          // 		  else
                          // 		    CurrentCrystal->SetPhotopeakCorrected(1,1);
                          CurrentCrystal->SetFitCorrected(gauss_corr);
                          // 		std::cout << "Photopeak Mean for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakPosition() << std::endl;
                          // 		std::cout << "Photopeak Sigma for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakSigma() << std::endl;
                          // 		std::cout << "Photopeak Energy Resolution FWHM for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakEnergyResolution() << std::endl;
                          //Compute the energy Tcut
                          std::stringstream streamEnergyCutCorrected;
                          streamEnergyCutCorrected << baseVar.str() << " > " << gauss_corr->GetParameter(1) - photopeakSigmasMin*std::abs(gauss_corr->GetParameter(2)) << " && " << baseVar.str() << " < " << gauss_corr->GetParameter(1) + photopeakSigmasMax*std::abs(gauss_corr->GetParameter(2));
                          PhotopeakEnergyCutCorrected = streamEnergyCutCorrected.str().c_str();

                          // then PhotopeakEnergyCut becomes PhotopeakEnergyCutCorrected

                          PhotopeakEnergyCut = PhotopeakEnergyCutCorrected;
                          // 		CurrentCrystal->SetSpectrum(*spectrum);
                          sname.str("");

                          //then prepare the highlighted spectrum and store it in the crystal
                          sname << "Hg Charge Spectrum Correctd - Crystal " << CurrentCrystal->GetID();
                          // 			var << "("  <<  SumChannels<< " ) - ( ( FloodZ - " <<  meanW20 << " ) * ( " << parM << ") ) >> " << sname.str();
                          var << baseVar.str() << " >> " << sname.str();
                          TH1F* spectrumChargeCorrectedHighlighted = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                          tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                          spectrumChargeCorrectedHighlighted->GetXaxis()->SetTitle("ADC Channels");
                          spectrumChargeCorrectedHighlighted->GetYaxis()->SetTitle("N");
                          CurrentCrystal->SetHighlightedSpectrumCorrected(spectrumChargeCorrectedHighlighted);
                          var.str("");
                          sname.str("");

                        }
                      }

                      // save the TCuts
                      TCut CrystalCutWithoutCutG = BasicCut + CutTrigger;               // useful for the calcDoiResWithCalibration external part
                      CrystalCut.SetName("CrystalCut");
                      PhotopeakEnergyCut.SetName("PhotopeakEnergyCut");
                      CrystalCutWithoutCutG.SetName("CrystalCutWithoutCutG");
                      CurrentCrystal->SetCrystalCut(CrystalCut);                        // this is BasicCut + CutTrigger + cutg[0] + cutg[1] - so "geometrical" without energy constrains
                      CurrentCrystal->SetPhotopeakEnergyCut(PhotopeakEnergyCut);        // energy cut, events in the photopeak
                      CurrentCrystal->SetCrystalCutWithoutCutG(CrystalCutWithoutCutG);  // this is BasicCut + CutTrigger


                      // Histogram 2d of time evolution
                      sname << "ADC channels vs. Time - Crystal " << CurrentCrystal->GetID();
                      var << SumChannels << ":ExtendedTimeTag >> " << sname.str();
                      TH2F* spectrum2dVersusTime = new TH2F(sname.str().c_str(),sname.str().c_str(),250,0,tree->GetMaximum("ExtendedTimeTag"),histo1Dbins,0,histo1Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut,"COLZ");
                      spectrum2dVersusTime->GetXaxis()->SetTitle("ExtendedTimeTag");
                      spectrum2dVersusTime->GetYaxis()->SetTitle("ADC channels");
                      CurrentCrystal->SetVersusTime(spectrum2dVersusTime);
                      var.str("");
                      sname.str("");

                      // Time evolution for W
                      sname << "W vs. Time - Crystal " << CurrentCrystal->GetID();
                      var << FloodZ.str() << ":ExtendedTimeTag >> " << sname.str();
                      TH2F* spectrum2dWversusTime = new TH2F(sname.str().c_str(),sname.str().c_str(),250,0,tree->GetMaximum("ExtendedTimeTag"),wHistogramsBins,histo3Dmin,histo3Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut,"COLZ");
                      spectrum2dWversusTime->GetXaxis()->SetTitle("ExtendedTimeTag");
                      spectrum2dWversusTime->GetYaxis()->SetTitle("W");
                      CurrentCrystal->SetWversusTime(spectrum2dWversusTime);
                      var.str("");
                      sname.str("");

                      //histogram of w versus adc channels - this time without the photopeak cut (so it looks nicer in the paper...)
                      sname << "Complete ADC channels vs. W - Crystal " << CurrentCrystal->GetID();
                      var << SumChannels << ":" << FloodZ.str() << " >> " << sname.str() ;
                      TH2F* spectrum2dADCversusWComplete = new TH2F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax,histo1Dbins,0,histo1Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut,"COLZ");
                      spectrum2dADCversusWComplete->GetXaxis()->SetTitle("W");
                      spectrum2dADCversusWComplete->GetYaxis()->SetTitle("ADC channels");
                      spectrum2dADCversusWComplete->SetName(sname.str().c_str());
                      CurrentCrystal->SetADCversusWComplete(spectrum2dADCversusWComplete);
                      var.str("");
                      sname.str("");


                      // if(correctingForDOI)
                      // {
                      //   // alternative method, even more mayhem
                      //   // manually split the 2d wistogram in slices and fit the individual slices, finding the
                      //   // photopeak. this does not require the w20percCut
                      //
                      //
                      //
                      //
                      //
                      // }

                      // w histogram with energy cut on the corrected spectrum and
                      // fit the w histogram with the thetaFunction -> m and q of correlation line
                      // fit on the rise of w function with gaussian -> sigma w
                      // no cut on the photopeak will be applied if it's a background run (since there will be no photopeak)
                      sname << "W histogram Corrected - Crystal " << CurrentCrystal->GetID();
                      var << "(" << FloodZ.str() << ") >> " << sname.str();
                      TH1F* spectrumHistoWCorrected = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                      spectrumHistoWCorrected->GetXaxis()->SetTitle("W");
                      spectrumHistoWCorrected->GetYaxis()->SetTitle("N");
                      var.str("");
                      sname.str("");

                      CurrentCrystal->SetHistoWCorrected(spectrumHistoWCorrected);

                      TGraph *calibrationGraphChosen = NULL;
                      float centralW = -1;
                      float beginW = -1;
                      float endW = -1;

                      if(!calcDoiResWithCalibration) // only if it's not calcDoiResWithCalibration, because in that case we load the calibration TGraph from another file
                      {
                        var.str("");
                        sname.str("");

                        //cumulative function
                        //first, take the w corrected and normalize to integral -> get the PDF
                        TH1F *pdfW = (TH1F*) spectrumHistoWCorrected->Clone();
                        sname << "PDF W histogram Corrected - Crystal " << CurrentCrystal->GetID();
                        pdfW->SetName(sname.str().c_str());
                        pdfW->SetTitle(sname.str().c_str());
                        double integral = pdfW->Integral();
                        pdfW->Scale(wHistogramsBins*(1.0/integral));
                        CurrentCrystal->SetPdfW(pdfW);
                        sname.str("");
                        // then do the comulative histo
                        //and the calibration Graph

                        sname << "Cumulative W histogram Corrected - Crystal " << CurrentCrystal->GetID();
                        // 		    var << "(ch" << channel << "/(" << SumChannels << ")) >> " << sname.str();
                        TH1F* cumulativeW = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax);
                        cumulativeW->GetXaxis()->SetTitle("W");
                        double sumPdf = 0;
                        // 		    std::ofstream wfile;
                        // 		    wfile.open("cumulative2_2.dat", std::ofstream::out);
                        std::vector<double> calibrationW;
                        std::vector<double> calibrationZ;
                        // 		    double lambda511 = 12.195; //everything in mm

                        for(int iPdfHisto = 0 ; iPdfHisto < wHistogramsBins; iPdfHisto++)
                        {
                          sumPdf += (pdfW->GetBinContent(iPdfHisto+1))/wHistogramsBins;
                          calibrationW.push_back(cumulativeW->GetBinCenter(iPdfHisto+1));
                          if(backgroundRun | lateralRun)
                          {
                            calibrationZ.push_back(-(sumPdf*crystalz) + crystalz); // if it's a background run, interaction probability is constant everywhere, so z is just a rescale of the cumulative
                          }
                          else
                          {
                            // if it's a far source run, interaction probability is exponential, and relation between z and w is give by
                            // z = L + l * ln( 1 - (1 - exp(-L/l))*integral_0^w(PDF(w)dw) )
                            // where
                            // L = crystal length
                            // l = interaction length of 511 gammas in lyso (from literature, 12.195 mm)
                            // integral_0^w(PDF(w)dw) = the cumulative of PDF(w) from 0 to w
                            calibrationZ.push_back( crystalz +  lambda511 * TMath::Log( 1.0 - (1.0 - TMath::Exp(-(crystalz/lambda511)) )* sumPdf )  );
                          }
                          // 		      wfile << cumulativeW->GetBinCenter(iPdfHisto+1) << " " << -(sumPdf*15) + 15.0 << std::endl;
                          cumulativeW->Fill(cumulativeW->GetBinCenter(iPdfHisto+1),sumPdf);
                        }

                        TGraph *calibrationGraph = new TGraph(calibrationW.size(),&calibrationW[0],&calibrationZ[0]);
                        sname.str("");
                        sname << "Calibration Plot - Crystal " << CurrentCrystal->GetID();
                        calibrationGraph->SetTitle(sname.str().c_str());
                        calibrationGraph->SetName(sname.str().c_str());
                        calibrationGraph->GetXaxis()->SetTitle("W");
                        calibrationGraph->GetYaxis()->SetTitle("Z [mm]");
                        sname.str("");
                        CurrentCrystal->SetCalibrationGraph(calibrationGraph);
                        calibrationGraphChosen = calibrationGraph;

                        //w(z) for compton analysis
                        TGraph *wzgraph = new TGraph(calibrationW.size(),&calibrationZ[0],&calibrationW[0]);
                        sname.str("");
                        sname << "w(z) Plot - Crystal " << CurrentCrystal->GetID();
                        wzgraph->SetTitle(sname.str().c_str());
                        wzgraph->SetName(sname.str().c_str());
                        wzgraph->GetXaxis()->SetTitle("Z [mm]");
                        wzgraph->GetYaxis()->SetTitle("W");
                        centralW = wzgraph->Eval(crystalz/2.0);
                        beginW = wzgraph->Eval(crystalz - marginWZgraph);
                        endW = wzgraph->Eval(marginWZgraph);
                        // std::cout << beginW << "\t" << endW << std::endl;

                        sname.str("");

                        CurrentCrystal->SetWZgraph(wzgraph);
                        CurrentCrystal->SetCumulativeW(cumulativeW);
                      }
                      else // look for the calibGraph in an external calibration file
                      {
                        //save the crystal identification
                        std::stringstream MppcDirStream;
                        MppcDirStream << "MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " - 0.0-" << iMppc << "." << jMppc;
                        std::stringstream CrystalDirStream;
                        CrystalDirStream << "Crystal " <<  CurrentCrystal->GetID();
                        std::stringstream directory;
                        directory << "Module 0.0/" << MppcDirStream.str() << "/" << CrystalDirStream.str();
                        bool dir_exists = calibrationFile->cd(directory.str().c_str());
                        TCanvas* C_graph = NULL;
                        TGraph *calibGraph = NULL;
                        if(dir_exists)
                        {
                          //take the graph
                          std::stringstream stream;
                          stream << "Calibration Plot - Crystal " << CurrentCrystal->GetID();
                          C_graph = (TCanvas*) gDirectory->Get(stream.str().c_str());
                          if(C_graph)
                            calibGraph = (TGraph*) C_graph->GetPrimitive(stream.str().c_str());

                          var.str("");
                          sname.str("");
                        }
                        CurrentCrystal->SetCalibrationGraph(calibGraph);
                        calibrationGraphChosen = calibGraph;
                      }


                      if(usingTaggingBench || taggingForTiming)
                      {
                        if( (digitizerType == 1) || (digitizerType == 2)) //only for digitizers with timing capabilities
                        {
                          //prepare three std::vector, timingChannel, meanDelta and sigmaTiming
                          // they will hold the t channel number, mean difference in t and sigma from crystalball fits

                          std::vector<int> tChannelsForPolishedCorrection;
                          std::vector<double> meanForPolishedCorrection;
                          std::vector<double> fwhmForPolishedCorrection;


                          // basic CTR with respect to an external tagging crystal, cutting on photopeak of this crystal
                          // var.str("");
                          // sname.str("");
                          CurrentCrystal->SetTimingChannel(detector[thisChannelID].timingChannel);

                          sname << "Basic CTR - Crystal " << CurrentCrystal->GetID();
                          var << "t" << detector[thisChannelID].timingChannel
                              << "- t" << taggingCrystalTimingChannel
                              << " >> " << sname.str();
                          // std::cout << sname.str() << std::endl
                          // << var.str()   << std::endl
                          // << CrystalCut+PhotopeakEnergyCut << std::endl;

                          // Always esclude from production of calibration histograms the events where the timing channels involved in the plot are both = 0
                          // these cuts won't need to be exported to the calibration file explicitly because their absence is reflected in the resulting histograms
                          // when Analyzing a dataset based on a calibration data, of course this dataset will have to be filtered to exclude explicitly these events
                          // although they probably self exclude because of the other cuts

                          TCut noZerosCut;
                          std::stringstream sNoZerosCut;
                          // noZerosCut for Basic CTR
                          sNoZerosCut <<  "(t" << detector[thisChannelID].timingChannel << "!= 0) && ("
                                      << "t" << taggingCrystalTimingChannel << "!= 0)";
                          noZerosCut = sNoZerosCut.str().c_str();
                          sNoZerosCut.str("");

                          TH1F* aSpectrum = new TH1F(sname.str().c_str(),sname.str().c_str(),CTRbins,CTRmin,CTRmax);
                          tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+noZerosCut);
                          aSpectrum->GetXaxis()->SetTitle("Time [S]");
                          aSpectrum->GetYaxis()->SetTitle("N");

                          var.str("");
                          sname.str("");

                          //do a preliminary fit with gauss
                          // TCanvas *cTemp  = new TCanvas("temp","temp");
                          // TF1 *gaussDummyBasic = new TF1("gaussDummy","gaus");
                          // aSpectrum->Fit(gaussDummyBasic,"QN");
                          //

                          //
                          // double f1min = aSpectrum->GetXaxis()->GetXmin();
                          // double f1max = aSpectrum->GetXaxis()->GetXmax();
                          // // std::cout << f1min << " " << f1max << std::endl;
                          // TF1* f1  = new TF1("f1","crystalball");
                          // f1->SetLineColor(kRed);
                          // f1->SetParameters(gaussDummyBasic->GetParameter(0),gaussDummyBasic->GetParameter(1),gaussDummyBasic->GetParameter(2),1,3);
                          // double fitMin = gaussDummyBasic->GetParameter(1)
                          //                 - fitPercMin*(gaussDummyBasic->GetParameter(2));
                          // double fitMax = gaussDummyBasic->GetParameter(1)
                          //                 + fitPercMax*(gaussDummyBasic->GetParameter(2));
                          // if(fitMin < f1min)
                          // {
                          //   fitMin = f1min;
                          // }
                          // if(fitMax > f1max)
                          // {
                          //   fitMax = f1max;
                          // }
                          // aSpectrum->Fit(f1,"Q","",fitMin,fitMax);
                          double fitPercMin = 5.0;
                          double fitPercMax = 6.0;
                          int divisions = 10000;
                          double res[2];
                          extractFromCrystalBall(aSpectrum,fitPercMin,fitPercMax,divisions,res);
                          tChannelsForPolishedCorrection.push_back(detector[thisChannelID].timingChannel);
                          meanForPolishedCorrection.push_back(res[0]);
                          fwhmForPolishedCorrection.push_back(res[1]);

                          CurrentCrystal->SetDeltaTimeWRTTagging(aSpectrum);

                          // scatter plot of the above delta and w (which correlates to doi in depolished matrix)
                          sname << "Delta T vs. W - Crystal " << CurrentCrystal->GetID();
                          var << "t" << detector[thisChannelID].timingChannel
                              << "- t" << taggingCrystalTimingChannel
                              << " : "
                              << FloodZ.str()
                              << " >> " << sname.str() ;

                          //noZerosCut
                          sNoZerosCut << "(t" << detector[thisChannelID].timingChannel << "!= 0) && ("
                                      << "t" << taggingCrystalTimingChannel << "!= 0)";
                          noZerosCut = sNoZerosCut.str().c_str();
                          sNoZerosCut.str("");

                          // limits of the 2d histo are derived from W histo and ctr histo
                          TH2F* spectrumCrystalDeltaTvsW = new TH2F(sname.str().c_str(),sname.str().c_str(),
                                                                    wHistogramsBins,
                                                                    spectrumHistoW->GetMean() - 3.0*spectrumHistoW->GetRMS(),
                                                                    spectrumHistoW->GetMean() + 3.0*spectrumHistoW->GetRMS(),
                                                                    CTRbins,
                                                                    aSpectrum->GetMean() - 3.0*aSpectrum->GetRMS(),
                                                                    aSpectrum->GetMean() + 3.0*aSpectrum->GetRMS());
                          tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+noZerosCut,"COLZ");
                          spectrumCrystalDeltaTvsW->GetXaxis()->SetTitle("W");
                          spectrumCrystalDeltaTvsW->GetYaxis()->SetTitle("T crystal - T tagging [S]");
                          CurrentCrystal->SetDeltaTvsW(spectrumCrystalDeltaTvsW);
                          var.str("");
                          sname.str("");

                          //delta t vs "amplitude" (i.e. total charge deposited) - so far it makes sense only for polished crystals
                          sname << "Delta T vs. ch" << detector[thisChannelID].digitizerChannel << " - Crystal " << CurrentCrystal->GetID();
                          var << "t" << detector[thisChannelID].timingChannel
                              << "- t" << taggingCrystalTimingChannel
                              << " : "
                              << thisChannel.string
                              << " >> " << sname.str() ;
                          //noZerosCut
                          sNoZerosCut << "(t" << detector[thisChannelID].timingChannel << "!= 0) && ("
                                      << "t" << taggingCrystalTimingChannel << "!= 0)";
                          noZerosCut = sNoZerosCut.str().c_str();
                          sNoZerosCut.str("");

                          TH2F* spectrumCrystalDeltaTvsCH = new TH2F(sname.str().c_str(),sname.str().c_str(),
                                                                     histo1Dbins,
                                                                     LScentralSpectrum->GetMean() - 3.0*LScentralSpectrum->GetRMS(),
                                                                     LScentralSpectrum->GetMean() + 3.0*LScentralSpectrum->GetRMS(),
                                                                     CTRbins,
                                                                     aSpectrum->GetMean() - 3.0*aSpectrum->GetRMS(),
                                                                     aSpectrum->GetMean() + 3.0*aSpectrum->GetRMS());
                          tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+noZerosCut,"COLZ");
                          sname.str("");
                          sname << "ch" << detector[thisChannelID].digitizerChannel << " [ADC channels]";
                          spectrumCrystalDeltaTvsCH->GetXaxis()->SetTitle(sname.str().c_str());
                          spectrumCrystalDeltaTvsCH->GetYaxis()->SetTitle("T crystal - T tagging [S]");
                          CurrentCrystal->SetDeltaTvsCH(spectrumCrystalDeltaTvsCH);
                          var.str("");
                          sname.str("");

                          // spectrumCrystalDeltaTvsW->FitSlicesY(0, 0, -1, 0, "QNRG5S");
                          // sname << spectrumCrystalDeltaTvsW->GetName() << "_1";
                          // TH1D *spectrumCrystalDeltaTvsW_1 = (TH1D*)gDirectory->Get(sname.str().c_str()); // _1 is the TH1D automatically created by ROOT when FitSlicesY is called, holding the TH1F of the mean values
                          // sname.str("");
                          // sname << "linearFit CTR - Crystal " << CurrentCrystal->GetID();
                          // TF1 *linearFitCTR = new TF1(sname.str().c_str(),  "[0]*x + [1]",histo3Dmin,histo3Dmax);
                          // spectrumCrystalDeltaTvsW_1->Fit(sname.str().c_str(),"QR");
                          // CurrentCrystal->SetDeltaTvsWSlice(spectrumCrystalDeltaTvsW_1);
                          // CurrentCrystal->SetDeltaTvsWFit(linearFitCTR);
                          // var.str("");
                          // sname.str("");



                          if(timingCorrectionForPolished)
                          {



                            // std::vector<int> DelayTimingChannelsNum;

                            //get a TGraph
                            std::vector<float> delta_X,delta_X_core;
                            std::vector<float> W_Y,W_Y_core;
                            std::vector<float> rmsBasic_Y,rmsBasic_Y_core;

                            std::cout << "Performing timing correction analysis"  << std::endl;
                            std::cout << "Analyzing interaction crystal..."  << std::endl;

                            if(timingCorrection)
                            {
                              //---------------------------------------------//
                              // Interaction crystal w slicing               //
                              //---------------------------------------------//

                              //Tgraph from mean and rms of "slices"
                              // std::cout << beginW << "\t" << endW << std::endl;
                              for(int iBin = 0; iBin < WrangeBinsForTiming; iBin++) // -1 and +1 are put to include the w limits
                              {
                                //slicing of central crystal ctr vs w

                                // define the Limits and mid point for this w slice
                                Float_t wLowerLimit = beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming);
                                Float_t wUpperLimit = beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming);
                                Float_t midW        = beginW + (((iBin+0.5)*(endW - beginW))/WrangeBinsForTiming);

                                std::stringstream sCut;
                                sCut << FloodZ.str() << " > "
                                << wLowerLimit
                                << " && " << FloodZ.str() << " < "
                                << wUpperLimit
                                << " && ("
                                << "t" << detector[thisChannelID].timingChannel
                                << "- t" << taggingCrystalTimingChannel << ") >"
                                << CTRmin
                                << " && ("
                                << "t" << detector[thisChannelID].timingChannel
                                << "- t" << taggingCrystalTimingChannel << ") < "
                                << CTRmax;
                                TCut wCut = sCut.str().c_str();

                                //noZerosCut
                                sNoZerosCut << "(t" << detector[thisChannelID].timingChannel << "!= 0) && ("
                                            << "t" << taggingCrystalTimingChannel << "!= 0)";
                                noZerosCut = sNoZerosCut.str().c_str();
                                sNoZerosCut.str("");

                                sname.str("");

                                var.str("");
                                sname.str("");
                                sname << "DeltaT_" << iBin << " - Crystal " << CurrentCrystal->GetID();
                                TH1F *tempHisto = new TH1F(sname.str().c_str(),sname.str().c_str(),CTRbins,CTRmin,CTRmax);
                                var << "t" << detector[thisChannelID].timingChannel
                                << "- t" << taggingCrystalTimingChannel
                                << "  >> "
                                << sname.str() ;

                                tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+wCut+noZerosCut);

                                double fitPercMin = 5.0;
                                double fitPercMax = 6.0;
                                int divisions = 10000;
                                double res[2];
                                extractFromCrystalBall(tempHisto,fitPercMin,fitPercMax,divisions,res);
                                // tChannelsForPolishedCorrection.push_back(detector[thisChannelID].timingChannel);
                                // meanForPolishedCorrection.push_back(res[0]);
                                // fwhmForPolishedCorrection.push_back(res[1]);

                                // TF1* crystalball  = new TF1("crystalball","crystalball");
                                // crystalball->SetParameters(tempHisto->GetMaximum(),tempHisto->GetMean(),tempHisto->GetRMS(),1,3);
                                // tempHisto->Fit(crystalball,"Q","",CTRmin,CTRmax);

                                delta_X_core.push_back(midW);
                                W_Y_core.push_back(res[0]);
                                rmsBasic_Y_core.push_back(res[1]);
                                // // ------- BEGIN OF MODS FOR AMPL CORRECTION
                                // CurrentCrystal->AddTvsQHistos(TvsQ);
                                CurrentCrystal->AddDeltaTHistos(tempHisto);
                                var.str("");
                                sname.str("");
                                // // ------- END OF MODS FOR AMPL CORRECTION

                              }

                              //add beginning and end points
                              delta_X.push_back(0.0);
                              W_Y.push_back(W_Y_core[0]);
                              rmsBasic_Y.push_back(rmsBasic_Y_core[0]);
                              for(unsigned int iCore = 0; iCore < delta_X_core.size(); iCore++)
                              {
                                delta_X.push_back(delta_X_core[iCore]);
                                W_Y.push_back(W_Y_core[iCore]);
                                rmsBasic_Y.push_back(rmsBasic_Y_core[iCore]);
                              }
                              delta_X.push_back(1.0);
                              W_Y.push_back(W_Y_core[delta_X_core.size()-1]);
                              rmsBasic_Y.push_back(rmsBasic_Y_core[delta_X_core.size()-1]);



                              TGraph *graphDeltaW = new TGraph(delta_X.size(),&delta_X[0],&W_Y[0]);
                              sname.str("");
                              sname << "DeltaW Graph - Crystal " << CurrentCrystal->GetID();
                              graphDeltaW->SetTitle(sname.str().c_str());
                              graphDeltaW->SetName(sname.str().c_str());
                              graphDeltaW->GetXaxis()->SetTitle("W");
                              graphDeltaW->GetYaxis()->SetTitle("T crystal - T tagging [S]");
                              sname.str("");
                              CurrentCrystal->SetGraphDeltaW(graphDeltaW);

                              TGraph *graphDeltaRMS = new TGraph(delta_X.size(),&delta_X[0],&rmsBasic_Y[0]);
                              sname.str("");
                              sname << "RMS DeltaW Graph - Crystal " << CurrentCrystal->GetID();
                              graphDeltaRMS->SetTitle(sname.str().c_str());
                              graphDeltaRMS->SetName(sname.str().c_str());
                              graphDeltaRMS->GetXaxis()->SetTitle("W");
                              graphDeltaRMS->GetYaxis()->SetTitle("T crystal - T tagging [S]");
                              sname.str("");
                              CurrentCrystal->SetGraphDeltaRMS(graphDeltaRMS);
                            }


                            //---------------------------------------------//
                            // Neighboring crystals analysis and w slicing //
                            //---------------------------------------------//
                            // plots for the neighbour channels (channels, NOT crystals!) of T cry - T neighbour
                            // one for each neighbour

                            std::vector<int> DelayTimingChannelsNum;
                            for(unsigned int iNeig = 0; iNeig < neighbours.size(); iNeig++)
                            {
                              //get the timingChannel and digitizerChannel from the neighbours[iNeig] value, i.e. the digitizerChannel of this iNeigh
                              int iNeighTimingChannel;
                              int iNeighDigitizerChannel;
                              for(unsigned int iDet = 0; iDet < detector.size(); iDet++)
                              {
                                if(detector[iDet].digitizerChannel == neighbours[iNeig])
                                {
                                  iNeighTimingChannel = detector[iDet].timingChannel;
                                  iNeighDigitizerChannel = detector[iDet].digitizerChannel;
                                }
                              }

                              std::cout << "Analyzing neighboring channel "
                                        << iNeig+1 << "/"
                                        << neighbours.size() << "... " << std::endl;

                              DelayTimingChannelsNum.push_back(iNeighTimingChannel);


                              // histogram of tCry - tNeighbour
                              sname <<  "T_Channel_" << neighbours[iNeig] << " - T_Crystal_" << CurrentCrystal->GetID();
                              var << "t" << iNeighTimingChannel
                                  << " - t" << detector[thisChannelID].timingChannel
                                  << " >> " << sname.str() ;
                              //noZerosCut
                              sNoZerosCut << "(t" << detector[thisChannelID].timingChannel << "!= 0) && ("
                                          << "t" << iNeighTimingChannel << "!= 0)";
                              noZerosCut = sNoZerosCut.str().c_str();
                              sNoZerosCut.str("");

                              TH1F* spectrumDeltaTcryTneig = new TH1F(sname.str().c_str(),sname.str().c_str(),DeltaTimeBins,DeltaTimeMin,DeltaTimeMax);
                              tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+noZerosCut);
                              spectrumDeltaTcryTneig->GetXaxis()->SetTitle("Time [S]");
                              spectrumDeltaTcryTneig->GetYaxis()->SetTitle("N");

                              //do a preliminary fit with gauss
                              // TCanvas *cTemp  = new TCanvas("temp","temp");
                              // TF1 *gaussDummyBasicNeigh = new TF1("gaussDummyBasicNeigh","gaus");
                              // spectrumDeltaTcryTneig->Fit(gaussDummyBasicNeigh,"QN");
                              //
                              // double fitPercMin = 5.0;
                              // double fitPercMax = 6.0;
                              //
                              // double f1min = spectrumDeltaTcryTneig->GetXaxis()->GetXmin();
                              // double f1max = spectrumDeltaTcryTneig->GetXaxis()->GetXmax();
                              // // std::cout << f1min << " " << f1max << std::endl;
                              // TF1* f1  = new TF1("f1","crystalball");
                              // f1->SetLineColor(kRed);
                              // f1->SetParameters(gaussDummyBasicNeigh->GetParameter(0),gaussDummyBasicNeigh->GetParameter(1),gaussDummyBasicNeigh->GetParameter(2),1,3);
                              // double fitMin = gaussDummyBasicNeigh->GetParameter(1)
                              //                 - fitPercMin*(gaussDummyBasicNeigh->GetParameter(2));
                              // double fitMax = gaussDummyBasicNeigh->GetParameter(1)
                              //                 + fitPercMax*(gaussDummyBasicNeigh->GetParameter(2));
                              // if(fitMin < f1min)
                              // {
                              //   fitMin = f1min;
                              // }
                              // if(fitMax > f1max)
                              // {
                              //   fitMax = f1max;
                              // }

                              double fitPercMin = 5.0;
                              double fitPercMax = 6.0;
                              int divisions = 10000;
                              double res[2];
                              extractFromCrystalBall(spectrumDeltaTcryTneig,fitPercMin,fitPercMax,divisions,res);
                              tChannelsForPolishedCorrection.push_back(iNeighTimingChannel);
                              meanForPolishedCorrection.push_back(res[0]);
                              fwhmForPolishedCorrection.push_back(res[1]);



                              // spectrumDeltaTcryTneig->Fit(f1,"Q","",fitMin,fitMax);



                              var.str("");
                              sname.str("");

                              sname << "T_Channel_" << neighbours[iNeig] << " - T_Crystal_" << CurrentCrystal->GetID() << " vs. W";
                              var << "t" << iNeighTimingChannel
                                  << " - t" << detector[thisChannelID].timingChannel
                                  << ":" << FloodZ.str() <<" >> " << sname.str() ;
                              //noZerosCut
                              sNoZerosCut << "(t" << detector[thisChannelID].timingChannel << "!= 0) && ("
                                          << "t" << iNeighTimingChannel << "!= 0)";
                              noZerosCut = sNoZerosCut.str().c_str();
                              sNoZerosCut.str("");
                              TH2F* spectrumCrystalDeltaT2vsW = new TH2F(sname.str().c_str(),
                                                                         sname.str().c_str(),
                                                                         wHistogramsBins,
                                                                         spectrumHistoW->GetMean() - 3.0*spectrumHistoW->GetRMS(),spectrumHistoW->GetMean() + 3.0*spectrumHistoW->GetRMS(),
                                                                         DeltaTimeBins,
                                                                         spectrumDeltaTcryTneig->GetMean() - 3.0*spectrumDeltaTcryTneig->GetRMS(),spectrumDeltaTcryTneig->GetMean() + 3.0*spectrumDeltaTcryTneig->GetRMS());
                              tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+noZerosCut,"COLZ");
                              spectrumCrystalDeltaT2vsW->GetXaxis()->SetTitle("W");
                              sname.str("");
                              sname << "T_channel_"<< neighbours[iNeig] << " - T_crystal " << CurrentCrystal->GetID() << ", [S]";
                              spectrumCrystalDeltaT2vsW->GetYaxis()->SetTitle(sname.str().c_str());
                              var.str("");
                              sname.str("");


                              int  neighID = 0;
                              for(unsigned int iNeigh = 0; iNeigh < neighbourChannels.size(); iNeigh++)
                              {
                                if( detector[neighbourChannels[iNeigh].detectorIndex].digitizerChannel == iNeighDigitizerChannel)
                                {
                                  neighID = iNeigh;
                                }
                              }

                              //---------------------------------------------//
                              // Save plots                                  //
                              //---------------------------------------------//
                              int plotPos = -1;
                              for(int iTimeMppc = 0 ; iTimeMppc < nmppcx ; iTimeMppc++)
                              {
                                for(int jTimeMppc = 0 ; jTimeMppc < nmppcy ; jTimeMppc++)
                                {
                                  if(mppc[(iModule*nmppcx)+iTimeMppc][(jModule*nmppcy)+jTimeMppc]->GetDigitizerChannel() == neighbours[iNeig])
                                  {
                                    plotPos = mppc[(iModule*nmppcx)+iTimeMppc][(jModule*nmppcy)+jTimeMppc]->GetCanvasPosition();
                                  }
                                }
                              }

                              // histogram of delta t versus "amplitude" (total charge deposited) - so far it makes sense only for polished crystals
                              sname << "T_Channel_" << neighbours[iNeig] << " - T_Crystal_" << CurrentCrystal->GetID() << " vs. ch" << iNeighDigitizerChannel;
                              var << "t" << iNeighTimingChannel
                                  << " - t" << detector[thisChannelID].timingChannel
                                  << ":"    << neighbourChannels[neighID].string <<" >> " << sname.str() ;
                              //noZerosCut
                              sNoZerosCut << "(t" << detector[thisChannelID].timingChannel << "!= 0) && ("
                                          << "t" << iNeighTimingChannel << "!= 0)";
                              noZerosCut = sNoZerosCut.str().c_str();
                              sNoZerosCut.str("");
                              //get the limits from the corresponding light sharing plot
                              double lsMin;
                              double lsMax;
                              for(unsigned int iLS = 0 ; iLS < CurrentCrystal->GetLSSpectra().size() ; iLS++)
                              {
                                if(CurrentCrystal->GetLSSpectra()[iLS].canvasPosition == plotPos)
                                {
                                  lsMin = CurrentCrystal->GetLSSpectra()[iLS].spectrum->GetMean() -
                                          3.0 * CurrentCrystal->GetLSSpectra()[iLS].spectrum->GetRMS();
                                  lsMax = CurrentCrystal->GetLSSpectra()[iLS].spectrum->GetMean() +
                                          3.0 * CurrentCrystal->GetLSSpectra()[iLS].spectrum->GetRMS();
                                }
                              }


                              TH2F* spectrumCrystalDeltaT2vsCH = new TH2F(sname.str().c_str(),sname.str().c_str(),
                                                                          histo1Dbins,
                                                                          lsMin,
                                                                          lsMax,
                                                                          DeltaTimeBins,
                                                                          spectrumDeltaTcryTneig->GetMean() - 3.0*spectrumDeltaTcryTneig->GetRMS(),spectrumDeltaTcryTneig->GetMean() + 3.0*spectrumDeltaTcryTneig->GetRMS());
                              tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+noZerosCut,"COLZ");
                              spectrumCrystalDeltaT2vsCH->GetXaxis()->SetTitle("W");
                              sname.str("");
                              sname << "T_channel_"<< neighbours[iNeig] << " - T_crystal " << CurrentCrystal->GetID() << ", [S]";
                              spectrumCrystalDeltaT2vsCH->GetYaxis()->SetTitle(sname.str().c_str());
                              sname.str("");
                              sname << "Ch"<< iNeighDigitizerChannel << " [ADC channels]";
                              spectrumCrystalDeltaT2vsCH->GetXaxis()->SetTitle(sname.str().c_str());
                              var.str("");
                              sname.str("");

                              std::vector<float> delay_X,delay_X_core;
                              std::vector<float> Wcoord_Y,Wcoord_Y_core;
                              std::vector<float> rms_Y,rms_Y_core;
                              TGraph *graphDelayW;
                              TGraph *graphDelayRMS;

                              if(timingCorrection)
                              {
                                //get TGraphs from building N th1f, in the range defined by beginW and endW previously found
                                // the range beginw-endW is spilt in WrangeBinsForTiming parts,
                                for(int iBin = 0; iBin < WrangeBinsForTiming; iBin++) // -1 and +1 are put to include the w limits
                                {
                                  // std::cout << "W Slice" << iBin << std::endl;
                                  var.str("");
                                  sname.str("");
                                  sname << "Delay ch_" << neighbours[iNeig] << "_t_" << iNeighTimingChannel << "_w_" << iBin;

                                  TH1F *tempHisto = new TH1F(sname.str().c_str(),sname.str().c_str(),DeltaTimeBins,DeltaTimeMin,DeltaTimeMax);
                                  var << "t" << iNeighTimingChannel
                                  << " - t" << detector[thisChannelID].timingChannel
                                  << " >> " << sname.str().c_str();

                                  //noZerosCut
                                  sNoZerosCut << "(t" << detector[thisChannelID].timingChannel << "!= 0) && ("
                                              << "t" << iNeighTimingChannel << "!= 0)";
                                  noZerosCut = sNoZerosCut.str().c_str();
                                  sNoZerosCut.str("");

                                  std::stringstream sCut;
                                  sCut << FloodZ.str() << " > "
                                  << beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming)
                                  << "&& " << FloodZ.str() << " < "
                                  << beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming)
                                  << "&& ("
                                  << "t" << iNeighTimingChannel
                                  << " - t" << detector[thisChannelID].timingChannel
                                  << ") > " << DeltaTimeMin
                                  << "&& ("
                                  << "t" << iNeighTimingChannel
                                  << " - t" << detector[thisChannelID].timingChannel
                                  << ") < " << DeltaTimeMax;
                                  TCut wCut = sCut.str().c_str();
                                  tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+wCut+noZerosCut);

                                  //do a preliminary fit with gauss
                                  // TCanvas *cTemp  = new TCanvas("temp","temp");
                                  // TF1 *gaussDummyBasicCorr = new TF1("gaussDummyBasicCorr","gaus");
                                  // tempHisto->Fit(gaussDummyBasicCorr,"QN");
                                  //
                                  // double fitPercMin = 5.0;
                                  // double fitPercMax = 6.0;
                                  //
                                  // double f1min = tempHisto->GetXaxis()->GetXmin();
                                  // double f1max = tempHisto->GetXaxis()->GetXmax();
                                  // // std::cout << f1min << " " << f1max << std::endl;
                                  // TF1* crystalball  = new TF1("f1","crystalball");
                                  // crystalball->SetLineColor(kRed);
                                  // crystalball->SetParameters(gaussDummyBasicCorr->GetParameter(0),gaussDummyBasicCorr->GetParameter(1),gaussDummyBasicCorr->GetParameter(2),1,3);
                                  // double fitMin = gaussDummyBasicCorr->GetParameter(1)
                                  //                 - fitPercMin*(gaussDummyBasicCorr->GetParameter(2));
                                  // double fitMax = gaussDummyBasicCorr->GetParameter(1)
                                  //                 + fitPercMax*(gaussDummyBasicCorr->GetParameter(2));
                                  // if(fitMin < f1min)
                                  // {
                                  //   fitMin = f1min;
                                  // }
                                  // if(fitMax > f1max)
                                  // {
                                  //   fitMax = f1max;
                                  // }
                                  // tempHisto->Fit(crystalball,"Q","",fitMin,fitMax);

                                  // TF1* crystalball  = new TF1("crystalball","crystalball");
                                  // crystalball->SetParameters(tempHisto->GetMaximum(),tempHisto->GetMean(),tempHisto->GetRMS(),1,3);
                                  // tempHisto->Fit(crystalball,"Q","",DeltaTimeMin,DeltaTimeMax);

                                  double fitPercMin = 5.0;
                                  double fitPercMax = 6.0;
                                  int divisions = 10000;
                                  double res[2];
                                  extractFromCrystalBall(tempHisto,fitPercMin,fitPercMax,divisions,res);
                                  // tChannelsForPolishedCorrection.push_back(iNeighTimingChannel);
                                  // meanForPolishedCorrection.push_back(res[0]);
                                  // fwhmForPolishedCorrection.push_back(res[1]);


                                  delay_X_core.push_back( beginW + (((iBin+0.5)*(endW - beginW))/WrangeBinsForTiming) );
                                  Wcoord_Y_core.push_back(res[0]);
                                  rms_Y_core.push_back(res[1]);
                                  CurrentCrystal->AddDelayTHistos(tempHisto);
                                  var.str("");
                                  sname.str("");
                                }

                                //add beginning and end points
                                delay_X.push_back(0.0);
                                Wcoord_Y.push_back(Wcoord_Y_core[0]);
                                rms_Y.push_back(rms_Y_core[0]);
                                for(unsigned int iCore = 0; iCore < delay_X_core.size(); iCore++)
                                {
                                  delay_X.push_back(delay_X_core[iCore]);
                                  Wcoord_Y.push_back(Wcoord_Y_core[iCore]);
                                  rms_Y.push_back(rms_Y_core[iCore]);
                                }
                                delay_X.push_back(1.0);
                                Wcoord_Y.push_back(Wcoord_Y_core[delay_X_core.size()-1]);
                                rms_Y.push_back(rms_Y_core[delay_X_core.size()-1]);

                                graphDelayW = new TGraph(delay_X.size(),&delay_X[0],&Wcoord_Y[0]);
                                sname.str("");
                                sname << "Graph Delay ch_" << neighbours[iNeig] << "_t_" << iNeighTimingChannel;
                                // sname << "DeltaW Graph - Crystal " << CurrentCrystal->GetID();
                                graphDelayW->SetTitle(sname.str().c_str());
                                graphDelayW->SetName(sname.str().c_str());
                                graphDelayW->GetXaxis()->SetTitle("W");
                                sname.str("");
                                sname << "T_channel_"<< neighbours[iNeig] << " - T_crystal " << CurrentCrystal->GetID() << ", [S]";
                                graphDelayW->GetYaxis()->SetTitle(sname.str().c_str());
                                sname.str("");

                                graphDelayRMS = new TGraph(delay_X.size(),&delay_X[0],&rms_Y[0]);
                                sname.str("");
                                sname << "RMS Graph Delay ch_" << neighbours[iNeig ]<< "_t_" << iNeighTimingChannel;
                                // sname << "DeltaW Graph - Crystal " << CurrentCrystal->GetID();
                                graphDelayRMS->SetTitle(sname.str().c_str());
                                graphDelayRMS->SetName(sname.str().c_str());
                                graphDelayRMS->GetXaxis()->SetTitle("W");
                                sname.str("");
                                sname << "RMS (T_channel_"<< neighbours[iNeig] << " - T_crystal " << CurrentCrystal->GetID() << "), [S]";
                                graphDelayRMS->GetYaxis()->SetTitle(sname.str().c_str());
                                sname.str("");
                              }




                              if(plotPos != -1)
                              {
                                CurrentCrystal->AddDeltaTcryTneig(spectrumDeltaTcryTneig,plotPos);
                                CurrentCrystal->AddDeltaT2vsW(spectrumCrystalDeltaT2vsW,plotPos);
                                CurrentCrystal->AddDeltaT2vsCH(spectrumCrystalDeltaT2vsCH,plotPos);
                                if(timingCorrection)
                                {
                                  CurrentCrystal->AddGraphDelayW(graphDelayW,plotPos);
                                  CurrentCrystal->AddGraphDelayRMS(graphDelayRMS,plotPos);
                                }
                              }
                            }
                            CurrentCrystal->SetDelayTimingChannels(DelayTimingChannelsNum);

                            CurrentCrystal->SetTChannelsForPolishedCorrection(tChannelsForPolishedCorrection);
                            CurrentCrystal->SetMeanForPolishedCorrection(meanForPolishedCorrection);
                            CurrentCrystal->SetFwhmForPolishedCorrection(fwhmForPolishedCorrection);
                          }
                        }
                      }


                      //Compton calibration plots
                      //TEST on single. MPPC is now C2, crystal 28. We want to see just MPPC C3
                      //the p_C3(E,w)
                      // create a 3x3 matrix of TGraph2D*
                      // if(comptonAnalysis)   // FIXME compton analysis is OFF for now!
                      // {
                      //   TGraph2D***  ComptonCalibation;
                      //   ComptonCalibation = new TGraph2D** [nmppcx];
                      //   for(int iCal = 0; iCal < nmppcx ; iCal++) ComptonCalibation[iCal] = new TGraph2D*[nmppcy];
                      //
                      //   TGraph2D***  ConvertedComptonCalibration;
                      //   ConvertedComptonCalibration = new TGraph2D** [nmppcx];
                      //   for(int iCal = 0; iCal < nmppcx ; iCal++) ConvertedComptonCalibration[iCal] = new TGraph2D*[nmppcy];
                      //
                      //
                      //   for(int iComptMppc = 0; iComptMppc < nmppcx ; iComptMppc++)
                      //   {
                      //     for(int jComptMppc = 0; jComptMppc < nmppcy ; jComptMppc++)
                      //     {
                      //       sname.str("");
                      //       var.str("");
                      //
                      //       // int channel = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
                      //       // std::cout << "Generating spectra for MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " ..." << std::endl;
                      //       //
                      //       // // this channel ID, i.e. the index of the array of struct detector[] that corresponds to this MPPC
                      //       // int thisChannelID;
                      //       // // get mppc channel ID in detector[]
                      //       // for(unsigned int iDet = 0; iDet < detector.size(); iDet++) // run on all detector[] entries
                      //       // {
                      //       //   if(detector[iDet].digitizerChannel == channel) // if the digitizer channel in detector[iDet] is the same as this mppc digitize channel
                      //       //   {
                      //       //     thisChannelID = iDet; // get the index of detector[] that identifies it
                      //       //   }
                      //       // }
                      //
                      //       int tempChannel = mppc[iComptMppc][jComptMppc]->GetDigitizerChannel();
                      //       int tempChannelID;
                      //       for(unsigned int iDet = 0; iDet < detector.size(); iDet++) // run on all detector[] entries
                      //       {
                      //         if(detector[iDet].digitizerChannel == tempChannel) // if the digitizer channel in detector[iDet] is the same as this mppc digitize channel
                      //         {
                      //           tempChannelID = iDet; // get the index of detector[] that identifies it
                      //         }
                      //       }
                      //       var << "(ch"
                      //           << tempChannel
                      //           << "):("
                      //           << SumChannels
                      //           << "):(" << FloodZ.str() <<")" ;
                      //       tree->Draw(var.str().c_str(),CrystalCut);
                      //       TGraph2D* tempGraph2d = new TGraph2D(tree->GetSelectedRows(),tree->GetV3(),tree->GetV2(), tree->GetV1());
                      //
                      //
                      //
                      //
                      //       // sname << "Histo_Pi(E,w)[" << iMppc + iNeighbour <<  "][" << jMppc + jNeighbour <<  "]_" << CurrentCrystal->GetID();
                      //       sname << "Pi_[" << iComptMppc <<  "][" << jComptMppc <<  "]_" << CurrentCrystal->GetID();
                      //       TH3I* ComptonCalibrationHistogram = new TH3I(
                      //         sname.str().c_str(),
                      //         sname.str().c_str(),
                      //         20,
                      //         tempGraph2d->GetXmin(),
                      //         tempGraph2d->GetXmax(),
                      //         20,
                      //         tempGraph2d->GetYmin(),
                      //         tempGraph2d->GetYmax(),
                      //         100,
                      //         tempGraph2d->GetZmin(),
                      //         tempGraph2d->GetZmax());
                      //
                      //
                      //         var.str("");
                      //
                      //         // ComptonCalibrationHistogram->SetTitle(sname.str().c_str());
                      //         // ComptonCalibrationHistogram->SetName(sname.str().c_str());
                      //         ComptonCalibrationHistogram->GetXaxis()->SetTitle("W");
                      //         ComptonCalibrationHistogram->GetYaxis()->SetTitle("Sum Charge [ADC ch]");
                      //         ComptonCalibrationHistogram->GetZaxis()->SetTitle("pi [ADC ch]");
                      //         var.str("");
                      //         var << "(ch"<< tempChannel << "):(" << SumChannels << "):(FloodZ)" << " >> " << sname.str();
                      //         tree->Draw(var.str().c_str(),CrystalCut);
                      //         // std::cout << ComptonCalibrationHistogram->GetEntries() << std::endl;
                      //         CurrentCrystal->SetComptonHistogram(ComptonCalibrationHistogram);
                      //         sname.str("");
                      //         sname << "Graph_Pi(q,w)[" << iComptMppc <<  "][" << jComptMppc <<  "]_" << CurrentCrystal->GetID();
                      //         std::vector<double> ex,ey,ez;
                      //         for(int iError = 0; iError < tree->GetSelectedRows(); iError++)
                      //         {
                      //           ex.push_back(0.05);
                      //           ey.push_back(20.0);
                      //           ez.push_back(20.0);
                      //         }
                      //         ComptonCalibation[iComptMppc][jComptMppc] = new TGraph2D(tree->GetSelectedRows(),tree->GetV3(),tree->GetV2(), tree->GetV1());
                      //         ComptonCalibation[iComptMppc][jComptMppc]->SetTitle(sname.str().c_str());
                      //         ComptonCalibation[iComptMppc][jComptMppc]->SetName(sname.str().c_str());
                      //         ComptonCalibation[iComptMppc][jComptMppc]->GetXaxis()->SetTitle("W");
                      //         ComptonCalibation[iComptMppc][jComptMppc]->GetYaxis()->SetTitle("Sum Charge [ADC ch]");
                      //         ComptonCalibation[iComptMppc][jComptMppc]->GetZaxis()->SetTitle("pi [ADC ch]");
                      //         // CurrentCrystal->SetComptonCalibration(testGraph);
                      //         var.str("");
                      //         sname.str("");
                      //
                      //
                      //         //produce
                      //
                      //         // TGraph2D* graph2d = (TGraph2D*) canv-> GetPrimitive(sstream.str().c_str()); //get the 3d graph
                      //         //
                      //         // // int Npoints = graph2d->GetN();
                      //         //
                      //         sname << "Pi(z,p_tot)[" << iComptMppc <<  "][" << jComptMppc <<  "]_" << CurrentCrystal->GetID();
                      //         ConvertedComptonCalibration[iComptMppc][jComptMppc] = new TGraph2D();
                      //         Double_t *x = ComptonCalibation[iComptMppc][jComptMppc]->GetX();
                      //         Double_t *y = ComptonCalibation[iComptMppc][jComptMppc]->GetY();
                      //         Double_t *z = ComptonCalibation[iComptMppc][jComptMppc]->GetZ();
                      //         for(int i = 0 ; i < ComptonCalibation[iComptMppc][jComptMppc]->GetN() ; i++){
                      //           ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetPoint(i,calibrationGraphChosen->Eval(x[i]),y[i],z[i]);
                      //           // ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetPointError(i,2.0,sqrt(y[i]),sqrt(z[i]));
                      //         }
                      //         ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetTitle(sname.str().c_str());
                      //         ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetName(sname.str().c_str());
                      //         ConvertedComptonCalibration[iComptMppc][jComptMppc]->GetXaxis()->SetTitle("W");
                      //         ConvertedComptonCalibration[iComptMppc][jComptMppc]->GetYaxis()->SetTitle("Sum Charge [ADC ch]");
                      //         ConvertedComptonCalibration[iComptMppc][jComptMppc]->GetZaxis()->SetTitle("pi [ADC ch]");
                      //         sname.str("");
                      //
                      //
                      //       }
                      //     }
                      //
                      //
                      //     //       }
                      //     //     }
                      //     //   }
                      //     // }
                      //     CurrentCrystal->SetComptonCalibration(ComptonCalibation);
                      //     CurrentCrystal->SetConvertedComptonCalibration(ConvertedComptonCalibration);
                      //     // CurrentCrystal->SetComptonCalibrationHistogram(ComptonCalibationHistogram);
                      //     // CurrentCrystal->SetInterpolationGraph(interpolationGraph);
                      //   }


                        if(usingRealSimData) // only if this is a sim dataset
                        {
                          long long int nPoints;
                          // a 2d plot of real vs. w, using the same cuts as before


                          int RealZBins = 100;

                          sname << "Real Z vs. W - Crystal " << CurrentCrystal->GetID();
                          var << "-(RealZ-"
                              << CurrentCrystal->GetDimensionZ()/2.0
                              << "):"
                              << FloodZ.str()
                              << " >> " << sname.str();
                          TH2F* spectrum2dSimDOIplot = new TH2F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax,RealZBins,0,CurrentCrystal->GetDimensionZ());
                          nPoints = tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut ,"COLZ");
                          spectrum2dSimDOIplot->GetXaxis()->SetTitle("W");
                          spectrum2dSimDOIplot->GetYaxis()->SetTitle("Z");
                          CurrentCrystal->SetSimDOIplot(spectrum2dSimDOIplot);
                          sname.str("");
                          var.str("");

                          // same plot but with a TGraph, to allow fitting
                          sname << "Graph Z vs. W - Crystal " << CurrentCrystal->GetID();
                          TGraph* graph = new TGraph(nPoints,tree->GetV2(),tree->GetV1()); // same but TGraph (so it can be fitted in 1D)
                          graph->SetName(sname.str().c_str());
                          graph->SetTitle(sname.str().c_str());
                          graph->GetXaxis()->SetTitle("W");
                          graph->GetYaxis()->SetTitle("Z");
                          graph->Draw("ap");
                          sname.str("");
                          // 		  TF1 *linear = new TF1("linear",  "[0]*x + [1]",0,1);
                          // 		      sname << "expfit - Crystal " << CurrentCrystal->GetID();
                          // 		      TF1 *expfit = new TF1(sname.str().c_str(),  "[0]*exp(-x/[1])",0,1);
                          // 		      // 		  linear->SetParameter(0,-100);
                          // 		      // 		  linear->SetParameter(1,50);
                          // 		      expfit->SetParameter(0,50);
                          // 		      expfit->SetParameter(1,0.1);
                          // 		      // 		  graph->SetStats(1);
                          // 		      graph->Fit(sname.str().c_str(),"Q","",0.1,0.7);
                          // 		      CurrentCrystal->SetSimFit(expfit);
                          CurrentCrystal->SetSimGraph(graph);
                          sname.str("");
                          var.str("");
                          // FitSlicesX so we can plot it together with the calibration graph
                          spectrum2dSimDOIplot->FitSlicesX(0, 1, RealZBins, 0, "QNRG5S");
                          sname << spectrum2dSimDOIplot->GetName() << "_1";
                          TH1D *spectrum2dSimDOIplot_1 = (TH1D*)gDirectory->Get(sname.str().c_str());
                          sname.str("");
                          sname << spectrum2dSimDOIplot->GetName() << "_2";
                          TH1D *spectrum2dSimDOIplot_2 = (TH1D*)gDirectory->Get(sname.str().c_str());
                          sname.str("");
                          //run on the two histograms and store the results in a TGraphErrors and in a TH1F

                          sname << "Sigma W from SIM - Crystal " << CurrentCrystal->GetID();
                          TH1F *sigmaSim = new TH1F(sname.str().c_str(),sname.str().c_str(),50,0,0.05);
                          std::vector<double> ySim,xSim,exSim,eySim;
                          for(int iSim = 1 ; iSim < spectrum2dSimDOIplot_1->GetNbinsX()  ; iSim++)
                          {
                            ySim.push_back(spectrum2dSimDOIplot_1->GetBinCenter(iSim));
                            xSim.push_back(spectrum2dSimDOIplot_1->GetBinContent(iSim));
                            eySim.push_back(0);
                            exSim.push_back(spectrum2dSimDOIplot_2->GetBinContent(iSim));
                            sigmaSim->Fill(spectrum2dSimDOIplot_2->GetBinContent(iSim));
                          }
                          TGraphErrors *gTot = new TGraphErrors(xSim.size(),&xSim[0],&ySim[0],&exSim[0],&eySim[0]);
                          sname.str("");
                          CurrentCrystal->SetSimZvsW(gTot);
                          CurrentCrystal->SetSimSigmaW(sigmaSim);






                        }
                      }



                  }
                  else
                  {
                    crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)]->SetCrystalOn(false);
                  }
                }
              }
            }
          }
        }
      }

      //set the global spectra for the module
      module[iModule][jModule]->SetFloodMap2D(spectrum2dModule);
      module[iModule][jModule]->SetFloodMap3D(spectrum3dModule);

    }
  }



  // MAIN LOOP v1
  // Loop on modules, mppcs and crystal
  //---------> QUI <--------------
  /*
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      //useful strings etc
      // First cut is a general, global cut
      // it will be composed of all cuts that will apply generally to all the events
      // For example
      // 1. A global XYZ cut that allows us to consider only the events with "reasonable" uvw coordinates
      // 2. A cut on photopeak of the tagging crystal, if there is any such crystal involved
      TCut BasicCut;
      std::stringstream CutXYZstream;
      CutXYZstream << "FloodX > " << -moduleLateralSideX << " && FloodX < " << moduleLateralSideX  << "&& FloodY > " << -moduleLateralSideY <<   " && FloodY < " << moduleLateralSideY << "&& FloodZ > " << histo3Dmin << " && FloodZ <  " << histo3Dmax;
      TCut CutXYZ = CutXYZstream.str().c_str();
      std::cout << "Generating global spectra..." << std::endl;
      TString nameModule;
      std::stringstream varModule;

      BasicCut += CutXYZ;  // add the module physical constrains the to global cut of accepted events

      //
      TCut taggingPhotopeakCut = "" ;

      // int Tagging =

      // GLOBAL SPECTRA
      // Flood histogram
      nameModule = "Flood Histogram 2D - " + module[iModule][jModule]->GetName();
      varModule << "FloodY:FloodX >> " << nameModule;
      TH2F *spectrum2dModule = new TH2F(nameModule,nameModule,histo2DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo2DglobalBins,-moduleLateralSideY,moduleLateralSideY);
      tree->Draw(varModule.str().c_str(),"","COLZ");
      spectrum2dModule->SetName(nameModule);
      spectrum2dModule->SetTitle(nameModule);
      spectrum2dModule->GetXaxis()->SetTitle("U");
      spectrum2dModule->GetYaxis()->SetTitle("V");
      module[iModule][jModule]->SetFloodMap2D(spectrum2dModule);
      varModule.str("");

      if(usingTaggingBench || calcDoiResWithCalibration || taggingForTiming)//trigger spectrum
      {

        std::stringstream tagStream;
        // tagStream << "ch" << taggingCrystalChannel; // if using original files
        tagStream << "Tagging " ; // if using tree

        TaggingCrystalSpectrum =  new TH1F("TaggingCrystalSpectrum","TaggingCrystalSpectrum",taggingCrystalBins,taggingSpectrumMin,taggingSpectrumMax);
        varModule << tagStream.str().c_str() << ">> TaggingCrystalSpectrum";
        // 	std::cout << nameModule << " ... ";
        TaggingCrystalSpectrum->SetName("TaggingCrystalSpectrum");
        TaggingCrystalSpectrum->GetXaxis()->SetTitle("ADC");
        TaggingCrystalSpectrum->GetYaxis()->SetTitle("Counts");
        TaggingCrystalSpectrum->SetTitle("Spectrum of Tagging Crystal");
        tree->Draw(varModule.str().c_str(),"");
        varModule.str("");


        //restrict the region where to look for peaks. Fix for tspectrum...
        TaggingCrystalSpectrum->GetXaxis()->SetRangeUser(taggingPeakMin,taggingPeakMax);
        //find peak in the tagging crystal
        TSpectrum *sTagCrystal;
        sTagCrystal = new TSpectrum(1);
        Int_t TagCrystalPeaksN = sTagCrystal->Search(TaggingCrystalSpectrum,1,"goff",0.5);
        Double_t *TagCrystalPeaks  = sTagCrystal->GetPositionX();
        Double_t *TagCrystalPeaksY = sTagCrystal->GetPositionY();
        // float saturationPeakFraction
        //delete s;
        // float distPeak = INFINITY;
        float maxPeak = 0.0;
        int peakID = 0;
        for (int peakCounter = 0 ; peakCounter < TagCrystalPeaksN ; peakCounter++ )
        {
          // if( fabs(CrystalPeaks[peakCounter] - 0.5*(saturationPeak[iSaturation].peakMin+saturationPeak[iSaturation].peakMax)) < distPeak)//take closest peak to the center of search range selected by the user
          if(TagCrystalPeaksY[peakCounter] > maxPeak)
          {
            // distPeak = CrystalPeaks[peakCounter];
            maxPeak = TagCrystalPeaksY[peakCounter];
            peakID = peakCounter;
          }
        }
        TF1 *gaussTag = new TF1("gaussTag", "gaus",TagCrystalPeaks[peakID] - tagFitLowerFraction*TagCrystalPeaks[peakID],TagCrystalPeaks[peakID] + tagFitUpperFraction*TagCrystalPeaks[peakID]);
        gaussTag->SetParameter(1,TagCrystalPeaks[peakID]);
        gaussTag->SetParameter(0,TagCrystalPeaksY[peakID]);
        TaggingCrystalSpectrum->Fit("gaussTag","QR","",TagCrystalPeaks[peakID] - tagFitLowerFraction*TagCrystalPeaks[peakID],TagCrystalPeaks[peakID] + tagFitUpperFraction*TagCrystalPeaks[peakID]);
        TaggingCrystalSpectrum->GetXaxis()->SetRangeUser(taggingSpectrumMin,taggingSpectrumMax);
        //define a TCut for this peak
        double tagPhotopeakMin = gaussTag->GetParameter(1) - TaggingPhotopeakSigmasMin*gaussTag->GetParameter(2);
        double tagPhotopeakMax = gaussTag->GetParameter(1) + TaggingPhotopeakSigmasMax*gaussTag->GetParameter(2);
        std::stringstream tagString;
        tagString << tagStream.str().c_str() << " > " << tagPhotopeakMin << "&& " << tagStream.str().c_str() <<" < " << tagPhotopeakMax;
        taggingPhotopeakCut = tagString.str().c_str();
        //highlighted spectrum
        TriggerSpectrumHighlight = new TH1F("TriggerSpectrumHighlight","",taggingCrystalBins,taggingSpectrumMin,taggingSpectrumMax);
        varModule.str("");
        varModule << tagStream.str().c_str() << " >> TriggerSpectrumHighlight";
        TriggerSpectrumHighlight->SetLineColor(3);
        TriggerSpectrumHighlight->SetFillColor(3);
        TriggerSpectrumHighlight->SetFillStyle(3001);
        tree->Draw(varModule.str().c_str(),taggingPhotopeakCut);
        varModule.str("");
        if(TagEdgeCalculation)
        {
          tagPeakHgEntries = TriggerSpectrumHighlight->GetEntries();
        }
        if(cuttingOnTagPhotopeak) // only if tagging on photopeak is chosen
        {
          BasicCut += taggingPhotopeakCut;  // add the cut on photopeak of the tagging crystal.
        }
        //but always save the tagging cut in the output file
        taggingPhotopeakCut.SetName("taggingPhotopeakCut");
        module[iModule][jModule]->SetTaggingPhotopeakCut(taggingPhotopeakCut);        // energy cut, events in the photopeak
        // 	delete sTagCrystal;
        delete gaussTag;


        // 	std::cout << " done" << std::endl;

      }




      nameModule = "Flood Histogram 3D - Module " + module[iModule][jModule]->GetName();
      varModule << "FloodZ:FloodY:FloodX >> " << nameModule;
      //       std::cout << nameModule << " ... ";
      TH3I* spectrum3dModule = new TH3I(nameModule,nameModule,histo3DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo3DglobalBins,-moduleLateralSideY,moduleLateralSideY,histo3DglobalBins,histo3Dmin,histo3Dmax);
      // if(usingTaggingBench | calcDoiResWithCalibration)
      // {
      tree->Draw(varModule.str().c_str(),BasicCut);
      // }
      // else{
        // tree->Draw(varModule.str().c_str(),BasicCut);
      // }
      spectrum3dModule->GetXaxis()->SetTitle("U");
      spectrum3dModule->GetYaxis()->SetTitle("V");
      spectrum3dModule->GetZaxis()->SetTitle("W");
      module[iModule][jModule]->SetFloodMap3D(spectrum3dModule);
      varModule.str("");

      if(usingRealSimData)
      {
        // GLOBAL SPECTRA but with events confined in one crystal (from sim data)
        // Flood histogram
        nameModule = "SIM - Crystal Hit = 1 - Flood Histogram 2D - " + module[iModule][jModule]->GetName();
        varModule << "FloodY:FloodX >> " << nameModule;
        //       std::cout << nameModule << " ... ";
        TH2F *spectrum2dModuleSingleCrystalHit = new TH2F(nameModule,nameModule,histo2DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo2DglobalBins,-moduleLateralSideY,moduleLateralSideY);
        tree->Draw(varModule.str().c_str(),SingleCrystalInteraction,"COLZ");
        spectrum2dModuleSingleCrystalHit->SetName(nameModule);
        spectrum2dModuleSingleCrystalHit->SetTitle(nameModule);
        spectrum2dModuleSingleCrystalHit->GetXaxis()->SetTitle("U");
        spectrum2dModuleSingleCrystalHit->GetYaxis()->SetTitle("V");
        module[iModule][jModule]->SetFloodMap2DSingleCrystalHit(spectrum2dModuleSingleCrystalHit);
        varModule.str("");


      }

      //spectra for each mppc
      // #pragma omp parallel for
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
        for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
        {

          //go on only if the MPPC is on for modular analysis
          if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForModular())
          {

            // different variables logic, trying to use directly the original files:
            // even global module plots need to be produced channel by channel (otherwise the
            // definition of neighbours is impossible)
            // so the variables in the original files are just:
            //
            // - BoardTimeTags ... ch0 ... chN , t0 ... tN
            //
            // and with the pre-processing we were genetating some additional variables:
            //
            // - TriggerChannel   = the channel with highest charge among all channels in the event
            // - FloodX,Y         = weighted average on charge of detector positions, with multiple coices (i.e. summing on all detectors or only on trigger channel + 8 neighbours)
            // - FloodZ           = ratio (charge in trigger channel)/(sum charge) where sum chage can be the sum of all channels or the sum of trigger channel + 8 neighbours
            // - BadEvent         = counter of events with charge higher of saturation parameter for that specific channel. these events are "bad" because given the logarithmic formula of saturation correction, it would be mathematically impossible to correct them. It's just a counter that provides a feedback to user (typically if all makes sense, badEvents need to be just a small fraction of all events)
            // - Tagging          = charge measured in the external tagging crystal, if any
            // - TaggingTimeStamp = time stamp of the external tagging crystal, if any
            // - ZPosition        = z position of the tagging bench, for DOI measurements, for this event
            //
            // These variables can be computed on the fly directly on the original data (with some limitations, see later):
            //
            // TriggerChannel     = "max(ch0,max(ch1,...max(chN))) == chI", for I == mppc channel and 0...N are all the channels in this module.


            // write down the channel collection

            // std::stringstream TriggerChannel;

            // TriggerChannel

            // 	    mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->PrintSpecific();
            // strings and stuff
            TString name;
            std::stringstream var,cut,sname;
            //FIXME careful, at the moment it works only because there's one module. but honestly, at this stage it is supposed to work only on one module.
            // it should be fixed for more modules by using the same mppc[][] logic used for the crystals, below
            // did i already fix it? i guess i'll test once we have to deal with more tha 1 module
            int channel = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
            std::cout << "Generating spectra for MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " ..." << std::endl;
            cut << "TriggerChannel == " << channel  ;
            TCut CutTrigger = cut.str().c_str();
            cut.str("");
            // same as the global ones, but selecting on TriggerChannel
            // raw spectrum
            name = "Raw Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() + " - Module " + module[iModule][jModule]->GetName();
            var << "ch" << channel << " >> " << name;
            TH1F* spectrumRaw = new TH1F(name,name,histo1Dbins,1,histo1Dmax);
            tree->Draw(var.str().c_str(),BasicCut);
            spectrumRaw->GetXaxis()->SetTitle("ADC Channels");
            spectrumRaw->GetYaxis()->SetTitle("N");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetRawSpectrum(spectrumRaw);
            var.str("");

            //trigger selected spectrum
            name = "Trigger Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            var << "ch" << channel << " >> " << name;
            TH1F* spectrumTrigger = new TH1F(name,name,histo1Dbins,1,histo1Dmax);
            tree->Draw(var.str().c_str(),CutTrigger+BasicCut);
            spectrumTrigger->GetXaxis()->SetTitle("ADC Channels");
            spectrumTrigger->GetYaxis()->SetTitle("N");

            //set a very broad cut on the trigger spectrum (single channel, not sum) to get rid of low energy events
            std::stringstream broadCutstream;
            if(!backgroundRun)
            {
              TSpectrum *sTrigger;
              sTrigger = new TSpectrum(20);
              Int_t TriggerCrystalPeaksN    = sTrigger->Search(spectrumTrigger,2,"goff",0.2);
              Double_t *TriggerCrystalPeaks  = sTrigger->GetPositionX();
              // Float_t *TriggerCrystalPeaksY = sTrigger->GetPositionY();
              //delete s;
              float TriggermaxPeak = 0.0;
              int TriggerpeakID    = 0;
              for (int TriggerpeakCounter = 0 ; TriggerpeakCounter < TriggerCrystalPeaksN ; TriggerpeakCounter++ )
              {
                if(TriggerCrystalPeaks[TriggerpeakCounter] > TriggermaxPeak)
                {
                  TriggermaxPeak = TriggerCrystalPeaks[TriggerpeakCounter];
                  TriggerpeakID = TriggerpeakCounter;
                }
              }
              // now if 511KeV or 662KeV correspond to TriggerCrystalPeaks[TriggerpeakID], it means that a broad cut, energy > 200-250KeV is approximately that divided by 2.5 (assuming 0 is 0 and scale is linear)

              broadCutstream << "ch" << channel << ">" << (TriggerCrystalPeaks[TriggerpeakID] / ((1000*sourceMeV)/thresholdKev) );
            }
            else //otherwise, for backgroundRun set the broadcut to the value chosen by the user (or default to 0 adc channels)
            {
              broadCutstream << "ch" << channel << ">" << userBroadCut;
            }
            TCut broadCut = broadCutstream.str().c_str();
            // to make things easier, we add this directly to the CutTrigger. FIXME Not a clean solution, but ehi...
            CutTrigger += broadCut;
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTriggerSpectrum(spectrumTrigger);
            var.str("");

            //prepare an highlighted plot to show the broad cut
            name = "Trigger Spectrum Hg - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            var << "ch" << channel << " >> " << name;
            TH1F* spectrumTriggerHighlighted = new TH1F(name,name,histo1Dbins,1,histo1Dmax);
            tree->Draw(var.str().c_str(),CutTrigger+BasicCut);
            spectrumTriggerHighlighted->GetXaxis()->SetTitle("ADC Channels");
            spectrumTriggerHighlighted->GetYaxis()->SetTitle("N");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTriggerSpectrumHighlighted(spectrumTriggerHighlighted);
            var.str("");


            //trigger spectrum of this mppc, expressed in charge, WITHOUT any consideration on where the crystal of interaction.
            //useful for getting the peaks of lu-176 emitted from OTHER crystals, that interact in crystal(s) coupled to this detector.
            //it's the same as the trigger spectrum done before, but expressed in charge
            // sname << "Trigger Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            name = "Trigger Charge Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            var << "ch" << channel << "*( "<< chargeBinningADC <<  ") >> " << name;
            TH1F* spectrumTriggerCharge = new TH1F(name,name,histoSingleChargeBin,0,histoSingleChargeMax);
            tree->Draw(var.str().c_str(),BasicCut+CutTrigger);
            spectrumTriggerCharge->GetXaxis()->SetTitle("Charge [C]");
            spectrumTriggerCharge->GetYaxis()->SetTitle("N");
            // sname.str("");
            var.str("");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTriggerChargeSpectrum(spectrumTriggerCharge);

            //even more basic, the raw spectrum of this channel, expressed in terms of charge
            name = "Raw Charge Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            // sname << "Raw Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            var << "ch" << channel << "*( "<< chargeBinningADC <<  ") >> " << name;
            TH1F* spectrumRawCharge = new TH1F(name,name,histoSingleChargeBin,0,histoSingleChargeMax);
            tree->Draw(var.str().c_str(),BasicCut);
            spectrumRawCharge->GetXaxis()->SetTitle("Charge [C]");
            spectrumRawCharge->GetYaxis()->SetTitle("N");
            sname.str("");
            var.str("");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetRawChargeSpectrum(spectrumRawCharge);



            // if(digitizerType == 1 || digitizerType == 2)
            // {
            //   //global time spectrum of this mppc
            //   name = "Time Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            //   var << "t" << channel << " >> " << name;
            //   TH1F* spectrumTime = new TH1F(name,name,100,0,20e-6);
            //   tree->Draw(var.str().c_str(),BasicCut+CutTrigger);
            //   spectrumTime->GetXaxis()->SetTitle("Time [S]");
            //   spectrumTime->GetYaxis()->SetTitle("N");
            //   // sname.str("");
            //   var.str("");
            //   mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTimeSpectrum(spectrumTime);
            //
            //   // //delta time spectrum of this mppc
            //   // name = "Delta Time WRT Tagging - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            //   // var << "(t" << channel << " - TaggingTimeStamp) >> " << name;
            //   // TH1F* spectrumDeltaTimeWRTTagging = new TH1F(name,name,500,-5e-9,5e-9);
            //   // tree->Draw(var.str().c_str(),BasicCut+CutTrigger);
            //   // spectrumDeltaTimeWRTTagging->GetXaxis()->SetTitle("Time [S]");
            //   // spectrumDeltaTimeWRTTagging->GetYaxis()->SetTitle("N");
            //   // // sname.str("");
            //   // var.str("");
            //   // mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetDeltaTimeWRTTagging(spectrumDeltaTimeWRTTagging);
            // }
            //standard 2d plot
            name = "Flood Histogram 2D - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            TH2F* spectrum2dMPPC = new TH2F(name,name,histo2DchannelBin,-moduleLateralSideX,moduleLateralSideX,histo2DchannelBin,-moduleLateralSideY,moduleLateralSideY);
            var << "FloodY:FloodX >> " << name;
            tree->Draw(var.str().c_str(),BasicCut+CutTrigger,"COLZ");
            spectrum2dMPPC->GetXaxis()->SetTitle("U");
            spectrum2dMPPC->GetYaxis()->SetTitle("V");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetFloodMap2D(spectrum2dMPPC);
            var.str("");

            // 3D spectrum for this mppc
            // little trick to try and use less bins
            // take the mean x and y and their sigma from previous 2dplot, define the limit of this 3dplot in x and y accordingly
            double minX3Dplot = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetMean(1) - 3.0*mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetRMS(1);
            double maxX3Dplot = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetMean(1) + 3.0*mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetRMS(1);
            double minY3Dplot = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetMean(2) - 3.0*mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetRMS(2);
            double maxY3Dplot = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetMean(2) + 3.0*mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetRMS(2);

            //DEBUG
            // 	  std::cout << "###### Main Program " << std::endl;
            // 	  std::cout << minX3Dplot << " " << maxX3Dplot << " " << minY3Dplot << " "<< maxY3Dplot << std::endl;
            // 	  std::cout << "----------------- " << std::endl;
            name = "Flood Histogram 3D - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            var << "FloodZ:FloodY:FloodX >> " << name;

            int histo3DchannelBin =  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetHisto3DchannelBin();
            TH3I* spectrum3dMPPC = new TH3I(name,name,histo3DchannelBin,minX3Dplot,maxX3Dplot,histo3DchannelBin,minY3Dplot,maxY3Dplot,histo3DchannelBin,histo3Dmin,histo3Dmax);//FIXME temp
            // if(usingTaggingBench | calcDoiResWithCalibration)
            // {
            tree->Draw(var.str().c_str(),BasicCut+CutTrigger);
            // }
            // else
            // {
              // tree->Draw(var.str().c_str(),BasicCut+CutTrigger);
            // }
            spectrum3dMPPC->GetXaxis()->SetTitle("U");
            spectrum3dMPPC->GetYaxis()->SetTitle("V");
            spectrum3dMPPC->GetZaxis()->SetTitle("W");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetFloodMap3D(spectrum3dMPPC);
            var.str("");
            // automatic crystal finder on the mppc
            // it runs always, unless the user has set onlyuserinput //FIXME onlyuserinput is not used anymore!
            TCutG**** cutg; // prepare the graphical cuts
            //const int numbOfCrystals = 4;
            cutg = new TCutG***[2]; // two planes of cuts, their intersection will create a 3d cut
            int right_ncrystalsx;
            if(usingTaggingBench)
            { //with tagging crystal setup, use only one row of crystals
              right_ncrystalsx =1;
            }
            else
            {
              right_ncrystalsx =ncrystalsx;
            }
            for(int iCut =0 ; iCut < 2 ; iCut++)
            {
              cutg[iCut] = new TCutG**[right_ncrystalsx];
              for(int iCry = 0; iCry < right_ncrystalsx ; iCry++)
              {
                cutg[iCut][iCry] = new TCutG*[ncrystalsy];
              }
            }
            bool found = false;
            if(usingTaggingBench) //if it's a tagging bench run, check first if the mppc is on for DOI bench measurement
            {
              if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForDoi())
              {
                found = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->FindCrystalCuts(cutg,1,ncrystalsy);
              }
            }
            else
            {
              found = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->FindCrystalCuts(cutg,ncrystalsx,ncrystalsy);
            }

            // run on all the possible crystals (i.e. all the crystals coupled to this mppc)
            for(int iCry = 0; iCry < right_ncrystalsx ; iCry++)
            {

              for(int jCry = 0; jCry < ncrystalsy ; jCry++)
              {
                if(crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)]->GetIsOnForModular())
                {
                  if(found)
                  {

                    // get a pointer to this crystal
                    Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
                    CurrentCrystal->SetCrystalOn(true);
                    //store the cutg in the crystal
                    CurrentCrystal->SetZXCut(cutg[0][iCry][jCry]);
                    CurrentCrystal->SetZYCut(cutg[1][iCry][jCry]);

                    // create "sum channels" string, a string to have the variable "sum of all channels" to be used later
                    // it will be either all channels selected as input
                    // or just trigger channel plus neighbourChannels

                    // first get the input digitizer channels
                    std::vector<std::string> digitizer_f;
                    config.split( digitizer_f, digitizer_s, "," );
                    std::vector<int> digitizer;
                    for(unsigned int i = 0 ; i < digitizer_f.size() ; i++)
                    {
                      config.trim(digitizer_f[i]);
                      digitizer.push_back(atoi(digitizer_f[i].c_str()));
                    }
                    // create the string
                    std::stringstream sSumChannels;
                    std::string SumChannels;
                    std::vector<int> neighbours = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetNeighbours(); //get neighbours anyway, they will be useful later
                    // now the 2 cases

                    if(usingAllChannelsForEnergySpectra) // add all channels in the digitizer list
                    {
                      sSumChannels << "ch" <<  digitizer[0];
                      for(unsigned int iDigi = 1 ; iDigi < digitizer.size() ; iDigi++)
                      {
                        sSumChannels << "+ch" << digitizer[iDigi];
                      }
                      SumChannels = sSumChannels.str();
                      sSumChannels.str("");
                    }
                    else // add only this channel and the neighbourChannels
                    {
                      sSumChannels << "ch" << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();

                      for(unsigned int iNeighbour = 0 ; iNeighbour < neighbours.size() ; iNeighbour++)
                      {
                        sSumChannels << "+ch" << neighbours[iNeighbour];
                      }
                      SumChannels = sSumChannels.str();
                      sSumChannels.str("");
                    }


                    // std::cout << SumChannels << std::endl;
                    //----------------------------------------------------------//

                    std::cout << "Generating spectra for crystal " << CurrentCrystal->GetID() << " ..." << std::endl;

                    // a cut for events falling in this crystal. This basically means
                    // 1. "Physical" constrains are satisfied (u and v are not beyond the broad physical x-y limits of the array, w is between 0 and 1). This is set at module level
                    // 2. If a tagging crystal is used, events have charge in the channel of tagging crystal in the photopeak of the tagging crystal. This is set at module level
                    // 3. The MPPC where the crystal is has the maximum signal in the event ("trigger" channel).
                    // 4. The "broad" cut on energy is applied, which means that the user decided to discard events that in the "trigger" spectrum of this mppc are "too low". This is set at mppc level
                    // 5. Events are in the density area of this crystal, defined by the two tcutg
                    //
                    // BasicCut = 1. and 2.
                    // CutTrigger = 3. and 4.
                    // CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName() = 5.

                    TCut CrystalCut = BasicCut + CutTrigger + CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName();

                    // a 3d historgram for this crystal, mainly to check the 3d cut
                    sname << "Flood Histogram 3D - Crystal " << CurrentCrystal->GetID();
                    var << "FloodZ:FloodY:FloodX >> " << sname.str();
                    // NB histo3DchannelBin is taken from the MPPC element earlier in the code
                    TH3I* spectrum3dCrystal = new TH3I(sname.str().c_str(),sname.str().c_str(),histo3DchannelBin,minX3Dplot,maxX3Dplot,histo3DchannelBin,minY3Dplot,maxY3Dplot,histo3DchannelBin,histo3Dmin,histo3Dmax);
                    tree->Draw(var.str().c_str(),CrystalCut);
                    spectrum3dCrystal->GetXaxis()->SetTitle("U");
                    spectrum3dCrystal->GetYaxis()->SetTitle("V");
                    spectrum3dCrystal->GetZaxis()->SetTitle("W");
                    CurrentCrystal->SetFloodMap3D(spectrum3dCrystal);
                    sname.str("");
                    var.str("");

                    // a 2d historgram for this crystal, mainly to check the 3d cut
                    sname << "Flood Histogram 2D - Crystal " << CurrentCrystal->GetID();
                    var << "FloodY:FloodX >> " << sname.str();
                    TH2F* spectrum2dCrystal = new TH2F(sname.str().c_str(),sname.str().c_str(),histo2DchannelBin,-moduleLateralSideX,moduleLateralSideX,histo2DchannelBin,-moduleLateralSideY,moduleLateralSideY);
                    tree->Draw(var.str().c_str(),CrystalCut);
                    spectrum2dCrystal->GetXaxis()->SetTitle("U");
                    spectrum2dCrystal->GetYaxis()->SetTitle("V");
                    CurrentCrystal->SetFloodMap2D(spectrum2dCrystal);
                    sname.str("");
                    var.str("");


                    if(saturationRun) //only if this is a saturation analysis run
                    {
                      //sum spectrum for events localized in this crystal, but expressed in Q
                      sname << "Sum Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      // var << SumChannels << " >> " << sname.str();
                      var << "(" << SumChannels << ")*( "<< chargeBinningADC <<  ") >> " << sname.str();
                      TH1F* spectrumSumCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histoSumChargeBin,0,histoSumChargeMax);
                      tree->Draw(var.str().c_str(),CrystalCut);
                      spectrumSumCharge->GetXaxis()->SetTitle("Charge [C]");
                      spectrumSumCharge->GetYaxis()->SetTitle("N");
                      sname.str("");
                      var.str("");
                      CurrentCrystal->SetSumChargeSpectrum(spectrumSumCharge);
                      // CurrentCrystal->SetSaturationFits(gaussFitSaturation);


                      //sum of all channels minus the trigger channel, for events confined in this crystal, expressed in Q
                      //this should be almost not saturated. useful for sources, to compare with background spectrum and get an energy measurement of relevant parts (to use in the
                      // identification of relevant part of bg spectrum in single channels spectrum of bg)
                      sname << "All-Minus-Trigger Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      // var << SumChannels << " >> " << sname.str();
                      var << "(" << SumChannels << "-ch" << channel << ")*( "<< chargeBinningADC <<  ") >> " << sname.str();
                      TH1F* spectrumAMTCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histoSumChargeBin,0,histoSumChargeMax);
                      tree->Draw(var.str().c_str(),CrystalCut);
                      spectrumAMTCharge->GetXaxis()->SetTitle("Charge [C]");
                      spectrumAMTCharge->GetYaxis()->SetTitle("N");
                      sname.str("");
                      var.str("");
                      CurrentCrystal->SetAMTChargeSpectrum(spectrumAMTCharge);
                      // CurrentCrystal->SetSaturationFits(gaussFitSaturation);


                      //single MPPC spectrum for this crystal. for the saturation Run the important quantity is the q_max, i.e. th maximum charge that can be
                      //"seen" by the SiPM, corresponding to a situation where all the pixels are firing. So these plots are generate not in ADCch but in charge, using the chargeBinningADC
                      sname << "Single Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      var << "ch" << channel << "*( "<< chargeBinningADC <<  ") >> " << sname.str();
                      TH1F* spectrumSingleCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histoSingleChargeBin,0,histoSingleChargeMax);
                      tree->Draw(var.str().c_str(),CrystalCut);
                      spectrumSingleCharge->GetXaxis()->SetTitle("Charge [C]");
                      spectrumSingleCharge->GetYaxis()->SetTitle("N");
                      sname.str("");
                      var.str("");


                      //spectrum of charge in this mppc, but when the trigger is NOT itself nor a neighbouring channel
                      sname << "Not-Neighbour Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      // std::vector<int> neighbours = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetNeighbours();
                      TCut NotNeighbours;
                      std::stringstream sneigh;
                      sneigh << "TriggerChannel != " <<   mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
                      NotNeighbours += sneigh.str().c_str();
                      for(unsigned int iNeig = 0 ;  iNeig < neighbours.size(); iNeig++)
                      {
                        std::stringstream sneigh2;
                        sneigh2.str("");
                        sneigh2 << "TriggerChannel != " <<   neighbours[iNeig];
                        NotNeighbours += sneigh2.str().c_str();
                      }
                      var << "ch" << channel << "*( "<< chargeBinningADC <<  ") >> " << sname.str();
                      TH1F* NotNeighboursSingleCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histoSingleChargeBin,0,histoSingleChargeMax);
                      // std::cout << NotNeighbours << std::endl;
                      tree->Draw(var.str().c_str(),BasicCut+NotNeighbours);
                      NotNeighboursSingleCharge->GetXaxis()->SetTitle("Charge [C]");
                      NotNeighboursSingleCharge->GetYaxis()->SetTitle("N");

                      sname.str("");
                      var.str("");
                      CurrentCrystal->SetNotNeighboursSingleCharge(NotNeighboursSingleCharge);




                      //look for peaks in the areas selected by the user
                      //here there could be more peaks the user is interested in (for example, in na22)
                      //and they could be anywhere, so it's pointless to implement a very complicated algorithm
                      //to find all possible peaks...
                      std::vector<TF1*> gaussFitSaturation;
                      TH1F* spectrumWhereToSearch;  //the spectrum where to search peaks is different if the run is backgroun saturation or not.

                      if(backgroundSaturationRun)
                      {
                        spectrumWhereToSearch = (TH1F*) NotNeighboursSingleCharge->Clone();
                      }
                      else
                      {
                        spectrumWhereToSearch = (TH1F*) spectrumSingleCharge->Clone();
                      }
                      // in background saturation runs, in fact, the spectrum where to search is the trigger spectrum without any crystal consideration, where the 202 Kev and 307 Kev peaks
                      // are much more clear. otherwise the single channel spectrum cut on the events in the crystal is used
                      if(performSaturationPeakSearch)
                      {
                        spectrumWhereToSearch->SetName("Investigated Spectrum");
                        float maxPeakHight = 0;
                        for(unsigned int iSaturation = 0 ; iSaturation < saturationPeak.size(); iSaturation++)
                        {
                          sname << "Peak " << saturationPeak[iSaturation].energy << " KeV - Crystal " << CurrentCrystal->GetID() << " - MPPC " <<  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();

                          spectrumWhereToSearch->GetXaxis()->SetRangeUser(saturationPeak[iSaturation].peakMin,saturationPeak[iSaturation].peakMax);

                          TSpectrum *s;
                          s = new TSpectrum(20);
                          // 		Input[i].SumSpectraCanvas->cd(j+1);
                          Int_t CrystalPeaksN = s->Search(spectrumWhereToSearch,1,"nobackground goff",0.3);
                          Double_t *CrystalPeaks = s->GetPositionX();
                          Double_t *CrystalPeaksY = s->GetPositionY();
                          // float saturationPeakFraction
                          //delete s;
                          // float distPeak = INFINITY;
                          float maxPeak = 0.0;
                          int peakID = 0;
                          for (int peakCounter = 0 ; peakCounter < CrystalPeaksN ; peakCounter++ )
                          {
                            // if( fabs(CrystalPeaks[peakCounter] - 0.5*(saturationPeak[iSaturation].peakMin+saturationPeak[iSaturation].peakMax)) < distPeak)//take closest peak to the center of search range selected by the user
                            if(CrystalPeaksY[peakCounter] > maxPeak)
                            {
                              // distPeak = CrystalPeaks[peakCounter];
                              maxPeak = CrystalPeaksY[peakCounter];
                              peakID = peakCounter;
                            }
                          }


                          TF1 *satGauss = new TF1(sname.str().c_str(),  "gaus",CrystalPeaks[peakID]-CrystalPeaks[peakID]*saturationPeakFractionLow,CrystalPeaks[peakID]+CrystalPeaks[peakID]*saturationPeakFractionHigh);
                          satGauss->SetParameter(1,CrystalPeaks[peakID]);
                          satGauss->SetParameter(0,CrystalPeaksY[peakID]);
                          spectrumWhereToSearch->Fit(sname.str().c_str(),"QN","",CrystalPeaks[peakID]-CrystalPeaks[peakID]*saturationPeakFractionLow,CrystalPeaks[peakID]+CrystalPeaks[peakID]*saturationPeakFractionHigh);
                          gaussFitSaturation.push_back(satGauss);
                          saturationFile << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << "\t"
                                         << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetI() << "\t"
                                         << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetJ() << "\t"
                                         << CurrentCrystal->GetID() << "\t"
                                         << CurrentCrystal->GetI() << "\t"
                                         << CurrentCrystal->GetJ() << "\t"
                                         << saturationPeak[iSaturation].energy << "\t"
                                         << satGauss->GetParameter(1) << "\t"
                                         << satGauss->GetParError(1)
                                         << std::endl;
                          spectrumWhereToSearch->GetXaxis()->SetRangeUser(0,histoSingleChargeMax);
                          if(CrystalPeaksY[peakID] > maxPeakHight )
                          {
                            maxPeakHight = CrystalPeaksY[peakID];
                          }

                          sname.str("");
                        }
                        spectrumWhereToSearch->GetYaxis()->SetRangeUser(0,maxPeakHight*2.0);
                        CurrentCrystal->SetSaturationFits(gaussFitSaturation);
                      }
                      CurrentCrystal->SetInvestigatedSpectrum(spectrumWhereToSearch);
                      CurrentCrystal->SetSingleChargeSpectrum(spectrumSingleCharge);

                    }
                    else
                    {

                      //-------------------------------------------------------------------------
                      //standard sum spectrum with cut on crystal events, xyz and trigger channel
                      //------------------------------------------------------------------------
                      //draw charge spectrum
                      sname << "Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      var << SumChannels << " >> " << sname.str();
                      TH1F* spectrumCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut);
                      spectrumCharge->GetXaxis()->SetTitle("ADC Channels");
                      spectrumCharge->GetYaxis()->SetTitle("N");
                      sname.str("");
                      var.str("");
                      if(lightYieldComputation)
                      {
                        //SumSpectrum in Ph/MeV -- CAREFUL this is not as accurate as measuring LY on PMTs
                        sname << "Light Yield @ "<< sourceMeV <<  " MeV - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                        var << "(" << SumChannels << ")* " << chargeBinningADC << "/(" << gainMPPC <<"*1.6e-19* " << qe << "*"<<  sourceMeV << ") >> " << sname.str();
                        TH1F* spectrumLY = new TH1F(sname.str().c_str(),sname.str().c_str(),histoLYbins,1,histoLYmax);
                        tree->Draw(var.str().c_str(),CrystalCut);
                        spectrumLY->GetXaxis()->SetTitle("[Ph/MeV]");
                        spectrumLY->GetYaxis()->SetTitle("N");
                        CurrentCrystal->SetLYSpectrum(spectrumLY);
                        sname.str("");
                        var.str("");

                        TSpectrum *s;
                        s = new TSpectrum(20);
                        // 		Input[i].SumSpectraCanvas->cd(j+1);
                        Int_t CrystalPeaksN = s->Search(spectrumLY,2,"goff",0.5);
                        Double_t *CrystalPeaks = s->GetPositionX();
                        Double_t *CrystalPeaksY = s->GetPositionY();
                        //delete s;
                        float maxPeak = 0.0;
                        int peakID = 0;
                        for (int peakCounter = 0 ; peakCounter < CrystalPeaksN ; peakCounter++ )
                        {
                          if(CrystalPeaks[peakCounter] > maxPeak)
                          {
                            maxPeak = CrystalPeaks[peakCounter];
                            peakID = peakCounter;
                          }
                        }
                        //fit the spectra - TODO use the gaussian plus fermi?
                        if (energyResolution == 0)
                        {
                          if (correctingSaturation)
                          energyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
                          else
                          energyResolution = ENERGY_RESOLUTION;
                        }
                        float par0 = CrystalPeaksY[peakID];
                        float par1 = CrystalPeaks[peakID];
                        float par2 = (CrystalPeaks[peakID]*energyResolution)/2.35;
                        float fitmin = par1-1.5*par2;
                        float fitmax = par1+1.8*par2;
                        sname << "gaussLY - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                        TF1 *gauss = new TF1(sname.str().c_str(),  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitmin,fitmax);
                        gauss->SetParameter(0,par0);
                        gauss->SetParameter(1,par1);
                        gauss->SetParameter(2,par2);
                        spectrumLY->Fit(sname.str().c_str(),"Q","",fitmin,fitmax);
                        //store the mean and sigma in the crystal
                        if(gauss->GetParameter(1) > 0) // otherwise the fit was very wrong..)
                        CurrentCrystal->SetLY(gauss->GetParameter(1),std::abs(gauss->GetParameter(2)));
                        CurrentCrystal->SetLYFit(gauss);


                      }
                      //photons*1.25e6*1.6e-19/156e-15

                      //prepare the photopeak cuts and stuff
                      TCut PhotopeakEnergyCutCorrected = "";
                      TCut PhotopeakEnergyCut  = "";
                      double EnergyCutMin;
                      double EnergyCutMax;
                      int bin3,bin4;
                      double wbin3,wbin4,meanW20;
                      //automatically look for the 511Kev peak to find the photopeak energy cut
                      //find peaks in each crystal spectrum, with TSpectrum
                      if(!backgroundRun)// do it only if this is NOT a background run
                      {
                        TSpectrum *s;
                        s = new TSpectrum(20);
                        // 		Input[i].SumSpectraCanvas->cd(j+1);
                        Int_t CrystalPeaksN = s->Search(spectrumCharge,2,"",0.5);
                        Double_t *CrystalPeaks =  s->GetPositionX();
                        Double_t *CrystalPeaksY = s->GetPositionY();
                        //delete s;
                        float maxPeak = 0.0;
                        int peakID = 0;
                        for (int peakCounter = 0 ; peakCounter < CrystalPeaksN ; peakCounter++ )
                        {
                          if(CrystalPeaks[peakCounter] > maxPeak)
                          {
                            maxPeak = CrystalPeaks[peakCounter];
                            peakID = peakCounter;
                          }
                        }
                        //fit the spectra - TODO use the gaussian plus fermi?
                        if (energyResolution == 0)
                        {
                          if (correctingSaturation)
                          energyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
                          else
                          energyResolution = ENERGY_RESOLUTION;
                        }
                        Float_t par0 = CrystalPeaksY[peakID];
                        Float_t par1 = CrystalPeaks[peakID];
                        Float_t par2 = (CrystalPeaks[peakID]*energyResolution)/2.35;
                        Float_t fitmin = par1-1.5*par2;
                        Float_t fitmax = par1+1.8*par2;
                        // std::cout << par0   << "\t"
                        //           << par1   << "\t"
                        //           << par2   << "\t"
                        //           << fitmin << "\t"
                        //           << fitmax << std::endl;

                        sname << "gaussCharge - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                        TF1 *gauss = new TF1(sname.str().c_str(),  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitmin,fitmax);
                        gauss->SetParameter(0,par0);
                        gauss->SetParameter(1,par1);
                        gauss->SetParameter(2,par2);
                        spectrumCharge->Fit(sname.str().c_str(),"Q","",fitmin,fitmax);

                        if(usingTaggingBench && TagEdgeCalculation)
                        {
                          std::cout << "Crystal "
                                    << CurrentCrystal->GetID()
                                    << " - MPPC "
                                    << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel()
                                    << ", Ratio = "
                                    << fabs(gauss->GetParameter(0)*
                                            gauss->GetParameter(2)*
                                            TMath::Sqrt(2.0*TMath::Pi())) /
                                            tagPeakHgEntries
                                    << std::endl;
                        }
                        //store the mean and sigma in the crystal
                        if(gauss->GetParameter(1) > 0) // otherwise the fit was very wrong..)
                        CurrentCrystal->SetPhotopeak(gauss->GetParameter(1),std::abs(gauss->GetParameter(2)));
                        CurrentCrystal->SetFit(gauss);
                        // 		std::cout << "Photopeak Mean for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakPosition() << std::endl;
                        // 		std::cout << "Photopeak Sigma for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakSigma() << std::endl;
                        // 		std::cout << "Photopeak Energy Resolution FWHM for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakEnergyResolution() << std::endl;
                        //Compute the energy Tcut
                        std::stringstream streamEnergyCut;
                        EnergyCutMin = gauss->GetParameter(1) - photopeakSigmasMin*std::abs(gauss->GetParameter(2));
                        EnergyCutMax = gauss->GetParameter(1) + photopeakSigmasMax*std::abs(gauss->GetParameter(2));
                        streamEnergyCut << SumChannels << " > " << EnergyCutMin << " && " << SumChannels << " < " << EnergyCutMax;
                        PhotopeakEnergyCut  = streamEnergyCut.str().c_str();
                        sname.str("");

                        // then prepare the highlighted spectrum and store it in the crystal
                        sname << "Hg Charge Spectrum - Crystal " << CurrentCrystal->GetID();
                        var << SumChannels << " >> " << sname.str();
                        TH1F* spectrumChargeHighlighted = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                        tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                        spectrumChargeHighlighted->GetXaxis()->SetTitle("ADC Channels");
                        spectrumChargeHighlighted->GetYaxis()->SetTitle("N");
                        CurrentCrystal->SetHighlightedSpectrum(spectrumChargeHighlighted);
                        var.str("");
                        sname.str("");
                      }
                      //-----------------------------------------------------------------------
                      CurrentCrystal->SetSpectrum(spectrumCharge);

                      //w histogram with cut on crystal events, xyz and trigger channel and cut on photopeak
                      //if it's a background run, it won't cut on photopeak (since there will be none and the string will be empty)
                      sname << "W histogram - Crystal " << CurrentCrystal->GetID();
                      //var << "(ch" << channel << "/(" << SumChannels << ")) >> " << sname.str();
                      var << "FloodZ >> " << sname.str();
                      TH1F* spectrumHistoW = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                      spectrumHistoW->GetXaxis()->SetTitle("W");
                      spectrumHistoW->GetYaxis()->SetTitle("N");
                      // 		int bin1 = spectrumHistoW->FindFirstBinAbove(spectrumHistoW->GetMaximum()/2.0);
                      // 		int bin2 = spectrumHistoW->FindLastBinAbove(spectrumHistoW->GetMaximum()/2.0);
                      bin3 = spectrumHistoW->FindFirstBinAbove(wThreshold*spectrumHistoW->GetMaximum());
                      bin4 = spectrumHistoW->FindLastBinAbove(wThreshold*spectrumHistoW->GetMaximum());
                      // for the linear fit, let's try to stay in the really linear part of the crystal (in terms of ADCch vs w)
                      // so the fit rang will not be between these two points but closer to the center,
                      // for the moment let's say the central 50% of the w width
                      double wmin = spectrumHistoW->GetBinCenter(bin3);
                      double wmax = spectrumHistoW->GetBinCenter(bin4);
                      //DEBUG
                      // 		    std::cout << bin3 << " " << bin4 << " " << wmin << " " << wmax << std::endl;
                      meanW20 = (wmax + wmin) / 2.0;
                      // 		    double WhalfWidth  = (wmax-wmin)/2.0;
                      wbin3 = wmin + (wmax-wmin)*energyCorrectionMin;
                      wbin4 = wmin + (wmax-wmin)*energyCorrectionMax;
                      // 		    wbin3 = meanW20 - energyCorrectionMin*WhalfWidth;
                      // 		    wbin4 = meanW20 + 0.25*WhalfWidth;
                      std::stringstream ssCut20w;
                      ssCut20w << "(FloodZ) > " << spectrumHistoW->GetBinCenter(bin3) << " && " << "(FloodZ) < "<<  spectrumHistoW->GetBinCenter(bin4);
                      TCut w20percCut = ssCut20w.str().c_str();  //cut for w to get only the "relevant" part - TODO find a reasonable way to define this


                      // 		//DEBUG
                      // 		std::cout << bin3 << " " << bin4 << " " << spectrumHistoW->GetBinCenter(bin3) << " " << spectrumHistoW->GetBinCenter(bin4) << " " << meanW20 << std::endl;

                      CurrentCrystal->SetW20percCut(w20percCut);
                      // 		double width20perc =spectrumHistoW->GetBinCenter(bin4) - spectrumHistoW->GetBinCenter(bin3);
                      // 		double fwhm = spectrumHistoW->GetBinCenter(bin2) - spectrumHistoW->GetBinCenter(bin1);
                      // 		double rms = spectrumHistoW->GetRMS();
                      CurrentCrystal->SetHistoW(spectrumHistoW);
                      // 		CurrentCrystal->SetHistoWfwhm(fwhm);
                      // 		CurrentCrystal->SetHistoWrms(rms);
                      // 		CurrentCrystal->SetHistoWwidth20perc(width20perc);
                      var.str("");
                      sname.str("");

                      if(usingTaggingBench || calcDoiResWithCalibration) //done on the standard W because there's no ADCvsW correction in tagging bench
                      {
                        sname << "gaussW histogram - Crystal " << CurrentCrystal->GetID();
                        TF1 *gaussW = new TF1(sname.str().c_str(),  "gaus",histo3Dmin,histo3Dmax);
                        int binmax = spectrumHistoW->GetMaximumBin();
                        double maximum = spectrumHistoW->GetXaxis()->GetBinCenter(binmax);
                        int nentries= spectrumHistoW->GetEntries();
                        gaussW->SetParameter(0,maximum);
                        gaussW->SetParameter(1,spectrumHistoW->GetMean());
                        gaussW->SetParameter(2,spectrumHistoW->GetRMS());
                        Int_t fitStatus = spectrumHistoW->Fit(sname.str().c_str(),"QR");
                        if(!fitStatus)
                        doiFile << CurrentCrystal->GetI() + doiColumnOffset << "\t" << CurrentCrystal->GetJ() << "\t" << gaussW->GetParameter(1) << "\t" << taggingPosition <<"\t" << gaussW->GetParameter(2)/TMath::Sqrt(nentries) <<"\t"<<TMath::Sqrt(nentries)<< std::endl;
                        CurrentCrystal->SetHistoWfit(gaussW);
                        sname.str("");

                        if(calcDoiResWithDelta)
                        {
                          for(unsigned int kAltDoi = 0; kAltDoi < inputDoi.size(); kAltDoi++)
                          {
                            int ik = inputDoi[kAltDoi].i;
                            int jk = inputDoi[kAltDoi].j;
                            if( ((CurrentCrystal->GetI() + doiColumnOffset ) == ik) && CurrentCrystal->GetJ() == jk)
                            {
                              sname << "Alternate DOI res - Crystal " << CurrentCrystal->GetID();
                              var << "(FloodZ * " << inputDoi[kAltDoi].m << "+" << inputDoi[kAltDoi].q << ") - (ZPosition  ) >> " << sname.str();
                              TH1F* spectrumAltDoiRes = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,-10,10); //FIXME boundaries hardcoded
                              tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                              spectrumAltDoiRes->GetXaxis()->SetTitle("Doi_calculated - DoiTagPosition [mm]");
                              spectrumAltDoiRes->GetYaxis()->SetTitle("N");
                              var.str("");
                              sname.str("");

                              sname << "gaussFit Alternate - Crystal " << CurrentCrystal->GetID();
                              TF1 *gaussFitAlt = new TF1(sname.str().c_str(),  "gaus",-10,10); //FIXME boundaries hardcoded
                              Int_t fitStatusAlt = spectrumAltDoiRes->Fit(sname.str().c_str(),"QR");
                              if(!fitStatusAlt)
                              {
                                AltDoiFile << ik << " " << jk << " " << taggingPosition << " " << gaussFitAlt->GetParameter(2)*2.355 << std::endl;
                              }


                              CurrentCrystal->SetHistoAltDoiRes(spectrumAltDoiRes);
                              CurrentCrystal->SetHistoAltDoiFit(gaussFitAlt);
                              sname.str("");
                              // for(int iik =0 ; iik < pointsFromDoi ; iik++)
                              // {
                              //   sigmaWdoiCentral->Fill(inputDoi[kAltDoi].sw[iik] * inputDoi[kAltDoi].sqrt_nentries[iik]);
                              // }
                            }
                          }
                        }
                      }



                      if(!backgroundRun)
                      {
                        //histogram of w versus adc channels
                        //it will be useful for doi correction
                        //long long int nPoints;
                        sname << "ADC channels vs. W - Crystal " << CurrentCrystal->GetID();
                        var << SumChannels << ":FloodZ >> " << sname.str();
                        TH2F* spectrum2dADCversusW = new TH2F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax,histo1Dbins,0,histo1Dmax);
                        tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+w20percCut,"COLZ");
                        spectrum2dADCversusW->GetXaxis()->SetTitle("W");
                        spectrum2dADCversusW->GetYaxis()->SetTitle("ADC channels");
                        double parM;
                        double parQ;
                        var.str("");
                        sname.str("");


                        if(correctingForDOI)
                        {
                          //mayhem
                          //i.e. use FitSlicesY to get the gaussian fit of each slices of the TH2F. Sliced in x (bin by bin)
                          // first define the gaussian function
                          //TF1 *gaussFitSlice = new TF1("gaussFitSlice","[0]*exp(-0.5*((x-[1])/[2])**2)",EnergyCutMin,EnergyCutMax);
                          //gaussFitSlice->SetParameter(1,(EnergyCutMin+EnergyCutMax)/2.0);
                          //gaussFitSlice->SetParameter(2,0.15*(EnergyCutMin+EnergyCutMax)/2.0);
                          // 		    gaussFitSlice->SetRange();
                          spectrum2dADCversusW->FitSlicesY(0, bin3, bin4, 0, "QNRG5S");
                          sname << spectrum2dADCversusW->GetName() << "_1";
                          TH1D *spectrum2d_1 = (TH1D*)gDirectory->Get(sname.str().c_str()); // _1 is the TH1D automatically created by ROOT when FitSlicesX is called, holding the TH1F of the mean values
                          sname.str("");
                          sname << "linearCrystal - Crystal " << CurrentCrystal->GetID();
                          TF1 *linearCrystal = new TF1(sname.str().c_str(),  "[0]*x + [1]",wbin3,wbin4);
                          spectrum2d_1->Fit(sname.str().c_str(),"QR");

                          parM = linearCrystal->GetParameter(0); // m parameter for the linear fit to correct energy res for DOI
                          parQ = linearCrystal->GetParameter(1); // q parameter for the linear fit to correct energy res for DOI
                          CurrentCrystal->SetSlicesMean(spectrum2d_1);
                          CurrentCrystal->SetSlicesMeanFit(linearCrystal);
                          sname.str("");
                        }

                        // 		spectrum2dADCversusW->SetName(sname.str().c_str());
                        CurrentCrystal->SetADCversusW(spectrum2dADCversusW);

                        if(correctingForDOI)
                        {
                          //spectrum corrected for DOI
                          sname << "Charge Spectrum Corrected - Crystal " << CurrentCrystal->GetID();
                          std::stringstream baseVar;

                          // 		    baseVar << "(("  <<  SumChannels<< " ) - ( ( FloodZ - " <<  meanW20 << " ) * ( " << parM << ") ))";
                          //calculate ADC(w0), it's always the same. w0=meanw20
                          // 		    double adcw0 = parM * meanW20 + parQ;

                          baseVar << "(("  <<  SumChannels<< " ) * ( (" << parM * meanW20 + parQ << ") / (" << parM << "* FloodZ + "<< parQ << " )))";
                          // 		    std::cout << std::endl;
                          // 		    std::cout << bin3 << " " << bin4 << " "  << baseVar.str() << std::endl;
                          var << baseVar.str() << " >> " << sname.str();
                          TH1F* spectrumChargeCorrected = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                          tree->Draw(var.str().c_str(),CrystalCut,"COLZ");
                          spectrumChargeCorrected->GetXaxis()->SetTitle("ADC channels");
                          spectrumChargeCorrected->GetYaxis()->SetTitle("Counts");
                          CurrentCrystal->SetCorrectedSpectrum(spectrumChargeCorrected);
                          var.str("");
                          sname.str("");


                          spectrumChargeCorrected->GetXaxis()->SetRangeUser(peakSearchRangeMin,peakSearchRangeMax);
                          TSpectrum *s_corr;
                          s_corr = new TSpectrum(20);
                          // 		Input[i].SumSpectraCanvas->cd(j+1);


                          Int_t CrystalPeaksN_corr = s_corr->Search(spectrumChargeCorrected,2,"goff",0.5);
                          Double_t *CrystalPeaks_corr =  s_corr->GetPositionX();
                          Double_t *CrystalPeaksY_corr = s_corr->GetPositionY();
                          // 		  delete s_corr;
                          float maxPeak_corr = 0.0;
                          int peakID_corr = 0;

                          // std::cout << "-----------------------------------------------" << std::endl;
                          for (int peakCounter = 0 ; peakCounter < CrystalPeaksN_corr ; peakCounter++ )
                          {
                            //DEBUG
                            // std::cout << peakCounter << " " << CrystalPeaks_corr[peakID_corr] << " " << CrystalPeaksY_corr[peakID_corr] << std::endl;
                            // if( (CrystalPeaks_corr[peakCounter] >= peakSearchRangeMin) && (CrystalPeaks_corr[peakCounter] <= peakSearchRangeMax) ) //look for the 511 peak only in the selected range (which could be all histogram if nothing is set)
                            // {
                            if(CrystalPeaksY_corr[peakCounter] > maxPeak_corr) // then look for tallest peak in the range
                            {
                              maxPeak_corr = CrystalPeaksY_corr[peakCounter];
                              peakID_corr = peakCounter;
                            }
                            // }
                          }
                          //DEBUG
                          // std::cout << "FINAL - " << CurrentCrystal->GetID() << " " << CrystalPeaks_corr[peakID_corr] << " " << CrystalPeaksY_corr[peakID_corr] << std::endl;
                          spectrumChargeCorrected->GetXaxis()->SetRangeUser(1,histo1Dmax);
                          //std::cout << CrystalPeaks[0] << std::endl;
                          //std::cout << CrystalPeaksY[0] << std::endl;
                          //fit the spectra - TODO use the gaussian plus fermi?
                          //float energyResolution;
                          if (energyResolution == 0)
                          {
                            if (correctingSaturation)
                            energyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
                            else
                            energyResolution = ENERGY_RESOLUTION;
                          }
                          float par0_corr = CrystalPeaksY_corr[peakID_corr];
                          float par1_corr = CrystalPeaks_corr[peakID_corr];
                          float par2_corr = (CrystalPeaks_corr[peakID_corr]*energyResolution)/2.35;
                          float fitmin_corr = par1_corr-1.4*par2_corr;
                          float fitmax_corr = par1_corr+1.5*par2_corr;

                          //DEBUG
                          // std::cout << fitmin_corr << " " << fitmax_corr << std::endl;

                          sname << "gauss_corr - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                          TF1 *gauss_corr = new TF1(sname.str().c_str(),  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitmin_corr,fitmax_corr);
                          gauss_corr->SetParameter(0,par0_corr);
                          gauss_corr->SetParameter(1,par1_corr);
                          gauss_corr->SetParameter(2,par2_corr);
                          spectrumChargeCorrected->Fit(sname.str().c_str(),"Q","",fitmin_corr,fitmax_corr);
                          //store the mean and sigma in the crystal
                          if(gauss_corr->GetParameter(1) > 0) // otherwise the fit was very wrong..)
                          CurrentCrystal->SetPhotopeakCorrected(gauss_corr->GetParameter(1),std::abs(gauss_corr->GetParameter(2)));
                          // 		  else
                          // 		    CurrentCrystal->SetPhotopeakCorrected(1,1);
                          CurrentCrystal->SetFitCorrected(gauss_corr);
                          // 		std::cout << "Photopeak Mean for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakPosition() << std::endl;
                          // 		std::cout << "Photopeak Sigma for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakSigma() << std::endl;
                          // 		std::cout << "Photopeak Energy Resolution FWHM for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakEnergyResolution() << std::endl;
                          //Compute the energy Tcut
                          std::stringstream streamEnergyCutCorrected;
                          streamEnergyCutCorrected << baseVar.str() << " > " << gauss_corr->GetParameter(1) - photopeakSigmasMin*std::abs(gauss_corr->GetParameter(2)) << " && " << baseVar.str() << " < " << gauss_corr->GetParameter(1) + photopeakSigmasMax*std::abs(gauss_corr->GetParameter(2));
                          PhotopeakEnergyCutCorrected = streamEnergyCutCorrected.str().c_str();

                          // then PhotopeakEnergyCut becomes PhotopeakEnergyCutCorrected

                          PhotopeakEnergyCut = PhotopeakEnergyCutCorrected;
                          // 		CurrentCrystal->SetSpectrum(*spectrum);
                          sname.str("");

                          //then prepare the highlighted spectrum and store it in the crystal
                          sname << "Hg Charge Spectrum Correctd - Crystal " << CurrentCrystal->GetID();
                          // 			var << "("  <<  SumChannels<< " ) - ( ( FloodZ - " <<  meanW20 << " ) * ( " << parM << ") ) >> " << sname.str();
                          var << baseVar.str() << " >> " << sname.str();
                          TH1F* spectrumChargeCorrectedHighlighted = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                          tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                          spectrumChargeCorrectedHighlighted->GetXaxis()->SetTitle("ADC Channels");
                          spectrumChargeCorrectedHighlighted->GetYaxis()->SetTitle("N");
                          CurrentCrystal->SetHighlightedSpectrumCorrected(spectrumChargeCorrectedHighlighted);
                          var.str("");
                          sname.str("");

                        }
                      }

                      // save the TCuts
                      TCut CrystalCutWithoutCutG = BasicCut + CutTrigger;               // useful for the calcDoiResWithCalibration external part
                      CrystalCut.SetName("CrystalCut");
                      PhotopeakEnergyCut.SetName("PhotopeakEnergyCut");
                      CrystalCutWithoutCutG.SetName("CrystalCutWithoutCutG");
                      CurrentCrystal->SetCrystalCut(CrystalCut);                        // this is BasicCut + CutTrigger + cutg[0] + cutg[1] - so "geometrical" without energy constrains
                      CurrentCrystal->SetPhotopeakEnergyCut(PhotopeakEnergyCut);        // energy cut, events in the photopeak
                      CurrentCrystal->SetCrystalCutWithoutCutG(CrystalCutWithoutCutG);  // this is BasicCut + CutTrigger


                      // Histogram 2d of time evolution
                      sname << "ADC channels vs. Time - Crystal " << CurrentCrystal->GetID();
                      var << SumChannels << ":ExtendedTimeTag >> " << sname.str();
                      TH2F* spectrum2dVersusTime = new TH2F(sname.str().c_str(),sname.str().c_str(),250,0,tree->GetMaximum("ExtendedTimeTag"),histo1Dbins,0,histo1Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut,"COLZ");
                      spectrum2dVersusTime->GetXaxis()->SetTitle("ExtendedTimeTag");
                      spectrum2dVersusTime->GetYaxis()->SetTitle("ADC channels");
                      CurrentCrystal->SetVersusTime(spectrum2dVersusTime);
                      var.str("");
                      sname.str("");

                      // Time evolution for W
                      sname << "W vs. Time - Crystal " << CurrentCrystal->GetID();
                      var << "FloodZ:ExtendedTimeTag >> " << sname.str();
                      TH2F* spectrum2dWversusTime = new TH2F(sname.str().c_str(),sname.str().c_str(),250,0,tree->GetMaximum("ExtendedTimeTag"),wHistogramsBins,histo3Dmin,histo3Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut,"COLZ");
                      spectrum2dWversusTime->GetXaxis()->SetTitle("ExtendedTimeTag");
                      spectrum2dWversusTime->GetYaxis()->SetTitle("W");
                      CurrentCrystal->SetWversusTime(spectrum2dWversusTime);
                      var.str("");
                      sname.str("");

                      //histogram of w versus adc channels - this time without the photopeak cut (so it looks nicer in the paper...)
                      sname << "Complete ADC channels vs. W - Crystal " << CurrentCrystal->GetID();
                      var << SumChannels << ":FloodZ >> " << sname.str() ;
                      TH2F* spectrum2dADCversusWComplete = new TH2F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax,histo1Dbins,0,histo1Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut,"COLZ");
                      spectrum2dADCversusWComplete->GetXaxis()->SetTitle("W");
                      spectrum2dADCversusWComplete->GetYaxis()->SetTitle("ADC channels");
                      spectrum2dADCversusWComplete->SetName(sname.str().c_str());
                      CurrentCrystal->SetADCversusWComplete(spectrum2dADCversusWComplete);
                      var.str("");
                      sname.str("");

                      // w histogram with energy cut on the corrected spectrum and
                      // fit the w histogram with the thetaFunction -> m and q of correlation line
                      // fit on the rise of w function with gaussian -> sigma w
                      // no cut on the photopeak will be applied if it's a background run (since there will be no photopeak)
                      sname << "W histogram Corrected - Crystal " << CurrentCrystal->GetID();
                      var << "(FloodZ) >> " << sname.str();
                      TH1F* spectrumHistoWCorrected = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax);
                      tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                      spectrumHistoWCorrected->GetXaxis()->SetTitle("W");
                      spectrumHistoWCorrected->GetYaxis()->SetTitle("N");
                      var.str("");
                      sname.str("");


                      CurrentCrystal->SetHistoWCorrected(spectrumHistoWCorrected);





                      TGraph *calibrationGraphChosen = NULL;
                      float centralW = -1;
                      float beginW = -1;
                      float endW = -1;

                      if(!calcDoiResWithCalibration) // only if it's not calcDoiResWithCalibration, because in that case we load the calibration TGraph from another file
                      {
                        var.str("");
                        sname.str("");

                        //cumulative function
                        //first, take the w corrected and normalize to integral -> get the PDF
                        TH1F *pdfW = (TH1F*) spectrumHistoWCorrected->Clone();
                        sname << "PDF W histogram Corrected - Crystal " << CurrentCrystal->GetID();
                        pdfW->SetName(sname.str().c_str());
                        pdfW->SetTitle(sname.str().c_str());
                        double integral = pdfW->Integral();
                        pdfW->Scale(wHistogramsBins*(1.0/integral));
                        CurrentCrystal->SetPdfW(pdfW);
                        sname.str("");
                        // then do the comulative histo
                        //and the calibration Graph

                        sname << "Cumulative W histogram Corrected - Crystal " << CurrentCrystal->GetID();
                        // 		    var << "(ch" << channel << "/(" << SumChannels << ")) >> " << sname.str();
                        TH1F* cumulativeW = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax);
                        cumulativeW->GetXaxis()->SetTitle("W");
                        double sumPdf = 0;
                        // 		    std::ofstream wfile;
                        // 		    wfile.open("cumulative2_2.dat", std::ofstream::out);
                        std::vector<double> calibrationW;
                        std::vector<double> calibrationZ;
                        // 		    double lambda511 = 12.195; //everything in mm

                        for(int iPdfHisto = 0 ; iPdfHisto < wHistogramsBins; iPdfHisto++)
                        {
                          sumPdf += (pdfW->GetBinContent(iPdfHisto+1))/wHistogramsBins;
                          calibrationW.push_back(cumulativeW->GetBinCenter(iPdfHisto+1));
                          if(backgroundRun | lateralRun)
                          {
                            calibrationZ.push_back(-(sumPdf*crystalz) + crystalz); // if it's a background run, interaction probability is constant everywhere, so z is just a rescale of the cumulative
                          }
                          else
                          {
                            // if it's a far source run, interaction probability is exponential, and relation between z and w is give by
                            // z = L + l * ln( 1 - (1 - exp(-L/l))*integral_0^w(PDF(w)dw) )
                            // where
                            // L = crystal length
                            // l = interaction length of 511 gammas in lyso (from literature, 12.195 mm)
                            // integral_0^w(PDF(w)dw) = the cumulative of PDF(w) from 0 to w
                            calibrationZ.push_back( crystalz +  lambda511 * TMath::Log( 1.0 - (1.0 - TMath::Exp(-(crystalz/lambda511)) )* sumPdf )  );
                          }
                          // 		      wfile << cumulativeW->GetBinCenter(iPdfHisto+1) << " " << -(sumPdf*15) + 15.0 << std::endl;
                          cumulativeW->Fill(cumulativeW->GetBinCenter(iPdfHisto+1),sumPdf);
                        }

                        TGraph *calibrationGraph = new TGraph(calibrationW.size(),&calibrationW[0],&calibrationZ[0]);
                        sname.str("");
                        sname << "Calibration Plot - Crystal " << CurrentCrystal->GetID();
                        calibrationGraph->SetTitle(sname.str().c_str());
                        calibrationGraph->SetName(sname.str().c_str());
                        calibrationGraph->GetXaxis()->SetTitle("W");
                        calibrationGraph->GetYaxis()->SetTitle("Z [mm]");
                        sname.str("");
                        CurrentCrystal->SetCalibrationGraph(calibrationGraph);
                        calibrationGraphChosen = calibrationGraph;

                        //w(z) for compton analysis
                        TGraph *wzgraph = new TGraph(calibrationW.size(),&calibrationZ[0],&calibrationW[0]);
                        sname.str("");
                        sname << "w(z) Plot - Crystal " << CurrentCrystal->GetID();
                        wzgraph->SetTitle(sname.str().c_str());
                        wzgraph->SetName(sname.str().c_str());
                        wzgraph->GetXaxis()->SetTitle("Z [mm]");
                        wzgraph->GetYaxis()->SetTitle("W");
                        centralW = wzgraph->Eval(crystalz/2.0);
                        beginW = wzgraph->Eval(crystalz - 0.01);
                        endW = wzgraph->Eval(0.01);
                        // std::cout << beginW << "\t" << endW << std::endl;

                        sname.str("");

                        CurrentCrystal->SetWZgraph(wzgraph);
                        CurrentCrystal->SetCumulativeW(cumulativeW);
                      }
                      else // look for the calibGraph in an external calibration file
                      {
                        //save the crystal identification
                        std::stringstream MppcDirStream;
                        MppcDirStream << "MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " - 0.0-" << iMppc << "." << jMppc;
                        std::stringstream CrystalDirStream;
                        CrystalDirStream << "Crystal " <<  CurrentCrystal->GetID();
                        std::stringstream directory;
                        directory << "Module 0.0/" << MppcDirStream.str() << "/" << CrystalDirStream.str();
                        bool dir_exists = calibrationFile->cd(directory.str().c_str());
                        TCanvas* C_graph = NULL;
                        TGraph *calibGraph = NULL;
                        if(dir_exists)
                        {
                          //take the graph
                          std::stringstream stream;
                          stream << "Calibration Plot - Crystal " << CurrentCrystal->GetID();
                          C_graph = (TCanvas*) gDirectory->Get(stream.str().c_str());
                          if(C_graph)
                            calibGraph = (TGraph*) C_graph->GetPrimitive(stream.str().c_str());

                          var.str("");
                          sname.str("");
                        }
                        CurrentCrystal->SetCalibrationGraph(calibGraph);
                        calibrationGraphChosen = calibGraph;
                      }


                      if(usingTaggingBench || taggingForTiming)
                      {


                        if( (digitizerType == 1) || (digitizerType == 2)) //only for digitizers with timing capabilities
                        {
                          // basic CTR with respect to an external tagging crystal, cutting on photopeak of this crystal
                          // var.str("");
                          // sname.str("");
                          sname << "Basic CTR - Crystal " << CurrentCrystal->GetID();
                          var << "t" << channel << "- TaggingTimeStamp  >> " << sname.str();
                          // std::cout << sname.str() << std::endl
                          // << var.str()   << std::endl
                          // << CrystalCut+PhotopeakEnergyCut << std::endl;

                          TH1F* aSpectrum = new TH1F(sname.str().c_str(),sname.str().c_str(),CTRbins,CTRmin,CTRmax);
                          tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                          aSpectrum->GetXaxis()->SetTitle("Time [S]");
                          aSpectrum->GetYaxis()->SetTitle("N");
                          CurrentCrystal->SetDeltaTimeWRTTagging(aSpectrum);
                          var.str("");
                          sname.str("");

                          // scatter plot of the above delta and w (which correlates to doi)
                          sname << "Delta T vs. W - Crystal " << CurrentCrystal->GetID();
                          var << "t" << channel << " - TaggingTimeStamp :FloodZ >> " << sname.str() ;
                          TH2F* spectrumCrystalDeltaTvsW = new TH2F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax,CTRbins,CTRmin,CTRmax);
                          tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut,"COLZ");
                          spectrumCrystalDeltaTvsW->GetXaxis()->SetTitle("W");
                          spectrumCrystalDeltaTvsW->GetYaxis()->SetTitle("T crystal - T tagging [S]");
                          CurrentCrystal->SetDeltaTvsW(spectrumCrystalDeltaTvsW);
                          var.str("");
                          sname.str("");

                          // spectrumCrystalDeltaTvsW->FitSlicesY(0, 0, -1, 0, "QNRG5S");
                          // sname << spectrumCrystalDeltaTvsW->GetName() << "_1";
                          // TH1D *spectrumCrystalDeltaTvsW_1 = (TH1D*)gDirectory->Get(sname.str().c_str()); // _1 is the TH1D automatically created by ROOT when FitSlicesY is called, holding the TH1F of the mean values
                          // sname.str("");
                          // sname << "linearFit CTR - Crystal " << CurrentCrystal->GetID();
                          // TF1 *linearFitCTR = new TF1(sname.str().c_str(),  "[0]*x + [1]",histo3Dmin,histo3Dmax);
                          // spectrumCrystalDeltaTvsW_1->Fit(sname.str().c_str(),"QR");
                          // CurrentCrystal->SetDeltaTvsWSlice(spectrumCrystalDeltaTvsW_1);
                          // CurrentCrystal->SetDeltaTvsWFit(linearFitCTR);
                          // var.str("");
                          // sname.str("");

                          if(timingCorrection)
                          {

                            //get a TGraph
                            std::vector<float> delta_X;
                            std::vector<float> W_Y;
                            std::vector<float> rmsBasic_Y;

                            //get a TGraph from FitSlicesY
                            // for(int iBin = 0; iBin < spectrumCrystalDeltaTvsW_1->GetNbinsX() ; iBin++)
                            // {
                            //   delta_X.push_back(spectrumCrystalDeltaTvsW_1->GetBinCenter(iBin+1));
                            //   W_Y.push_back(spectrumCrystalDeltaTvsW_1->GetBinContent(iBin+1));
                            // }


                            //Tgraph from mean and rms of "slices"
                            for(int iBin = -1; iBin < WrangeBinsForTiming+1; iBin++) // -1 and +1 are put to include the w limits
                            {

                              // first do an "amplitude" vs. ctr histo, for this w range, and use it to correct
                              // the CTRs for of this slice before generating the spectrum from which to extract the
                              // delta_X and rmsBasic_Y

                              // define the Limits and mid point for this w slice
                              Float_t wLowerLimit = beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming);
                              Float_t wUpperLimit = beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming);
                              Float_t midW        = beginW + (((iBin+0.5)*(endW - beginW))/WrangeBinsForTiming);
                              //write the TCut for this data
                              std::stringstream sCut;
                              sCut << "FloodZ > "
                              << wLowerLimit
                              << " && FloodZ < "
                              << wUpperLimit
                              << " && (t"
                              << channel
                              << " - TaggingTimeStamp) > "
                              << CTRmin
                              << " && (t"
                              << channel
                              << " - TaggingTimeStamp) < "
                              << CTRmax;
                              TCut wCut = sCut.str().c_str();

                              sname.str("");



                              // // ------- BEGIN OF MODS FOR AMPL CORRECTION
                              // sname << "T_vs_Q_" << iBin << " - Crystal " << CurrentCrystal->GetID();
                              // TH2F *TvsQ = new TH2F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,0,histo1Dmax,CTRbins,CTRmin,CTRmax);
                              // var << "t" << channel << " - TaggingTimeStamp : ch" << channel <<  " >>  " << sname.str() ;
                              // tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+wCut);
                              // // ------- END OF MODS FOR AMPL CORRECTION





                              var.str("");
                              sname.str("");
                              sname << "DeltaT_" << iBin << " - Crystal " << CurrentCrystal->GetID();
                              TH1F *tempHisto = new TH1F(sname.str().c_str(),sname.str().c_str(),CTRbins,CTRmin,CTRmax);
                              var << "t" << channel << " - TaggingTimeStamp  >> " << sname.str() ;


                              tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+wCut);
                              delta_X.push_back(midW);
                              W_Y.push_back(tempHisto->GetMean());
                              rmsBasic_Y.push_back(tempHisto->GetRMS());
                              // std::cout << iBin << "\t"
                              //           << beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming) << "\t"
                              //           << beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming) << "\t"
                              //           << beginW +  (((iBin+0.5)*(endW - beginW))/WrangeBinsForTiming) << std::endl;
                              var.str("");
                              sname.str("");


                              // // ------- BEGIN OF MODS FOR AMPL CORRECTION
                              // CurrentCrystal->AddTvsQHistos(TvsQ);
                              // CurrentCrystal->AddDeltaTHistos(tempHisto);
                              // // ------- END OF MODS FOR AMPL CORRECTION




                            }


                            TGraph *graphDeltaW = new TGraph(delta_X.size(),&delta_X[0],&W_Y[0]);
                            sname.str("");
                            sname << "DeltaW Graph - Crystal " << CurrentCrystal->GetID();
                            graphDeltaW->SetTitle(sname.str().c_str());
                            graphDeltaW->SetName(sname.str().c_str());
                            graphDeltaW->GetXaxis()->SetTitle("W");
                            graphDeltaW->GetYaxis()->SetTitle("T crystal - T tagging [S]");
                            sname.str("");
                            CurrentCrystal->SetGraphDeltaW(graphDeltaW);

                            TGraph *graphDeltaRMS = new TGraph(delta_X.size(),&delta_X[0],&rmsBasic_Y[0]);
                            sname.str("");
                            sname << "RMS DeltaW Graph - Crystal " << CurrentCrystal->GetID();
                            graphDeltaRMS->SetTitle(sname.str().c_str());
                            graphDeltaRMS->SetName(sname.str().c_str());
                            graphDeltaRMS->GetXaxis()->SetTitle("W");
                            graphDeltaRMS->GetYaxis()->SetTitle("T crystal - T tagging [S]");
                            sname.str("");
                            CurrentCrystal->SetGraphDeltaRMS(graphDeltaRMS);

                            //calc w for z=(crystalz/2)

                            // float centralCTR = linearFitCTR->GetParameter(0)*centralW + linearFitCTR->GetParameter(1);
                            // std::cout << "centralCTR = " << centralCTR << std::endl;
                            // sname << "CTR central correction - Crystal " << CurrentCrystal->GetID();
                            // var << "((t" << channel << " - ((" << linearFitCTR->GetParameter(0) << "* FloodZ) + " << linearFitCTR->GetParameter(1) << ") - ("<< centralCTR <<")) - TaggingTimeStamp) >> " << sname.str() ;
                            // std::cout << "var = " << var.str() << std::endl;
                            // TH1F* spectrumCTRcentralCorrection = new TH1F(sname.str().c_str(),sname.str().c_str(),CTRbins,CTRmin,CTRmax);
                            // tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                            // spectrumCTRcentralCorrection->GetXaxis()->SetTitle("Time [S]");
                            // spectrumCTRcentralCorrection->GetYaxis()->SetTitle("N");
                            // CurrentCrystal->SetCTRcentralCorrection(spectrumCTRcentralCorrection);
                            // var.str("");
                            // sname.str("");

                            // plots for the neighbour channels (channels, NOT crystals!) of T cry - T neighbour
                            // one for each neighbour
                            for(unsigned int iNeig = 0; iNeig < neighbours.size(); iNeig++)
                            {
                              // histogram of tCry - tNeighbour
                              sname <<  "T_Channel_" << neighbours[iNeig] << " - T_Crystal_" << CurrentCrystal->GetID();
                              var << "t" << neighbours[iNeig] << " - t" << channel << " >> " << sname.str() ;
                              TH1F* spectrumDeltaTcryTneig = new TH1F(sname.str().c_str(),sname.str().c_str(),DeltaTimeBins,DeltaTimeMin,DeltaTimeMax);
                              tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut);
                              spectrumDeltaTcryTneig->GetXaxis()->SetTitle("Time [S]");
                              spectrumDeltaTcryTneig->GetYaxis()->SetTitle("N");
                              var.str("");
                              sname.str("");

                              sname << "T_Channel_" << neighbours[iNeig] << " - T_Crystal_" << CurrentCrystal->GetID() << " vs. W";
                              var << "t" << neighbours[iNeig] << " - t" << channel << ":FloodZ >> " << sname.str() ;
                              TH2F* spectrumCrystalDeltaT2vsW = new TH2F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax,DeltaTimeBins,DeltaTimeMin,DeltaTimeMax);
                              tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut,"COLZ");
                              spectrumCrystalDeltaT2vsW->GetXaxis()->SetTitle("W");
                              sname.str("");
                              sname << "T_channel_"<< neighbours[iNeig] << " - T_crystal " << CurrentCrystal->GetID() << ", [S]";
                              spectrumCrystalDeltaT2vsW->GetYaxis()->SetTitle(sname.str().c_str());
                              var.str("");
                              sname.str("");

                              std::vector<float> delay_X;
                              std::vector<float> Wcoord_Y;
                              std::vector<float> rms_Y;


                              //get TGraphs from building N th1f, in the range defined by beginW and endW previously found
                              // the range beginw-endW is spilt in WrangeBinsForTiming parts,
                              for(int iBin = -1; iBin < WrangeBinsForTiming+1; iBin++) // -1 and +1 are put to include the w limits
                              {
                                TH1F *tempHisto = new TH1F("tempHisto","tempHisto",DeltaTimeBins,DeltaTimeMin,DeltaTimeMax);
                                var << "t" << neighbours[iNeig] << " - t" << channel << " >> tempHisto";
                                std::stringstream sCut;
                                sCut << "FloodZ > "
                                << beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming)
                                << "&& FloodZ < "
                                << beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming)
                                << "&& (t" << neighbours[iNeig] << " - t" << channel << ") > " << DeltaTimeMin
                                << "&& (t" << neighbours[iNeig] << " - t" << channel << ") < " << DeltaTimeMax;
                                TCut wCut = sCut.str().c_str();
                                tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut+wCut);

                                delay_X.push_back( beginW + (((iBin+0.5)*(endW - beginW))/WrangeBinsForTiming) );
                                Wcoord_Y.push_back(tempHisto->GetMean());
                                rms_Y.push_back(tempHisto->GetRMS());
                                // std::cout << iBin << "\t"
                                //           << beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming) << "\t"
                                //           << beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming) << "\t"
                                //           << beginW +  (((iBin+0.5)*(endW - beginW))/WrangeBinsForTiming) << std::endl;
                                var.str("");
                              }

                              // Fit Slice method to derive the correction graphs above.
                              // in this case we use FitSlicesY to get directly the fit of a slice (possibly the sum of more slices in w) using a function that is gauss+exp
                              // spectrumCrystalDeltaT2vsW->FitSlicesY(0, 0, -1, 0, "QNR");
                              // sname << spectrumCrystalDeltaT2vsW->GetName() << "_1";
                              // TH1D *spectrumCrystalDeltaT2vsW_1 = (TH1D*)gDirectory->Get(sname.str().c_str()); // _1 is the TH1D automatically created by ROOT when FitSlicesY is called, holding the TH1F of the mean values
                              // sname.str("");
                              // sname << spectrumCrystalDeltaT2vsW->GetName() << "_2";
                              // TH1D *spectrumCrystalDeltaT2vsW_2 = (TH1D*)gDirectory->Get(sname.str().c_str());
                              // sname.str("");
                              //
                              // for(int iBin = 0; iBin < spectrumCrystalDeltaT2vsW_1->GetNbinsX() ; iBin++)
                              // {
                              //   if(spectrumCrystalDeltaT2vsW_2->GetBinContent(iBin) )
                              //   delay_X.push_back( spectrumCrystalDeltaT2vsW_1->GetBinCenter(iBin) );
                              //   Wcoord_Y.push_back( spectrumCrystalDeltaT2vsW_1->GetBinContent(iBin) );
                              //   var.str("");
                              // }

                              TGraph *graphDelayW = new TGraph(delay_X.size(),&delay_X[0],&Wcoord_Y[0]);
                              sname.str("");
                              sname << "Graph Delay ch_" << neighbours[iNeig];
                              // sname << "DeltaW Graph - Crystal " << CurrentCrystal->GetID();
                              graphDelayW->SetTitle(sname.str().c_str());
                              graphDelayW->SetName(sname.str().c_str());
                              graphDelayW->GetXaxis()->SetTitle("W");
                              sname.str("");
                              sname << "T_channel_"<< neighbours[iNeig] << " - T_crystal " << CurrentCrystal->GetID() << ", [S]";
                              graphDelayW->GetYaxis()->SetTitle(sname.str().c_str());
                              sname.str("");


                              TGraph *graphDelayRMS = new TGraph(delay_X.size(),&delay_X[0],&rms_Y[0]);
                              sname.str("");
                              sname << "RMS Graph Delay ch_" << neighbours[iNeig];
                              // sname << "DeltaW Graph - Crystal " << CurrentCrystal->GetID();
                              graphDelayRMS->SetTitle(sname.str().c_str());
                              graphDelayRMS->SetName(sname.str().c_str());
                              graphDelayRMS->GetXaxis()->SetTitle("W");
                              sname.str("");
                              sname << "RMS (T_channel_"<< neighbours[iNeig] << " - T_crystal " << CurrentCrystal->GetID() << "), [S]";
                              graphDelayRMS->GetYaxis()->SetTitle(sname.str().c_str());
                              sname.str("");




                              int plotPos = -1;
                              for(int iTimeMppc = 0 ; iTimeMppc < nmppcx ; iTimeMppc++)
                              {
                                for(int jTimeMppc = 0 ; jTimeMppc < nmppcy ; jTimeMppc++)
                                {
                                  if(mppc[(iModule*nmppcx)+iTimeMppc][(jModule*nmppcy)+jTimeMppc]->GetDigitizerChannel() == neighbours[iNeig])
                                  {
                                    plotPos = mppc[(iModule*nmppcx)+iTimeMppc][(jModule*nmppcy)+jTimeMppc]->GetCanvasPosition();
                                  }
                                }
                              }

                              if(plotPos != -1)
                              {
                                CurrentCrystal->AddDeltaTcryTneig(spectrumDeltaTcryTneig,plotPos);
                                CurrentCrystal->AddDeltaT2vsW(spectrumCrystalDeltaT2vsW,plotPos);;
                                CurrentCrystal->AddGraphDelayW(graphDelayW,plotPos);
                                CurrentCrystal->AddGraphDelayRMS(graphDelayRMS,plotPos);
                              }

                            }
                          }
                        }
                      }


                      //Compton calibration plots
                      //TEST on single. MPPC is now C2, crystal 28. We want to see just MPPC C3
                      //the p_C3(E,w)
                      // create a 3x3 matrix of TGraph2D*
                      if(comptonAnalysis)
                      {
                        TGraph2D***  ComptonCalibation;
                        ComptonCalibation = new TGraph2D** [nmppcx];
                        for(int iCal = 0; iCal < nmppcx ; iCal++) ComptonCalibation[iCal] = new TGraph2D*[nmppcy];

                        TGraph2D***  ConvertedComptonCalibration;
                        ConvertedComptonCalibration = new TGraph2D** [nmppcx];
                        for(int iCal = 0; iCal < nmppcx ; iCal++) ConvertedComptonCalibration[iCal] = new TGraph2D*[nmppcy];


                        for(int iComptMppc = 0; iComptMppc < nmppcx ; iComptMppc++)
                        {
                          for(int jComptMppc = 0; jComptMppc < nmppcy ; jComptMppc++)
                          {
                            sname.str("");
                            var.str("");

                            int tempChannel = mppc[iComptMppc][jComptMppc]->GetDigitizerChannel();
                            var << "(ch"<< tempChannel << "):(" << SumChannels << "):(FloodZ)" ;
                            tree->Draw(var.str().c_str(),CrystalCut);
                            TGraph2D* tempGraph2d = new TGraph2D(tree->GetSelectedRows(),tree->GetV3(),tree->GetV2(), tree->GetV1());




                            // sname << "Histo_Pi(E,w)[" << iMppc + iNeighbour <<  "][" << jMppc + jNeighbour <<  "]_" << CurrentCrystal->GetID();
                            sname << "Pi_[" << iComptMppc <<  "][" << jComptMppc <<  "]_" << CurrentCrystal->GetID();
                            TH3I* ComptonCalibrationHistogram = new TH3I(
                              sname.str().c_str(),
                              sname.str().c_str(),
                              20,
                              tempGraph2d->GetXmin(),
                              tempGraph2d->GetXmax(),
                              20,
                              tempGraph2d->GetYmin(),
                              tempGraph2d->GetYmax(),
                              100,
                              tempGraph2d->GetZmin(),
                              tempGraph2d->GetZmax());


                              var.str("");

                              // ComptonCalibrationHistogram->SetTitle(sname.str().c_str());
                              // ComptonCalibrationHistogram->SetName(sname.str().c_str());
                              ComptonCalibrationHistogram->GetXaxis()->SetTitle("W");
                              ComptonCalibrationHistogram->GetYaxis()->SetTitle("Sum Charge [ADC ch]");
                              ComptonCalibrationHistogram->GetZaxis()->SetTitle("pi [ADC ch]");
                              var.str("");
                              var << "(ch"<< tempChannel << "):(" << SumChannels << "):(FloodZ)" << " >> " << sname.str();
                              tree->Draw(var.str().c_str(),CrystalCut);
                              // std::cout << ComptonCalibrationHistogram->GetEntries() << std::endl;
                              CurrentCrystal->SetComptonHistogram(ComptonCalibrationHistogram);
                              sname.str("");
                              sname << "Graph_Pi(q,w)[" << iComptMppc <<  "][" << jComptMppc <<  "]_" << CurrentCrystal->GetID();
                              std::vector<double> ex,ey,ez;
                              for(int iError = 0; iError < tree->GetSelectedRows(); iError++)
                              {
                                ex.push_back(0.05);
                                ey.push_back(20.0);
                                ez.push_back(20.0);
                              }
                              ComptonCalibation[iComptMppc][jComptMppc] = new TGraph2D(tree->GetSelectedRows(),tree->GetV3(),tree->GetV2(), tree->GetV1());
                              ComptonCalibation[iComptMppc][jComptMppc]->SetTitle(sname.str().c_str());
                              ComptonCalibation[iComptMppc][jComptMppc]->SetName(sname.str().c_str());
                              ComptonCalibation[iComptMppc][jComptMppc]->GetXaxis()->SetTitle("W");
                              ComptonCalibation[iComptMppc][jComptMppc]->GetYaxis()->SetTitle("Sum Charge [ADC ch]");
                              ComptonCalibation[iComptMppc][jComptMppc]->GetZaxis()->SetTitle("pi [ADC ch]");
                              // CurrentCrystal->SetComptonCalibration(testGraph);
                              var.str("");
                              sname.str("");


                              //produce

                              // TGraph2D* graph2d = (TGraph2D*) canv-> GetPrimitive(sstream.str().c_str()); //get the 3d graph
                              //
                              // // int Npoints = graph2d->GetN();
                              //
                              sname << "Pi(z,p_tot)[" << iComptMppc <<  "][" << jComptMppc <<  "]_" << CurrentCrystal->GetID();
                              ConvertedComptonCalibration[iComptMppc][jComptMppc] = new TGraph2D();
                              Double_t *x = ComptonCalibation[iComptMppc][jComptMppc]->GetX();
                              Double_t *y = ComptonCalibation[iComptMppc][jComptMppc]->GetY();
                              Double_t *z = ComptonCalibation[iComptMppc][jComptMppc]->GetZ();
                              for(int i = 0 ; i < ComptonCalibation[iComptMppc][jComptMppc]->GetN() ; i++){
                                ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetPoint(i,calibrationGraphChosen->Eval(x[i]),y[i],z[i]);
                                // ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetPointError(i,2.0,sqrt(y[i]),sqrt(z[i]));
                              }
                              ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetTitle(sname.str().c_str());
                              ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetName(sname.str().c_str());
                              ConvertedComptonCalibration[iComptMppc][jComptMppc]->GetXaxis()->SetTitle("W");
                              ConvertedComptonCalibration[iComptMppc][jComptMppc]->GetYaxis()->SetTitle("Sum Charge [ADC ch]");
                              ConvertedComptonCalibration[iComptMppc][jComptMppc]->GetZaxis()->SetTitle("pi [ADC ch]");
                              sname.str("");


                            }
                          }


                          //       }
                          //     }
                          //   }
                          // }
                          CurrentCrystal->SetComptonCalibration(ComptonCalibation);
                          CurrentCrystal->SetConvertedComptonCalibration(ConvertedComptonCalibration);
                          // CurrentCrystal->SetComptonCalibrationHistogram(ComptonCalibationHistogram);
                          // CurrentCrystal->SetInterpolationGraph(interpolationGraph);
                        }


                        if(usingRealSimData) // only if this is a sim dataset
                        {
                          long long int nPoints;
                          // a 2d plot of real vs. w, using the same cuts as before


                          int RealZBins = 100;

                          sname << "Real Z vs. W - Crystal " << CurrentCrystal->GetID();
                          var << "-(RealZ-" << CurrentCrystal->GetDimensionZ()/2.0 << "):FloodZ >> " << sname.str();
                          TH2F* spectrum2dSimDOIplot = new TH2F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,histo3Dmin,histo3Dmax,RealZBins,0,CurrentCrystal->GetDimensionZ());
                          nPoints = tree->Draw(var.str().c_str(),CrystalCut+PhotopeakEnergyCut ,"COLZ");
                          spectrum2dSimDOIplot->GetXaxis()->SetTitle("W");
                          spectrum2dSimDOIplot->GetYaxis()->SetTitle("Z");
                          CurrentCrystal->SetSimDOIplot(spectrum2dSimDOIplot);
                          sname.str("");
                          var.str("");

                          // same plot but with a TGraph, to allow fitting
                          sname << "Graph Z vs. W - Crystal " << CurrentCrystal->GetID();
                          TGraph* graph = new TGraph(nPoints,tree->GetV2(),tree->GetV1()); // same but TGraph (so it can be fitted in 1D)
                          graph->SetName(sname.str().c_str());
                          graph->SetTitle(sname.str().c_str());
                          graph->GetXaxis()->SetTitle("W");
                          graph->GetYaxis()->SetTitle("Z");
                          graph->Draw("ap");
                          sname.str("");
                          // 		  TF1 *linear = new TF1("linear",  "[0]*x + [1]",0,1);
                          // 		      sname << "expfit - Crystal " << CurrentCrystal->GetID();
                          // 		      TF1 *expfit = new TF1(sname.str().c_str(),  "[0]*exp(-x/[1])",0,1);
                          // 		      // 		  linear->SetParameter(0,-100);
                          // 		      // 		  linear->SetParameter(1,50);
                          // 		      expfit->SetParameter(0,50);
                          // 		      expfit->SetParameter(1,0.1);
                          // 		      // 		  graph->SetStats(1);
                          // 		      graph->Fit(sname.str().c_str(),"Q","",0.1,0.7);
                          // 		      CurrentCrystal->SetSimFit(expfit);
                          CurrentCrystal->SetSimGraph(graph);
                          sname.str("");
                          var.str("");
                          // FitSlicesX so we can plot it together with the calibration graph
                          spectrum2dSimDOIplot->FitSlicesX(0, 1, RealZBins, 0, "QNRG5S");
                          sname << spectrum2dSimDOIplot->GetName() << "_1";
                          TH1D *spectrum2dSimDOIplot_1 = (TH1D*)gDirectory->Get(sname.str().c_str());
                          sname.str("");
                          sname << spectrum2dSimDOIplot->GetName() << "_2";
                          TH1D *spectrum2dSimDOIplot_2 = (TH1D*)gDirectory->Get(sname.str().c_str());
                          sname.str("");
                          //run on the two histograms and store the results in a TGraphErrors and in a TH1F

                          sname << "Sigma W from SIM - Crystal " << CurrentCrystal->GetID();
                          TH1F *sigmaSim = new TH1F(sname.str().c_str(),sname.str().c_str(),50,0,0.05);
                          std::vector<double> ySim,xSim,exSim,eySim;
                          for(int iSim = 1 ; iSim < spectrum2dSimDOIplot_1->GetNbinsX()  ; iSim++)
                          {
                            ySim.push_back(spectrum2dSimDOIplot_1->GetBinCenter(iSim));
                            xSim.push_back(spectrum2dSimDOIplot_1->GetBinContent(iSim));
                            eySim.push_back(0);
                            exSim.push_back(spectrum2dSimDOIplot_2->GetBinContent(iSim));
                            sigmaSim->Fill(spectrum2dSimDOIplot_2->GetBinContent(iSim));
                          }
                          TGraphErrors *gTot = new TGraphErrors(xSim.size(),&xSim[0],&ySim[0],&exSim[0],&eySim[0]);
                          sname.str("");
                          CurrentCrystal->SetSimZvsW(gTot);
                          CurrentCrystal->SetSimSigmaW(sigmaSim);






                        }
                      }
                    }
                    else
                    {
                      crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)]->SetCrystalOn(false);
                    }
                  }
                }
              }
            } // --end of OnForModular
          }  //--end of MPPC loop j
        }    //--end of MPPC loop i
      }
    }
    //----------------------------------------------------------//
    */
















  //----------------------------------------------------------//
  // SAVE DATA TO FILE                                        //
  //----------------------------------------------------------//



  std::cout << "Saving data to " << outputFileName << " ..." << std::endl;

  //----------------------------------------------------------//
  // Prepare Canvases                                         //
  //----------------------------------------------------------//
  TCanvas* C_spectrum;
  TCanvas* C_multi;
  TCanvas* C_multi_2d;
  // TCanvas* C_TaggingCrystalSpectrum = new TCanvas("TaggingCrystalSpectrum","TaggingCrystalSpectrum",1200,800);
  // TCanvas* RawCanvas = new TCanvas("RawSpectra","Rawspectra",1200,800);
  // TCanvas* TriggerCanvas = new TCanvas("TriggerSpectra","TriggerSpectra",1200,800);
  // TCanvas* BigLYCanvas = new TCanvas("BigLY","",800,800);
  // // TCanvas* FloodHistoCanvas = new TCanvas("FloodHisto","FloodHisto",800,800);
  // // TCanvas* FloodHisto3DCanvas = new TCanvas("FloodHisto3D","FloodHisto3D",800,800);
  // TCanvas* GlobalFlood2D = new TCanvas("Flood Histogram 2D","Flood Histogram 2D",800,800);
  // // TCanvas* GlobalFlood2DClean = new TCanvas("Flood Histogram 2D Clean","",800,800);
  // TCanvas* GlobalFlood3D = new TCanvas("Flood Histogram 3D","Flood Histogram 3D",800,800);
  // RawCanvas->Divide(nmppcx,nmppcy);
  // TriggerCanvas->Divide(nmppcx,nmppcy);
  // // plot general 2d and 3d plots of module, fill and draw multicanvas of raw and trigger spectra
  // for(int iModule = 0; iModule < nmodulex ; iModule++)
  // {
  //   for(int jModule = 0; jModule < nmoduley ; jModule++)
  //   {
  //     // GlobalFlood2DClean->cd();
  //     // module[iModule][jModule]->GetFloodMap2D()->SetTitle(""); // FIXME temporary removed title, for the poster...
  //     // module[iModule][jModule]->GetFloodMap2D()->Draw("COLZ");
  //     GlobalFlood2D->cd();
  //     module[iModule][jModule]->GetFloodMap2D()->Draw("COLZ");
  //     GlobalFlood3D->cd();
  //     module[iModule][jModule]->GetFloodMap3D()->Draw();
  //
  //     for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
  //     {
  //       for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
  //       {
  //         if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForModular())
  //         {
  //           RawCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
  //           RawCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition())->SetLogy();
  //           mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetRawSpectrum()->Draw();
  //
  //           TriggerCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
  //           TriggerCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition())->SetLogy();
  //           mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrum()->Draw();
  //           if( (thresholdKev > 0.0) || ( (backgroundRun) && (userBroadCut > 0))  )
  //           {
  //             mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrumHighlighted()->SetFillColor(3);
  //             mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrumHighlighted()->Draw("same");
  //           }
  //         }
  //       }
  //     }
  //   }
  // }


  //----------------------------------------------------------//
  // Produce some Canvases                                    //
  //----------------------------------------------------------//
  TCanvas* RawCanvas = new TCanvas("RawSpectra","Rawspectra",1200,800);
  TCanvas* TriggerCanvas = new TCanvas("TriggerSpectra","TriggerSpectra",1200,800);
  TCanvas* FloodHistoCanvas = new TCanvas("FloodHisto","FloodHisto",800,800);
  TCanvas* FloodHisto3DCanvas = new TCanvas("FloodHisto3D","FloodHisto3D",800,800);
  RawCanvas->Divide(nmppcx,nmppcy);
  TriggerCanvas->Divide(nmppcx,nmppcy);
  FloodHistoCanvas->Divide(nmppcx,nmppcy);
  FloodHisto3DCanvas->Divide(nmppcx,nmppcy);
  //canvases for the global plots
  TCanvas* C_TaggingCrystalSpectrum= new TCanvas("TaggingCrystalSpectrum","TaggingCrystalSpectrum",1200,800);
  TCanvas* GlobalFlood2D = new TCanvas("Flood Histogram 2D","Flood Histogram 2D",800,800);
  TCanvas* GlobalFlood2DClean = new TCanvas("Flood Histogram 2D Clean","",800,800);
  TCanvas* GlobalFlood3D = new TCanvas("Flood Histogram 3D","Flood Histogram 3D",800,800);
  TCanvas* BigSpectraCanvas = new TCanvas("BigSpectra","",800,800);
  TCanvas* BigLYCanvas = new TCanvas("BigLY","",800,800);
  BigSpectraCanvas->Divide(ncrystalsx*nmppcx*nmodulex,ncrystalsy*nmppcy*nmoduley);
  BigLYCanvas->Divide(ncrystalsx*nmppcx*nmodulex,ncrystalsy*nmppcy*nmoduley);

  // fill and draw multi canvas of BIG SPECTRA
  int canvascounter = 1;
  for(int jCrystal = ncrystalsy*nmppcy*nmoduley -1 ; jCrystal >=0  ; jCrystal--)
  {
    for(int iCrystal = 0 ; iCrystal < ncrystalsx*nmppcx*nmodulex   ; iCrystal++)
    {
      BigSpectraCanvas->cd(canvascounter);
      if(crystal[iCrystal][jCrystal]->GetIsOnForModular())
      {
        if(crystal[iCrystal][jCrystal]->CrystalIsOn())
        {
          if(saturationRun)
          {
            if(backgroundSaturationRun)
            {
              crystal[iCrystal][jCrystal]->GetInvestigatedSpectrum()->SetFillStyle(3001);
              crystal[iCrystal][jCrystal]->GetInvestigatedSpectrum()->SetFillColor(kBlue);
              crystal[iCrystal][jCrystal]->GetInvestigatedSpectrum()->Draw();
              if(performSaturationPeakSearch)
              {
                std::vector<TF1*> saturationFits = crystal[iCrystal][jCrystal]->GetSaturationFits();
                for(unsigned int iSaturation = 0; iSaturation < saturationFits.size(); iSaturation++)
                {
                  saturationFits[iSaturation]->Draw("same");
                }
              }
            }
            else
            {
              crystal[iCrystal][jCrystal]->GetSingleChargeSpectrum()->SetFillStyle(3001);
              crystal[iCrystal][jCrystal]->GetSingleChargeSpectrum()->SetFillColor(kBlue);
              crystal[iCrystal][jCrystal]->GetSingleChargeSpectrum()->Draw();
              if(performSaturationPeakSearch)
              {
                std::vector<TF1*> saturationFits = crystal[iCrystal][jCrystal]->GetSaturationFits();
                for(unsigned int iSaturation = 0; iSaturation < saturationFits.size(); iSaturation++)
                {
                  saturationFits[iSaturation]->Draw("same");
                }
              }
            }

          }
          else
          {
            if(backgroundRun)
            {
              crystal[iCrystal][jCrystal]->GetSpectrum()->SetFillStyle(3001);
              crystal[iCrystal][jCrystal]->GetSpectrum()->SetFillColor(kBlue);
              crystal[iCrystal][jCrystal]->GetSpectrum()->Draw();

            }
            else
            {
              if(correctingForDOI)
              {
                crystal[iCrystal][jCrystal]->GetCorrectedSpectrum()->SetFillStyle(3001);
                crystal[iCrystal][jCrystal]->GetCorrectedSpectrum()->SetFillColor(kBlue);
                crystal[iCrystal][jCrystal]->GetCorrectedSpectrum()->Draw();

                crystal[iCrystal][jCrystal]->GetHighlightedSpectrumCorrected()->SetFillStyle(3001);
                crystal[iCrystal][jCrystal]->GetHighlightedSpectrumCorrected()->SetFillColor(kGreen);
                crystal[iCrystal][jCrystal]->GetHighlightedSpectrumCorrected()->Draw("same");


              }
              else
              {
                crystal[iCrystal][jCrystal]->GetSpectrum()->SetFillStyle(3001);
                crystal[iCrystal][jCrystal]->GetSpectrum()->SetFillColor(kBlue);
                crystal[iCrystal][jCrystal]->GetSpectrum()->Draw();
              }
            }
          }
        }
      }
      BigLYCanvas->cd(canvascounter);
      if(lightYieldComputation)
      {
        crystal[iCrystal][jCrystal]->GetLYSpectrum()->SetFillStyle(3001);
        crystal[iCrystal][jCrystal]->GetLYSpectrum()->SetFillColor(kBlue);
        crystal[iCrystal][jCrystal]->GetLYSpectrum()->Draw();
      }
      canvascounter++;
    }
  }



  // plot general 2d and 3d plots of module, fill and draw multicanvas of raw and trigger spectra
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      // GlobalFlood2DClean->cd();
      // module[iModule][jModule]->GetFloodMap2D()->SetTitle(""); // FIXME temporary removed title, for the poster...
      // module[iModule][jModule]->GetFloodMap2D()->Draw("COLZ");
      GlobalFlood2D->cd();
      module[iModule][jModule]->GetFloodMap2D()->Draw("COLZ");
      GlobalFlood3D->cd();
      module[iModule][jModule]->GetFloodMap3D()->Draw();

      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
        for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
        {
          if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForModular())
          {
            RawCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
            RawCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition())->SetLogy();
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetRawSpectrum()->Draw();

            TriggerCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
            TriggerCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition())->SetLogy();
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrum()->Draw();
            if( (thresholdKev > 0.0) || ( (backgroundRun) && (userBroadCut > 0))  )
            {
              mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrumHighlighted()->SetFillColor(3);
              mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrumHighlighted()->Draw("same");
            }
          }
        }
      }
    }
  }
  //----------------------------------------------------------//



  // //----------------------------------------------------------//
  // // Global distributions                                     //
  // //----------------------------------------------------------//
  //plots to summarize the values of relevant variable found on each crystal
  //--Distribution of photopeak positions, in ADC channels
  //histogram
  TH1F *PeakPositionDistro = new TH1F("Photopeak position","Distribution photopeak positions",100,0,12000);
  PeakPositionDistro->GetXaxis()->SetTitle("ADC Channels");
  PeakPositionDistro->GetYaxis()->SetTitle("N");
  PeakPositionDistro->SetStats(1);
  //2d histogram

  TH1F *LYDistro = new TH1F("Light Yield position","Light Yield Distribution - All Crystals",250,0,40000);
  LYDistro->GetXaxis()->SetTitle("[Ph/MeV]");
  LYDistro->GetYaxis()->SetTitle("N");
  LYDistro->SetStats(1);

  TH1F *LYDistroCentral = new TH1F("Light Yield position - Central","Light Yield Distribution - Central Crystals",250,0,40000);
  LYDistroCentral->GetXaxis()->SetTitle("[Ph/MeV]");
  LYDistroCentral->GetYaxis()->SetTitle("N");
  LYDistroCentral->SetStats(1);

  //histogram
  TH1F *PeakEnergyResolutionDistro = new TH1F("Energy res FWHM","Distribution photopeak energy resolutions FWHM",200,0,1);
  PeakEnergyResolutionDistro->GetXaxis()->SetTitle("Energy Resolution FWHM");
  PeakEnergyResolutionDistro->GetYaxis()->SetTitle("N");
  PeakEnergyResolutionDistro->SetStats(1);

  TH1F *PeakPositionDistroCentral = new TH1F("Central Photopeak Position","Distribution photopeak positions - Central Crystals",200,0,histo1Dmax);
  PeakPositionDistroCentral->GetXaxis()->SetTitle("ADC Channels");
  PeakPositionDistroCentral->GetYaxis()->SetTitle("N");
  PeakPositionDistroCentral->SetStats(1);

  TH1F *PeakEnergyResolutionDistroCentral = new TH1F("Central Energy Resolution","Distribution photopeak energy resolutions FWHM - Central Crystals",200,0,1);
  PeakEnergyResolutionDistroCentral->GetXaxis()->SetTitle("Energy Resolution FWHM");
  PeakEnergyResolutionDistroCentral->GetYaxis()->SetTitle("N");
  PeakEnergyResolutionDistroCentral->SetStats(1);


  TH1F *PeakEnergyResolutionDistro_corr = new TH1F("Corrected Energy res FWHM","Distribution photopeak energy resolutions FWHM",200,0,1);
  PeakEnergyResolutionDistro_corr->GetXaxis()->SetTitle("Energy Resolution FWHM");
  PeakEnergyResolutionDistro_corr->GetYaxis()->SetTitle("N");
  PeakEnergyResolutionDistro_corr->SetStats(1);


  TH1F *PeakEnergyResolutionDistroCentral_corr = new TH1F("Corrected Central En Res","Distribution photopeak energy resolutions FWHM - Central Crystals",200,0,1);
  PeakEnergyResolutionDistroCentral_corr->GetXaxis()->SetTitle("Energy Resolution FWHM");
  PeakEnergyResolutionDistroCentral_corr->GetYaxis()->SetTitle("N");
  PeakEnergyResolutionDistroCentral_corr->SetStats(1);

  //LEGO PLOTS RELEVANT
  TH2F *PeakPositionVsIJ = new TH2F("Photopeak positions vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  PeakPositionVsIJ->GetXaxis()->SetTitle("i (U axis)");
  PeakPositionVsIJ->GetYaxis()->SetTitle("j (V axis)");
  PeakPositionVsIJ->GetZaxis()->SetTitle("ADC Channels");
  PeakPositionVsIJ->GetXaxis()->SetTitleOffset(1.8);
  PeakPositionVsIJ->GetYaxis()->SetTitleOffset(1.8);
  PeakPositionVsIJ->GetZaxis()->SetTitleOffset(2.2);
  PeakPositionVsIJ->GetZaxis()->SetRangeUser(0,PeakPositionVsIJmax);
  //2d histogram
  TH2F *EnergyResolutionVsIJ = new TH2F("Energy res FWHM vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  EnergyResolutionVsIJ->GetXaxis()->SetTitle("i (U axis)");
  EnergyResolutionVsIJ->GetYaxis()->SetTitle("j (V axis)");
  EnergyResolutionVsIJ->GetZaxis()->SetTitle("Energy Resolution FWHM");
  EnergyResolutionVsIJ->GetXaxis()->SetTitleOffset(1.8);
  EnergyResolutionVsIJ->GetYaxis()->SetTitleOffset(1.8);
  EnergyResolutionVsIJ->GetZaxis()->SetTitleOffset(2.2);
  EnergyResolutionVsIJ->GetZaxis()->SetRangeUser(0,EnergyResolutionVsIJmax);

  //LYVsIJ
  TH2F *LYVsIJ = new TH2F("Light Yield vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  LYVsIJ->GetXaxis()->SetTitle("i (U axis)");
  LYVsIJ->GetYaxis()->SetTitle("j (V axis)");
  LYVsIJ->GetZaxis()->SetTitle("Light Yield [Ph/MeV]");
  LYVsIJ->GetXaxis()->SetTitleOffset(1.8);
  LYVsIJ->GetYaxis()->SetTitleOffset(1.8);
  LYVsIJ->GetZaxis()->SetTitleOffset(2.2);
  LYVsIJ->GetZaxis()->SetRangeUser(0,LYvsIJmax);

  TH2F *EnergyResolutionVsIJ_corr = new TH2F("Corrected Energy res FWHM vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  EnergyResolutionVsIJ_corr->GetXaxis()->SetTitle("i (U axis)");
  EnergyResolutionVsIJ_corr->GetYaxis()->SetTitle("j (V axis)");
  EnergyResolutionVsIJ_corr->GetZaxis()->SetTitle("Energy Resolution FWHM");
  EnergyResolutionVsIJ_corr->GetXaxis()->SetTitleOffset(1.8);
  EnergyResolutionVsIJ_corr->GetYaxis()->SetTitleOffset(1.8);
  EnergyResolutionVsIJ_corr->GetZaxis()->SetTitleOffset(2.2);
  EnergyResolutionVsIJ_corr->GetZaxis()->SetRangeUser(0,0.3);
  //----------------------------------------------------------//


  //----------------------------------------------------------//
  // Write output to root file(s)                             //
  //----------------------------------------------------------//
  //write output plots
  TFile* fPlots = new TFile(outputFileName.c_str(),"recreate");
  fPlots->cd();
  //TDirectory
  TDirectory ****directory;
  directory = new TDirectory***[nmodulex*nmoduley];
  for(int i = 0; i < nmodulex*nmoduley+1 ; i++)
  {
    directory[i] = new TDirectory** [nmppcx*nmppcy+1];
    for(int j = 0; j < nmppcx*nmppcy+1 ; j++)
    {
      directory[i][j] = new TDirectory* [(nmppcx*ncrystalsx)*(nmppcy*ncrystalsy)+1];
    }
  }

  //run on all the elements, make directories in the output file, save the plots
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      std::stringstream ModuleDirStream;
      ModuleDirStream << "Module " << iModule << "." << jModule;
      ModuleDirStream.str();
      directory[iModule+jModule][0][0] = fPlots->mkdir(ModuleDirStream.str().c_str());
      directory[iModule+jModule][0][0]->cd();
      GlobalFlood2D->Write();
      GlobalFlood3D->Write();
      //
      // if(usingRealSimData)
      // {
      //   C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
      //   C_spectrum->SetName(module[iModule][jModule]->GetFloodMap2DSingleCrystalHit()->GetName());
      //   C_spectrum->cd();
      //   module[iModule][jModule]->GetFloodMap2DSingleCrystalHit()->Draw("COLZ");
      //   C_spectrum->Write();
      //   delete C_spectrum;
      //
      // }
      //
      if(usingTaggingBench || calcDoiResWithCalibration || taggingForTiming)
      {
        C_TaggingCrystalSpectrum->cd();
        TaggingCrystalSpectrum->Draw();
        TriggerSpectrumHighlight->Draw("same");
        C_TaggingCrystalSpectrum->Write();

        std::stringstream TagNum;
        TagNum.str("");
        TagNum << module[iModule][jModule]->GetTaggingTimingChannel();
        TNamed tNum("taggingCrystalTimingChannel",TagNum.str().c_str());
        tNum.Write();
        TagNum.str("");
      }
      //
      RawCanvas->Write();
      TriggerCanvas->Write();
      BigSpectraCanvas->Write();
      if(usingTaggingBench || calcDoiResWithCalibration || taggingForTiming) module[iModule][jModule]->GetTaggingPhotopeakCut().Write();
      if(lightYieldComputation) BigLYCanvas->Write();
      //
      if(smearTaggingTime)
      {
        std::stringstream sseed;
        sseed << module[iModule][jModule]->GetSeed();
        TNamed SeedD("Seed",sseed.str().c_str());
        SeedD.Write();
      }

      std::vector<int> detChData = module[iModule][jModule]->GetDetChData();
      std::vector<float> saturationData = module[iModule][jModule]->GetSaturationData();
      std::vector<float> pedestalData = module[iModule][jModule]->GetPedestalData();
      gDirectory->WriteObject(&detChData, "channels");
      gDirectory->WriteObject(&saturationData, "saturation");
      gDirectory->WriteObject(&pedestalData, "pedestal");
      //
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
        for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
        {
          if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForModular())
          {
            //directory
            std::stringstream MppcDirStream;
            MppcDirStream << "MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " - " <<  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetExtendedID();
            directory[iModule+jModule][(iMppc+jMppc)+1][0] = directory[iModule+jModule][0][0]->mkdir(MppcDirStream.str().c_str());
            directory[iModule+jModule][(iMppc+jMppc)+1][0]->cd();

            //2d flood map
            C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
            C_spectrum->SetName(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->GetName());
            C_spectrum->cd();
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->Draw("COLZ");
            C_spectrum->Write();
            delete C_spectrum;

            // 3d flood map
            C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
            C_spectrum->SetName(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap3D()->GetName());
            C_spectrum->cd();
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap3D()->Draw();
            C_spectrum->Write();
            delete C_spectrum;


           // C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
           // C_spectrum->SetName(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerChargeSpectrum()->GetName());
           // C_spectrum->cd();
           // mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerChargeSpectrum()->Draw();
           // C_spectrum->Write();
           // delete C_spectrum;
           //
           // C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
           // C_spectrum->SetName(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetRawChargeSpectrum()->GetName());
           // C_spectrum->cd();
           // mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetRawChargeSpectrum()->Draw();
           // C_spectrum->Write();
           // delete C_spectrum;

           //write digitizerChannel for this crystal
           // std::stringstream sChNum;
           // sChNum << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
           // TNamed ChNum("digitizerChannel",sChNum.str().c_str());
           // ChNum.Write();
           //
           // sChNum.str("");
           // sChNum << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetSaturation();
           // TNamed satTnamed("saturation",sChNum.str().c_str());
           // satTnamed.Write();


            //prepare a Canvas for the nxn 3D cuts
            C_multi = new TCanvas("C_multi","C_multi",1200,1200);
            TString nameMppc = "3D Cuts - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            C_multi->SetName(nameMppc);
            C_multi->SetTitle(nameMppc);
            TLegend *legend = new TLegend(0.7,0.75,0.893,0.89,"");
            legend->SetFillStyle(0);
            int counter = 1;

            C_multi_2d = new TCanvas("C_multi_2d","C_multi_2d",1200,1200);
            TString nameMppc_2d = "3D Cuts in 2D plots - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            C_multi_2d->SetName(nameMppc_2d);
            C_multi_2d->SetTitle(nameMppc_2d);
            TLegend *legend_2d = new TLegend(0.7,0.75,0.893,0.89,"");
            legend_2d->SetFillStyle(0);


            for(int iCry = 0; iCry < ncrystalsx ; iCry++)
            {
              for(int jCry = 0; jCry < ncrystalsy ; jCry++)
              {

                //directory
                std::stringstream CrystalDirStream;
                //create a pointer for the current crystal (mainly to make the code more readable)
                Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
                if(CurrentCrystal->GetIsOnForModular())
                {
                  CrystalDirStream << "Crystal " <<  CurrentCrystal->GetID();
                  directory[iModule+jModule][(iMppc+jMppc)+1][(iCry+jCry)+1] = directory[iModule+jModule][(iMppc+jMppc)+1][0]->mkdir(CrystalDirStream.str().c_str());
                  directory[iModule+jModule][(iMppc+jMppc)+1][(iCry+jCry)+1]->cd();
                  if(CurrentCrystal->CrystalIsOn() /*| usingRealSimData*/) // save data only if the crystal was found
                  {

                    //fill the global distributions histograms
                    //draw also the 3d cuts in the common canvas of this mppc
                    C_multi->cd();
                    CurrentCrystal->GetFloodMap3D()->SetMarkerColor(counter);
                    CurrentCrystal->GetFloodMap3D()->SetFillColor(counter);
                    legend->AddEntry(CurrentCrystal->GetFloodMap3D(),CurrentCrystal->GetFloodMap3D()->GetName(),"f");
                    CurrentCrystal->GetFloodMap3D()->Draw("same");
                    legend->Draw();
                    // 		counter++;

                    // same as above, but the 2D maxPeak
                    C_multi_2d->cd();
                    CurrentCrystal->GetFloodMap2D()->SetMarkerColor(counter);
                    CurrentCrystal->GetFloodMap2D()->SetFillColor(counter);
                    legend_2d->AddEntry(CurrentCrystal->GetFloodMap2D(),CurrentCrystal->GetFloodMap2D()->GetName(),"f");
                    if(counter == 1)
                      CurrentCrystal->GetFloodMap2D()->Draw();
                    else
                      CurrentCrystal->GetFloodMap2D()->Draw("same");
                    legend_2d->Draw();
                    counter++;
      //
      //
                    if(saturationRun)
                    {

                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName(CurrentCrystal->GetSingleChargeSpectrum()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetSingleChargeSpectrum()->Draw();
                      C_spectrum->Write();
                      delete C_spectrum;

                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName(CurrentCrystal->GetSumChargeSpectrum()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetSumChargeSpectrum()->Draw();
                      C_spectrum->Write();
                      delete C_spectrum;

                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName(CurrentCrystal->GetAMTChargeSpectrum()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetAMTChargeSpectrum()->Draw();
                      C_spectrum->Write();
                      delete C_spectrum;

                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName(CurrentCrystal->GetInvestigatedSpectrum()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetInvestigatedSpectrum()->Draw();
                      std::vector<TF1*> saturationFits = CurrentCrystal->GetSaturationFits();
                      for(unsigned int iSaturation = 0; iSaturation < saturationFits.size(); iSaturation++)
                      {
                        saturationFits[iSaturation]->Draw("same");
                      }
                      C_spectrum->Write();
                      delete C_spectrum;


                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName(CurrentCrystal->GetNotNeighboursSingleCharge()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetNotNeighboursSingleCharge()->Draw();
                      C_spectrum->Write();
                      delete C_spectrum;

                    }
                    else
                    {
                      if(!backgroundRun)
                      {
                        PeakPositionDistro->Fill(CurrentCrystal->GetPhotopeakPosition());
                        PeakEnergyResolutionDistro->Fill(CurrentCrystal->GetPhotopeakEnergyResolution());

                        if(correctingForDOI)
                        {
                          PeakEnergyResolutionDistro_corr->Fill(CurrentCrystal->GetPhotopeakEnergyResolutionCorrected());
                        }

                        PeakPositionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakPosition());
                        EnergyResolutionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakEnergyResolution());
                        EnergyResolutionVsIJ_corr->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakEnergyResolutionCorrected());
                      }

                      if(lightYieldComputation)
                      {
                        LYDistro->Fill(CurrentCrystal->GetLY());
                        LYVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetLY());
                      }
                      //
                      //
                      //
                      //
                      if( ((iModule*nmppcx)+iMppc) > 0 && (((iModule*nmppcx)+iMppc) < nmppcx -1) && ((jModule*nmppcy)+jMppc) > 0 && (((jModule*nmppcy)+jMppc) < nmppcy -1 ))
                      {
                        if(!backgroundRun)
                        {
                          PeakPositionDistroCentral->Fill(CurrentCrystal->GetPhotopeakPosition());
                          PeakEnergyResolutionDistroCentral->Fill(CurrentCrystal->GetPhotopeakEnergyResolution());

                          if(correctingForDOI)
                          {
                            PeakEnergyResolutionDistroCentral_corr->Fill(CurrentCrystal->GetPhotopeakEnergyResolutionCorrected());
                          }
                        }
                        if(lightYieldComputation)
                        {
                          LYDistroCentral->Fill(CurrentCrystal->GetLY());
                        }
                      }
                      //
                      // spectrum
                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName(CurrentCrystal->GetSpectrum()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetSpectrum()->Draw();
                      if(!backgroundRun)
                      {
                        CurrentCrystal->GetHighlightedSpectrum()->SetFillColor(3);
                        CurrentCrystal->GetHighlightedSpectrum()->Draw("same");
                        CurrentCrystal->GetFit()->Draw("same");
                      }
                      C_spectrum->Write();
                      delete C_spectrum;

                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName("Light Sharing");
                      C_spectrum->Divide(nmppcx,nmppcy);
                      //first plot the Delta Tcry - Ttagging for this crystal
                      C_spectrum->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()) ;
                      CurrentCrystal->GetLScentralSpectrum()->Draw();
                      for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetLSSpectra().size() ; iNeig++)
                      {
                        C_spectrum->cd(CurrentCrystal->GetLSSpectra()[iNeig].canvasPosition);
                        CurrentCrystal->GetLSSpectra()[iNeig].spectrum->Draw();
                      }
                      C_spectrum->Write();
                      delete C_spectrum;
                      //
                      if(lightYieldComputation)
                      {
                        //LY spectrum
                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                        C_spectrum->SetName(CurrentCrystal->GetLYSpectrum()->GetName());
                        C_spectrum->cd();
                        CurrentCrystal->GetLYSpectrum()->Draw();
                        CurrentCrystal->GetLYFit()->Draw("same");
                        C_spectrum->Write();
                        delete C_spectrum;
                      }

                      //w histo
                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName(CurrentCrystal->GetHistoW()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetHistoW()->Draw();
                      if(usingTaggingBench || calcDoiResWithCalibration) CurrentCrystal->GetHistoWfit()->Draw("same");
                      C_spectrum->Write();
                      delete C_spectrum;

                      // adc versus time
                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName(CurrentCrystal->GetVersusTime()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetVersusTime()->Draw("COLZ");
                      C_spectrum->Write();
                      delete C_spectrum;

                      //W versus time
                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName(CurrentCrystal->GetWversusTime()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetWversusTime()->Draw("COLZ");
                      C_spectrum->Write();
                      delete C_spectrum;

                      // adc versus w
                      if(!backgroundRun)
                      {
                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
                        C_spectrum->SetName(CurrentCrystal->GetADCversusW()->GetName());
                        C_spectrum->cd();
                        CurrentCrystal->GetADCversusW()->Draw("COLZ");
                        C_spectrum->Write();
                        delete C_spectrum;
                      }
                      // adc versus w complete
                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
                      C_spectrum->SetName(CurrentCrystal->GetADCversusWComplete()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetADCversusWComplete()->Draw("COLZ");
                      C_spectrum->Write();
                      delete C_spectrum;

                      if(!backgroundRun)
                      {
                        if(correctingForDOI)
                        {
                          C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                          std::stringstream sname;
                          sname << "Energy vs. W Correction - Crystal " <<  CurrentCrystal->GetID();
                          C_spectrum->SetName(sname.str().c_str());
                          C_spectrum->cd();
                          CurrentCrystal->GetSlicesMean()->Draw();
                          CurrentCrystal->GetSlicesMeanFit()->Draw("same");
                          C_spectrum->Write();
                          delete C_spectrum;

                          C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                          C_spectrum->SetName(CurrentCrystal->GetCorrectedSpectrum()->GetName());
                          C_spectrum->cd();
                          CurrentCrystal->GetCorrectedSpectrum()->Draw();
                          CurrentCrystal->GetHighlightedSpectrumCorrected()->SetFillColor(3);
                          CurrentCrystal->GetHighlightedSpectrumCorrected()->Draw("same");
                          CurrentCrystal->GetFitCorrected()->Draw("same");
                          C_spectrum->Write();
                          delete C_spectrum;

                        }
                      }

                      //w histo corrected
                      C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                      C_spectrum->SetName(CurrentCrystal->GetHistoWCorrected()->GetName());
                      C_spectrum->cd();
                      CurrentCrystal->GetHistoWCorrected()->Draw();
                      C_spectrum->Write();
                      delete C_spectrum;


                      if(usingTaggingBench || taggingForTiming)
                      {

                        if(digitizerType == 1 || digitizerType == 2)
                        {
                          //write timingChannel for this crystal
                          std::stringstream sChNum;
                          sChNum.str("");
                          sChNum << CurrentCrystal->GetTimingChannel();
                          TNamed tCryNum("timingChannel",sChNum.str().c_str());
                          tCryNum.Write();
                          sChNum.str("");


                          if(timingCorrectionForPolished) // only if timing correction is performed
                          {


                            C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                            C_spectrum->SetName("Delta Tcry - TNeig");
                            C_spectrum->Divide(nmppcx,nmppcy);
                            //first plot the Delta Tcry - Ttagging for this crystal
                            C_spectrum->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()) ;
                            CurrentCrystal->GetDeltaTimeWRTTagging()->Draw();
                            for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetDeltaTcryTneig().size() ; iNeig++)
                            {
                              C_spectrum->cd(CurrentCrystal->GetDeltaTcryTneig()[iNeig].canvasPosition);
                              CurrentCrystal->GetDeltaTcryTneig()[iNeig].spectrum->Draw();
                            }
                            C_spectrum->Write();
                            delete C_spectrum;

                            C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                            C_spectrum->SetName("Delta Tcry - TNeig vs W");
                            C_spectrum->Divide(nmppcx,nmppcy);
                            //first plot the Delta Tcry - Ttagging for this crystal
                            C_spectrum->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()) ;
                            CurrentCrystal->GetDeltaTvsW()->Draw("COLZ");
                            for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetDeltaT2vsW().size() ; iNeig++)
                            {
                              C_spectrum->cd(CurrentCrystal->GetDeltaT2vsW()[iNeig].canvasPosition);
                              CurrentCrystal->GetDeltaT2vsW()[iNeig].spectrum->Draw("COLZ");
                            }
                            C_spectrum->Write();
                            delete C_spectrum;

                            C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                            C_spectrum->SetName("Delta Tcry - TNeig vs CH");
                            C_spectrum->Divide(nmppcx,nmppcy);
                            //first plot the Delta Tcry - Ttagging for this crystal
                            C_spectrum->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()) ;
                            CurrentCrystal->GetDeltaTvsCH()->Draw("COLZ");
                            for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetDeltaT2vsCH().size() ; iNeig++)
                            {
                              C_spectrum->cd(CurrentCrystal->GetDeltaT2vsCH()[iNeig].canvasPosition);
                              CurrentCrystal->GetDeltaT2vsCH()[iNeig].spectrum->Draw("COLZ");
                            }
                            C_spectrum->Write();
                            delete C_spectrum;

                            if(timingCorrection)
                            {
                              C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                              C_spectrum->SetName("Correction Graphs");
                              C_spectrum->Divide(nmppcx,nmppcy);
                              C_spectrum->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()) ;
                              CurrentCrystal->GetGraphDeltaW()->Draw("AL");
                              for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetGraphDelayW().size() ; iNeig++)
                              {
                                C_spectrum->cd(CurrentCrystal->GetGraphDelayW()[iNeig].canvasPosition);
                                CurrentCrystal->GetGraphDelayW()[iNeig].spectrum->Draw("AL");
                              }
                              C_spectrum->Write();
                              delete C_spectrum;




                              C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                              C_spectrum->SetName("RMS Correction Graphs");
                              C_spectrum->Divide(nmppcx,nmppcy);
                              C_spectrum->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()) ;
                              CurrentCrystal->GetGraphDeltaRMS()->Draw("AL");
                              for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetGraphDelayRMS().size() ; iNeig++)
                              {
                                C_spectrum->cd(CurrentCrystal->GetGraphDelayRMS()[iNeig].canvasPosition);
                                CurrentCrystal->GetGraphDelayRMS()[iNeig].spectrum->Draw("AL");
                              }
                              C_spectrum->Write();
                              delete C_spectrum;

                              //create TimeCorrection folder
                              TDirectory *corrDir = directory[iModule+jModule][(iMppc+jMppc)+1][(iCry+jCry)+1]->mkdir("TimeCorrection");
                              corrDir->cd();

                              // //create the w subfolders
                              // for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetGraphDelayW().size() ; iNeig++)
                              // {
                              //
                              // }

                              //write delta graphs to file
                              CurrentCrystal->GetGraphDeltaW()->Write();
                              for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetGraphDelayW().size() ; iNeig++)
                              {
                                CurrentCrystal->GetGraphDelayW()[iNeig].spectrum->Write();
                              }

                              //write delta RMS graphs to file
                              CurrentCrystal->GetGraphDeltaRMS()->Write();
                              for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetGraphDelayW().size() ; iNeig++)
                              {
                                CurrentCrystal->GetGraphDelayRMS()[iNeig].spectrum->Write();
                              }

                              //save also the individual histograms
                              for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetDeltaTHistos().size() ; iNeig++)
                              {
                                CurrentCrystal->GetDeltaTHistos()[iNeig]->Write();
                              }

                              for(unsigned int iNeig = 0 ; iNeig < CurrentCrystal->GetDelayHistos().size() ; iNeig++)
                              {
                                CurrentCrystal->GetDelayHistos()[iNeig]->Write();
                              }



                              //write time Channels for neighbouring crystals
                              std::vector<int> DelayTimingChannelsNum = CurrentCrystal->GetDelayTimingChannels();
                              gDirectory->WriteObject(&DelayTimingChannelsNum, "delayTimingChannels");
                              corrDir->cd("..");
                            }

                            // std::vector<int>    GetTChannelsForPolishedCorrection(){return tChannelsForPolishedCorrection;};
                            // std::vector<double> GetMeanForPolishedCorrection(){return meanForPolishedCorrection;};
                            // std::vector<double> GetFwhmForPolishedCorrection(){return fwhmForPolishedCorrection;};

                            std::vector<int> tChannelsForPolishedCorrection = CurrentCrystal->GetTChannelsForPolishedCorrection();
                            gDirectory->WriteObject(&tChannelsForPolishedCorrection, "tChannelsForPolishedCorrection");
                            std::vector<double> meanForPolishedCorrection = CurrentCrystal->GetMeanForPolishedCorrection();
                            gDirectory->WriteObject(&meanForPolishedCorrection, "meanForPolishedCorrection");
                            std::vector<double> fwhmForPolishedCorrection = CurrentCrystal->GetFwhmForPolishedCorrection();
                            gDirectory->WriteObject(&fwhmForPolishedCorrection, "fwhmForPolishedCorrection");

                          }
                          else // otherwise plot only the CTR plot
                          {
                            C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                            C_spectrum->SetName(CurrentCrystal->GetDeltaTimeWRTTagging()->GetName());
                            CurrentCrystal->GetDeltaTimeWRTTagging()->Draw();
                            C_spectrum->Write();
                            delete C_spectrum;

                            C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
                            C_spectrum->SetName(CurrentCrystal->GetDeltaTvsW()->GetName());
                            CurrentCrystal->GetDeltaTvsW()->Draw("COLZ");
                            C_spectrum->Write();
                            delete C_spectrum;

                            C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
                            C_spectrum->SetName(CurrentCrystal->GetDeltaTvsCH()->GetName());
                            CurrentCrystal->GetDeltaTvsCH()->Draw("COLZ");
                            C_spectrum->Write();
                            delete C_spectrum;


                          }
                        }
                      }

                      if(!calcDoiResWithCalibration)
                      {
                        //pdf
                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                        C_spectrum->SetName(CurrentCrystal->GetPdfW()->GetName());
                        C_spectrum->cd();
                        CurrentCrystal->GetPdfW()->Draw();
                        C_spectrum->Write();
                        delete C_spectrum;

                        //cumulative
                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                        C_spectrum->SetName(CurrentCrystal->GetCumulativeW()->GetName());
                        C_spectrum->cd();
                        CurrentCrystal->GetCumulativeW()->Draw();
                        C_spectrum->Write();
                        delete C_spectrum;

                        //compton w(z)
                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                        C_spectrum->SetName(CurrentCrystal->GetWZgraph()->GetName());
                        C_spectrum->cd();
                        CurrentCrystal->GetWZgraph()->Draw("AL");
                        C_spectrum->Write();
                        delete C_spectrum;
                      }

                      //calibration graph
                      //this one needs to be checked, especially in case of calcDoiResWithCalibration but it doesn't hurt to always check..
                      TGraph *calibrationGraph = NULL;
                      calibrationGraph = CurrentCrystal->GetCalibrationGraph();
                      if(calibrationGraph)
                      {
                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                        C_spectrum->SetName(CurrentCrystal->GetCalibrationGraph()->GetName());
                        C_spectrum->cd();
                        CurrentCrystal->GetCalibrationGraph()->Draw("AL");
                        C_spectrum->Write();
                        delete C_spectrum;
                      }

                      CurrentCrystal->GetZYCut()->Write();
                      CurrentCrystal->GetZXCut()->Write();
                      CurrentCrystal->GetCrystalCut().Write();
                      CurrentCrystal->GetCrystalCutWithoutCutG().Write();
                      CurrentCrystal->GetPhotopeakEnergyCut().Write();

                      if(usingTaggingBench || calcDoiResWithCalibration)
                      {
                        if(calcDoiResWithDelta)
                        {
                          C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                          C_spectrum->SetName(CurrentCrystal->GetHistoAltDoiRes()->GetName());
                          C_spectrum->cd();
                          CurrentCrystal->GetHistoAltDoiRes()->Draw();
                          CurrentCrystal->GetHistoAltDoiFit()->Draw("same");
                          C_spectrum->Write();
                          delete C_spectrum;
                        }
                      }

                      //write digitizerChannel for this crystal
                      std::stringstream sChNum;
                      sChNum << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
                      TNamed ChNum("digitizerChannel",sChNum.str().c_str());
                      ChNum.Write();

                      //write charge channels used for W calculation for this crystal
                      std::vector<int> channelsNumRelevantForW = CurrentCrystal->GetRelevantForW();
                      gDirectory->WriteObject(&channelsNumRelevantForW, "channelsNumRelevantForW");





                      // // compton calibration
                      // if(comptonAnalysis)
                      // {
                      //   for(int iComptMppc = 0; iComptMppc < nmppcx ; iComptMppc++)
                      //   {
                      //     for(int jComptMppc = 0; jComptMppc < nmppcy ; jComptMppc++)
                      //     {
                      //       TGraph2D ***tempGraph = CurrentCrystal->GetComptonCalibration();
                      //       tempGraph[iComptMppc][jComptMppc]->Write();
                      //       TGraph2D ***tempCorrGraph = CurrentCrystal->GetConvertedComptonCalibration();
                      //       tempCorrGraph[iComptMppc][jComptMppc]->Write();
                      //     }
                      //   }
                      //   for(int iCompt = 0 ; iCompt < CurrentCrystal->GetNumOfComptonHisto(); iCompt++)
                      //   {
                      //     C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
                      //     C_spectrum->SetName(CurrentCrystal->GetComptonHistogram(iCompt)->GetName());
                      //     C_spectrum->cd();
                      //     CurrentCrystal->GetComptonHistogram(iCompt)->Draw();
                      //     C_spectrum->Write();
                      //     delete C_spectrum;
                      //   }
                      // }
                      //
                      if(usingRealSimData)
                      {

                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
                        C_spectrum->SetName(CurrentCrystal->GetSimDOIplot()->GetName());
                        C_spectrum->cd();
                        CurrentCrystal->GetSimDOIplot()->Draw("COLZ");
                        C_spectrum->Write();
                        delete C_spectrum;

                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                        std::stringstream sname;
                        sname << "Sim vs. Calibration - Crystal " <<  CurrentCrystal->GetID();
                        C_spectrum->SetName(sname.str().c_str());
                        C_spectrum->cd();
                        TLegend *legend_sim = new TLegend(0.6,0.75,0.893,0.89,"");
                        legend_sim->SetFillStyle(0);
                        CurrentCrystal->GetSimZvsW()->SetTitle(sname.str().c_str());
                        legend_sim->AddEntry(CurrentCrystal->GetSimZvsW(),"Real Z","l");
                        legend_sim->AddEntry(CurrentCrystal->GetCalibrationGraph(),"Calibration from W","l");
                        TMultiGraph* mg = new TMultiGraph();
                        CurrentCrystal->GetCalibrationGraph()->SetLineColor(kRed);
                        mg->Add(CurrentCrystal->GetSimZvsW(),"p");
                        mg->Add(CurrentCrystal->GetCalibrationGraph(),"l");
                        //                       CurrentCrystal->GetSimZvsW()->Draw();
                        //                       CurrentCrystal->GetCalibrationGraph()->Draw("same");
                        mg->Draw("a");
                        legend_sim->Draw();
                        C_spectrum->Write();
                        delete C_spectrum;

                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
                        C_spectrum->SetName(CurrentCrystal->GetSimSigmaW()->GetName());
                        C_spectrum->cd();
                        CurrentCrystal->GetSimSigmaW()->GetXaxis()->SetTitle("W");
                        CurrentCrystal->GetSimSigmaW()->Draw();
                        C_spectrum->Write();
                        delete C_spectrum;

                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                        C_spectrum->SetName(CurrentCrystal->GetSimGraph()->GetName());
                        C_spectrum->cd();
                        CurrentCrystal->GetSimGraph()->Draw("ap");
                        // 		      CurrentCrystal->GetSimFit()->Draw("same");
                        C_spectrum->Write();
                        delete C_spectrum;
                      }
                    }
                  }
                }
              }
            }




      //
            directory[iModule+jModule][(iMppc+jMppc)+1][0]->cd();
      //
            C_multi->Write();
            C_multi_2d->Write();
      //
            delete C_multi_2d;
            delete C_multi;
          }
        }
      }


      directory[iModule+jModule][0][0]->cd(); // go back to main directory


      if(!backgroundRun)
      {
        TCanvas *C_PeakPositionVsIJ = new TCanvas("C_PeakPositionVsIJ","C_PeakPositionVsIJ",800,800);
        C_PeakPositionVsIJ->SetName(PeakPositionVsIJ->GetName());
        C_PeakPositionVsIJ->cd();
        PeakPositionVsIJ->Draw("LEGO2");
        C_PeakPositionVsIJ->SetLeftMargin(0.15);
        C_PeakPositionVsIJ->Write();

        if(correctingForDOI)
        {
          TCanvas *C_EnergyResolutionVsIJ_corr = new TCanvas("C_EnergyResolutionVsIJ_corr","C_EnergyResolutionVsIJ_corr",800,800);
          C_EnergyResolutionVsIJ_corr->SetName(EnergyResolutionVsIJ_corr->GetName());
          C_EnergyResolutionVsIJ_corr->cd();
          EnergyResolutionVsIJ_corr->Draw("LEGO2");
          C_EnergyResolutionVsIJ_corr->SetLeftMargin(0.15);
          C_EnergyResolutionVsIJ_corr->Write();
        }

        TCanvas *C_EnergyResolutionVsIJ = new TCanvas("C_EnergyResolutionVsIJ","C_EnergyResolutionVsIJ",800,800);
        C_EnergyResolutionVsIJ->SetName(EnergyResolutionVsIJ->GetName());
        C_EnergyResolutionVsIJ->cd();
        EnergyResolutionVsIJ->Draw("LEGO2");
        C_EnergyResolutionVsIJ->SetLeftMargin(0.15);
        C_EnergyResolutionVsIJ->Write();
      }

      if(lightYieldComputation)
      {
        TCanvas *C_LYVsIJ = new TCanvas("C_LYVsIJ","C_LYVsIJ",800,800);
        C_LYVsIJ->SetName(LYVsIJ->GetName());
        C_LYVsIJ->cd();
        LYVsIJ->Draw("LEGO2");
        C_LYVsIJ->SetLeftMargin(0.15);
        C_LYVsIJ->Write();
      }

      if(!backgroundRun)
      {
        PeakPositionDistro->Write();
        PeakPositionDistroCentral->Write();
        if(correctingForDOI)
        {
          PeakEnergyResolutionDistro_corr->Write();
          PeakEnergyResolutionDistroCentral_corr->Write();
        }
        PeakEnergyResolutionDistro->Write();
        PeakEnergyResolutionDistroCentral->Write();

      }

      if(lightYieldComputation)
      {
        LYDistro->Write();
        LYDistroCentral->Write();
      }
    }
  }

  fPlots->cd();
  TDirectory *configDir = fPlots->mkdir("Configuration");
  configDir->cd();
  TNamed HostNameD("Hostname",HostNameString.c_str());
  TNamed PWDNameD("PWD",PWDstring.c_str());
  TNamed FilesNameD("InputFiles",InputFiles.str().c_str());
  TNamed ConfigNameD("ConfigFile",streamConfigFile.str().c_str());
  HostNameD.Write();
  PWDNameD.Write();
  FilesNameD.Write();
  ConfigNameD.Write();
  //----------------------------------------------------------//






  if(saveAnalysisTree) // save the TTree created for the analysis, if the user requires it in the config file
  {
    std::cout << "Saving analysis TTree to file " <<  saveAnalysisTreeName.c_str() << std::endl;
    TFile* fFile = new TFile(saveAnalysisTreeName.c_str(),"recreate");
    fFile->cd();
    tree->Write();
    fFile->Close();
  }
  fPlots->Close();
  //----------------------------------------------------------//
  if(usingTaggingBench)
  {
    doiFile.close();
    if(calcDoiResWithDelta)
    {
      AltDoiFile.close();
      altDoiCalcFile.close();
    }
  }

  if(saturationRun)
  {
    saturationFile.close();
  }

  delete crystal;
  delete mppc;
  delete module;

  //     std::cout << "-------------------------------------" << std::endl;

  return 0;
}
