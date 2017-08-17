//------------------------------------------------------------//
//                                                            //
//  PROGRAM FOR matrix calibration USING DT5740               //
//                                                            //
//------------------------------------------------------------//

// compile with
// g++ -o ModuleCalibration ModuleCalibration.cpp `root-config --cflags --glibs` -lSpectrum -lMLP -lTreePlayer

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

#include "ConfigFile.h"
#include "InputFile.h"
#include "Element.h"
#include "Crystal.h"
#include "Module.h"
#include "Mppc.h"
#include "ModuleCalibration.h"

// #include <omp.h>
//test

#define ENERGY_RESOLUTION 0.12
#define ENERGY_RESOLUTION_SATURATION_CORRECTION 0.25



// Double_t thetaFunction(Double_t *x, Double_t *par);

int main (int argc, char** argv)
{

  gStyle->SetOptStat(0);

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
  std::cout<<"#           New Clear PEM module calibration              #"<<std::endl;
  std::cout<<"#                                                         #"<<std::endl;
  std::cout<<"###########################################################"<<std::endl;
  std::cout<<"\n\n"<<std::endl;
  std::cout<<"=====>   C O N F I G U R A T I O N   <====\n"<<std::endl;

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
  ConfigFile config(ConfigFileName); // create a ConfigFile object
  InputFile input(argc,argv,config); // read the input chain of root files, passing the inputs and the config object

  //decide whether to load the ttree from file or produce it, and whether to save it or not to file
  std::string loadAnalysisTreeName         = config.read<std::string>("loadAnalysisTreeName","0");       // look for input analysis ttree in config file. if no input, it will produced by this program
  std::string saveAnalysisTreeName         = config.read<std::string>("saveAnalysisTreeName","0");       // look for filename to save analysis ttree in config file. if no filename, analysis tree won't be saved

  bool loadAnalysisTree;
  bool saveAnalysisTree          = config.read<bool>("saveAnalysisTree",0);            // choice to save or not the analysis TTree, in a file (name chosen above)

  if(loadAnalysisTreeName.compare("0") == 0) // no file provided to load ttree
    loadAnalysisTree = false;
  else  // file provided
    loadAnalysisTree = true;

  //   else // analysis tree is produced by the analysis, not loaded
  //   {
  if(!saveAnalysisTree) // user didn't set to save ttree
  {
    if(saveAnalysisTreeName.compare("0") != 0) // but gave a name to save the ttree
    {
      saveAnalysisTree = true;  // then set the ttree to be saved
    }
  }
  else //user set the flag to save ttree
  {
    if(saveAnalysisTreeName.compare("0") == 0) // but didn't give a name to save the ttree
    {
      saveAnalysisTreeName = "temp.root"; // then set the ttree to be saved into temp.root file
    }
  }
  //   }

  if(loadAnalysisTree) // anyway there's no point saving the analysis tree if it's not produced by this analysis but just loaded from file
  {
    saveAnalysisTree = false;
  }
  //     std::cout << "***** Analysis TTree is input by the user, it won't be overwritten or saved to another file *****" << std::endl;
  //   }
  //   else
  //   {
  //     std::cout << "Analysis TTree will be produced by this program" << std::endl;
  //   }

  //check

  std::cout<<"\n"<<std::endl;
  std::cout<<"###########################################################"<<std::endl;
  std::cout<<"#                                                         #"<<std::endl;
  if(loadAnalysisTree)
    std::cout << "# Analysis TTree loaded from file " << loadAnalysisTreeName.c_str() << std::endl;
  else
    std::cout << "# Analysis TTree produced by this program " << std::endl;
  if(saveAnalysisTree)
    std::cout << "# Analysis TTree will be saved to file " << saveAnalysisTreeName.c_str() << std::endl;
  else
    std::cout << "# No analysis TTree file will be saved " << std::endl;
  std::cout<<"#                                                         #"<<std::endl;
  std::cout<<"###########################################################"<<std::endl;
  std::cout<<"\n"<<std::endl;


  if(!loadAnalysisTree)
  {
    input.ImportTChain(argc,argv);
    input.PrepareTTree();
    input.FillTree();                        // create the TTree that will be used in analysis


  }
  else    // otherwise load it from the indicated file
  {
    TFile *fTemp = new TFile(loadAnalysisTreeName.c_str());
    TTree *TempTree =  (TTree*) fTemp->Get("adc");
    input.SetTree(TempTree);
  }

  int ncrystalsx                = config.read<int>("ncrystalsx",2);                 // number of crystals in x direction per mppc - default to 2 if the key is not found in the config file
  int ncrystalsy                = config.read<int>("ncrystalsy",2);                 // number of crystals in y direction per mppc - default to 2 if the key is not found in the config file
  int nmppcx                    = config.read<int>("nmppcx",2);                     // number of mppc in x direction per mppc - default to 2 if the key is not found in the config file
  int nmppcy                    = config.read<int>("nmppcy",2);                     // number of mppc in y direction per mppc - default to 2 if the key is not found in the config file
  int nmodulex                  = config.read<int>("nmodulex",1);                   // number of modules in x direction per mppc - default to 1 if the key is not found in the config file
  int nmoduley                  = config.read<int>("nmoduley",1);                   // number of modules in y direction per mppc - default to 1 if the key is not found in the config file
  int histo1Dmax                = config.read<int>("histo1Dmax");                   // max of the 1D charge histograms (in ADC channels)
  int histo1Dbins               = config.read<int>("histo1Dbins");                  // number of bins of the 1D charge histograms
  int histoLYmax                = config.read<int>("histoLYmax",20000);                   // max of the 1D charge histograms (in ADC channels)
  int histoLYbins               = config.read<int>("histoLYbins",500);                  // number of bins of the 1D charge histograms
  double qe                     = config.read<double>("qe",0.35);                   // mppc QE
  double gainMPPC               = config.read<double>("gainMPPC",1.25e6);           // mppc gain
  double chargeBinningADC       = config.read<double>("chargeBinningADC",156e-15);  // adc charge binning
  double sourceMeV              = config.read<double>("sourceMeV",0.511);           // gamma peak in MeV
  int histo2DchannelBin         = config.read<int>("histo2DchannelBin");            // number of bins of the 2D flood histograms, for single channels
  int histo2DglobalBins         = config.read<int>("histo2DglobalBins");            // number of bins of the 2D flood histograms, for entire module
  //   int histo3DchannelBin         = config.read<int>("histo3DchannelBin",100);            // number of bins of the 3D flood histograms, for single channels
  int histo3DglobalBins         = config.read<int>("histo3DglobalBins");            // number of bins of the 3D flood histograms, for entire module
  int taggingPeakMin            = config.read<int>("taggingPeakMin",8000);          // min range of tagging crystal photopeak, in ADC channels - to help TSpectrum
  int taggingPeakMax            = config.read<int>("taggingPeakMax",12000);         // max range of tagging crystal photopeak, in ADC channels - to help TSpectrum
  int clusterLevelPrecision     = config.read<int>("clusterLevelPrecision",10);     // precision of the level search when separating the cluster of 3D points

  float taggingPosition         = config.read<float>("taggingPosition");            // position of the tagging bench in mm
  bool usingTaggingBench        = config.read<bool>("usingTaggingBench");           // true if the input is using tagging bench, false if not
  int taggingCrystalChannel     = config.read<int>("taggingCrystalChannel");        // input channel where the tagging crystal information is stored
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
  float userBroadCut            = config.read<float>("userBroadCut",1750.0);            // if in backgroundRun, cut to get rid of low energy events is not done on photopeak search but by user input (default 1750ch)
  float thresholdKev            = config.read<float>("thresholdKev",1.0);
  float wThreshold              = config.read<float>("wThreshold",0.1);                 // Threshold for w plots limits
  double crystalz               = config.read<double>("crystalz",15);
  double qCalVsIJmax            = config.read<double>("qCalVsIJmax",100);               // max of the 2d qCal values plot (starts from 0)
  double RealmCalVsIJmin        = config.read<double>("RealmCalVsIJmin",-400);          // min of the 2d mCal values plot
  double RealmCalVsIJmax        = config.read<double>("RealmCalVsIJmax",400);           // max of the 2d mCal values plot
  double mCalVsIJmax            = config.read<double>("mCalVsIJmax",400);               // max of the 2d abs(mCal) values plot (starts from 0)
  double DoiResolutionVsIJmax   = config.read<double>("DoiResolutionVsIJmax",10);       // max of the 2d DoiResolution values plot (starts from 0) - it's mm
  double EnergyResolutionVsIJmax= config.read<double>("EnergyResolutionVsIJmax",0.3);   // max of the 2d EnergyResolution values plot (starts from 0)
  double LYvsIJmax              = config.read<double>("LYvsIJmax",40000);
  double PeakPositionVsIJmax    = config.read<double>("PeakPositionVsIJmax",12000);     // max of the 2d PeakPosition values plot (starts from 0)  - it's ADC channels
  int wHistogramsBins           = config.read<int>("wHistogramsBins",250);
  double userSigmaW             = config.read<double>("userSigmaW",-1);                 // sigma of w distros for "pointlike" excitation, fixed externally by the user. If nothing specified in the config file, it will be calculated by fitting the rise of w histogram
  int doiColumnOffset           = config.read<int>("doiColumnOffset",0);                // for DOI output, fix the column i by adding this quantity. if not stated, 0 by default
  double energyCorrectionMin    = config.read<double>("energyCorrectionMin",0.25);      // once the wmin and wmax are found for each w histo, choose at which point to start and to stop the linear fitting
  double energyCorrectionMax    = config.read<double>("energyCorrectionMax",0.75);      // (as percentage from min to max)
  double lambda511              = config.read<double>("lambda511",12.195); //everything in mm
  bool wAllChannels             = config.read<bool>("wAllChannels",0);                  // whether we use the sum of all channels to compute w (true = 1) of just the neighbours (false = 0). Deafult to false.
  bool comptonAnalysis          = config.read<bool>("comptonAnalysis",0);               //wheter to perform or not the compton recovery analysis part. Default to false
  bool lightYieldComputation    = config.read<bool>("lightYieldComputation",0);         //wheter to perform or not the light yield calculation. Default to false
  // float taggingPosition        = config.read<float>("taggingPosition");
  // --- paramenters for roto-translations to separate the nXn peaks
  // lateral, not corners
  //   double base_lateralQ1         = config.read<double>("lateralQ1",0.905);           // right and left
  //   double base_lateralQ2         = config.read<double>("lateralQ2",1.1);             // top and bottom
  //   double base_lateralDeltaU     = config.read<double>("lateralDeltaU",1);           // used for right and left
  //   double base_lateralDeltaV     = config.read<double>("lateralDeltaV",1);           // used for top and bottom
  //   double base_lateralRescaleRL  = config.read<double>("lateralRescaleRL",1.5);      // used for right and left
  //   double base_lateralRescaleTB  = config.read<double>("lateralRescaleTB",2);        // used for top and bottom
  // corners
  //   double base_cornerQ1          = config.read<double>("cornerQ1",0.675);            // rotation around Z
  //   double base_cornerQ2          = config.read<double>("cornerQ2",1.41);             // rotation around X
  //   double base_cornerDeltaU      = config.read<double>("cornerDeltaU",3);            // translations
  //   double base_cornerDeltaV      = config.read<double>("cornerDeltaV",2.1);          // translations
  //   double base_cornerRescale     = config.read<double>("cornerRescale",4);           // rescale factor
  bool   onlyuserinput          = config.read<double>("onlyuserinput",0);           // ignore 2d automatic fitting
  // set output file name
  std::string outputFileName = config.read<std::string>("output");
  outputFileName += ".root";
  // create "sum channels" string, a string to have the variable "sum of all channels" to be used later
  // first get the input digitizer channels
  std::string digitizer_s   = config.read<std::string>("digitizer");
  std::vector<std::string> digitizer_f;
  config.split( digitizer_f, digitizer_s, "," );
  std::vector<int> digitizer;
  for(int i = 0 ; i < digitizer_f.size() ; i++)
  {
    config.trim(digitizer_f[i]);
    digitizer.push_back(atoi(digitizer_f[i].c_str()));
  }
  // create the string
  std::stringstream sSumChannels;
  std::string SumChannels;
  sSumChannels << "ch" <<  digitizer[0];
  for(int i = 1 ; i < digitizer.size() ; i++)
    sSumChannels << "+ch" << digitizer[i];
  SumChannels = sSumChannels.str();
  //   std::cout << SumChannels << std::endl;
  //----------------------------------------------------------//


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
  //  Load electronic calibration values                      //
  //----------------------------------------------------------//
  std::string calibFileName = config.read<std::string>("doiComparison","0");
  std::vector<Double_t> m_tag;
  //   std::vector<Double_t> m_cal;
  std::vector<Double_t> q_tag;
  //   std::vector<Double_t> q_cal;

  if(calibFileName != "0")
  {
    std::ifstream inFile_calib;
    inFile_calib.open(calibFileName.c_str(),std::ios::in);
    while(!inFile_calib.eof())
    {
      Double_t foo, a, b, c, d;
      inFile_calib >> foo >> a >> b;
      //       m_cal.push_back(a);
      //       q_cal.push_back(b);
      m_tag.push_back(c);
      q_tag.push_back(d);
      //       inFile_calib >> foo >>m_cal[i]>>q_cal[i]>>m_tag[i]>>q_tag[i];
      //       i++;
    }
    inFile_calib.close();
  }


  //----------------------------------------------------------//
  //  Plots and spectra                                       //
  //----------------------------------------------------------//
  // get the TTree, to plot the spectra


  TTree* tree = input.GetTree();

  // doi bench specific part
  std::ofstream doiFile;
  std::ofstream AltDoiFile;
  std::ifstream altDoiCalcFile;
  TCut triggerPhotopeakCut = "" ;
  TH1F* TaggingCrystalSpectrum;
  TH1F *TriggerSpectrumHighlight;

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

  // MAIN LOOP
  // Loop on modules, mppcs and crystal
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      //usegul string etc
      std::stringstream CutXYZstream;
      CutXYZstream << "FloodX > " << -moduleLateralSideX << " && FloodX < " << moduleLateralSideX  << "&& FloodY > " << -moduleLateralSideY <<   " && FloodY < " << moduleLateralSideY << "&& FloodZ > 0 && FloodZ < 1";
      TCut CutXYZ = CutXYZstream.str().c_str();
      std::cout << "Generating global spectra..." << std::endl;
      TString nameModule;
      std::stringstream varModule;

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

      //3D plot
      // <<<<<<< HEAD
      //       nameModule = "Flood Histogram 3D - Module " + module[iModule][jModule]->GetName();
      //       varModule << "FloodZ:FloodY:FloodX >> " << nameModule;
      //       TH3I* spectrum3dModule = new TH3I(nameModule,nameModule,histo3DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo3DglobalBins,-moduleLateralSideY,moduleLateralSideY,histo3DglobalBins,0,1);
      //       tree->Draw(varModule.str().c_str(),CutXYZ);
      //       spectrum3dModule->GetXaxis()->SetTitle("U");
      //       spectrum3dModule->GetYaxis()->SetTitle("V");
      //       spectrum3dModule->GetZaxis()->SetTitle("W");
      //       module[iModule][jModule]->SetFloodMap3D(spectrum3dModule);
      //       varModule.str("");
      // =======
      //      nameModule = "Flood Histogram 3D - Module " + module[iModule][jModule]->GetName();
      //      varModule << "FloodZ:FloodY:FloodX >> " << nameModule;
      //       std::cout << nameModule << " ... ";
      //      TH3I* spectrum3dModule = new TH3I(nameModule,nameModule,histo3DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo3DglobalBins,-moduleLateralSideY,moduleLateralSideY,histo3DglobalBins,0,1);
      //      tree->Draw(varModule.str().c_str(),CutXYZ);
      //      spectrum3dModule->GetXaxis()->SetTitle("U");
      //      spectrum3dModule->GetYaxis()->SetTitle("V");
      //      spectrum3dModule->GetZaxis()->SetTitle("W");
      //      module[iModule][jModule]->SetFloodMap3D(spectrum3dModule);
      //      varModule.str("");
      //       std::cout << " done" << std::endl;
      //       delete spectrum3dModule;
      // >>>>>>> doiTag

      if(usingTaggingBench)//trigger spectrum
      {
        TaggingCrystalSpectrum =  new TH1F("TaggingCrystalSpectrum","TaggingCrystalSpectrum",1200,0,12000);
        varModule << "Tagging >> TaggingCrystalSpectrum";
        // 	std::cout << nameModule << " ... ";
        TaggingCrystalSpectrum->SetName("TaggingCrystalSpectrum");
        TaggingCrystalSpectrum->GetXaxis()->SetTitle("ADC");
        TaggingCrystalSpectrum->GetYaxis()->SetTitle("Counts");
        TaggingCrystalSpectrum->SetTitle("Spectrum of Tagging Crystal");
        tree->Draw(varModule.str().c_str(),"");
        //restrict the region where to look for peaks. Fix for tspectrum...
        TaggingCrystalSpectrum->GetXaxis()->SetRangeUser(taggingPeakMin,taggingPeakMax);
        //find peak in the tagging crystal
        TSpectrum *sTagCrystal;
        sTagCrystal = new TSpectrum(1);
        Int_t TagCrystalPeaksN = sTagCrystal->Search(TaggingCrystalSpectrum,1,"",0.5);
        Float_t *TagCrystalPeaks = sTagCrystal->GetPositionX();
        TF1 *gaussTag = new TF1("gaussTag", "gaus");
        TaggingCrystalSpectrum->Fit("gaussTag","NQ","",TagCrystalPeaks[0] - 0.075*TagCrystalPeaks[0],TagCrystalPeaks[0] + 0.075*TagCrystalPeaks[0]);
        TaggingCrystalSpectrum->GetXaxis()->SetRangeUser(0,12000);
        //define a TCut for this peak
        double tagPhotopeakMin = gaussTag->GetParameter(1) - 1.5*gaussTag->GetParameter(2);
        double tagPhotopeakMax = gaussTag->GetParameter(1) + 2.0*gaussTag->GetParameter(2);
        std::stringstream tagString;
        tagString << "Tagging > " << tagPhotopeakMin << "&& Tagging < " << tagPhotopeakMax;
        triggerPhotopeakCut = tagString.str().c_str();
        //highlighted spectrum
        TriggerSpectrumHighlight = new TH1F("TriggerSpectrumHighlight","",1200,0,12000);
        varModule.str("");
        varModule << "Tagging >> TriggerSpectrumHighlight";
        TriggerSpectrumHighlight->SetLineColor(3);
        TriggerSpectrumHighlight->SetFillColor(3);
        TriggerSpectrumHighlight->SetFillStyle(3001);
        tree->Draw(varModule.str().c_str(),triggerPhotopeakCut);
        varModule.str("");
        // 	delete sTagCrystal;
        delete gaussTag;
        // 	std::cout << " done" << std::endl;

      }

      nameModule = "Flood Histogram 3D - Module " + module[iModule][jModule]->GetName();
      varModule << "FloodZ:FloodY:FloodX >> " << nameModule;
      //       std::cout << nameModule << " ... ";
      TH3I* spectrum3dModule = new TH3I(nameModule,nameModule,histo3DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo3DglobalBins,-moduleLateralSideY,moduleLateralSideY,histo3DglobalBins,0,1);
      if(usingTaggingBench)
      {
        tree->Draw(varModule.str().c_str(),CutXYZ+triggerPhotopeakCut);
      }
      else{
        tree->Draw(varModule.str().c_str(),CutXYZ);
      }
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
      //and i can't parallelize it because ROOT coders live in the 80s...
      //#pragma omp parallel for
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
        for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
        {
          //these loop are supposed to run in parallel, so some variables need to exist only inside the loop
          //in fact, i can't make them run in parallel because root is not thread safe

          //go on only if the MPPC is on for modular analysis
          if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForModular())
          {
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
            tree->Draw(var.str().c_str(),"");
            spectrumRaw->GetXaxis()->SetTitle("ADC Channels");
            spectrumRaw->GetYaxis()->SetTitle("N");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetRawSpectrum(spectrumRaw);
            var.str("");

            //trigger selected spectrum
            name = "Trigger Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            var << "ch" << channel << " >> " << name;
            TH1F* spectrumTrigger = new TH1F(name,name,histo1Dbins,1,histo1Dmax);
            tree->Draw(var.str().c_str(),CutTrigger);
            spectrumTrigger->GetXaxis()->SetTitle("ADC Channels");
            spectrumTrigger->GetYaxis()->SetTitle("N");

            //set a very broad cut on the trigger spectrum (single channel, not sum) to get rid of low energy events
            std::stringstream broadCutstream;
            if(!backgroundRun)
            {
              TSpectrum *sTrigger;
              sTrigger = new TSpectrum(20);
              Int_t TriggerCrystalPeaksN    = sTrigger->Search(spectrumTrigger,2,"",0.2);
              Float_t *TriggerCrystalPeaks  = sTrigger->GetPositionX();
              Float_t *TriggerCrystalPeaksY = sTrigger->GetPositionY();
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

              broadCutstream << "ch" << channel << ">" << (TriggerCrystalPeaks[TriggerpeakID] / (511.0/thresholdKev) );
            }
            else //otherwise, for backgroundRun set the broadcut to the value chosen by the user (or default to 1750 adc channels)
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
            tree->Draw(var.str().c_str(),CutTrigger);
            spectrumTriggerHighlighted->GetXaxis()->SetTitle("ADC Channels");
            spectrumTriggerHighlighted->GetYaxis()->SetTitle("N");
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTriggerSpectrumHighlighted(spectrumTriggerHighlighted);
            var.str("");

            //standard 2d plot
            name = "Flood Histogram 2D - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            TH2F* spectrum2dMPPC = new TH2F(name,name,histo2DchannelBin,-moduleLateralSideX,moduleLateralSideX,histo2DchannelBin,-moduleLateralSideY,moduleLateralSideY);
            var << "FloodY:FloodX >> " << name;
            tree->Draw(var.str().c_str(),CutXYZ+CutTrigger,"COLZ");
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
            TH3I* spectrum3dMPPC = new TH3I(name,name,histo3DchannelBin,minX3Dplot,maxX3Dplot,histo3DchannelBin,minY3Dplot,maxY3Dplot,histo3DchannelBin,0,1);//FIXME temp
            if(usingTaggingBench){
              tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+triggerPhotopeakCut);
            }
            else{
              tree->Draw(var.str().c_str(),CutXYZ+CutTrigger);
            }
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
            if(usingTaggingBench){
              right_ncrystalsx =1;
            }
            else{
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
            // 	  if(!onlyuserinput)
            // 	  {
            if(usingTaggingBench) //if it's a tagging bench run, check first if the mppc is on for DOI bench measurement
            {
              if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForDoi())
              {
                found = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->FindCrystalCuts(cutg/*,histo3DchannelBin,clusterLevelPrecision*/,1,ncrystalsy);
                // 		mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->Find2Dpeaks(ncrystalsx*ncrystalsy,mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D());
              }
            }
            else
            {
              found = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->FindCrystalCuts(cutg/*,histo3DchannelBin,clusterLevelPrecision*/,ncrystalsx,ncrystalsy);
              // 	      mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->Find2Dpeaks(ncrystalsx*ncrystalsy,mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D());
            }
            // 	  }

            // run on all the possible crystals (i.e. all the crystals coupled to this mppc)
            // <<<<<<< HEAD
            // 	  for(int iCry = 0; iCry < ncrystalsx ; iCry++)
            // =======

            for(int iCry = 0; iCry < right_ncrystalsx ; iCry++)
              // >>>>>>> doiTag
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
                    std::cout << "Generating spectra for crystal " << CurrentCrystal->GetID() << " ..." << std::endl;

                    //-------------------------------------------------------------------------
                    //standard sum spectrum with cut on crystal events, xyz and trigger channel
                    //------------------------------------------------------------------------
                    //draw charge spectrum
                    sname << "Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                    var << SumChannels << " >> " << sname.str();
                    TH1F* spectrumCharge = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                    tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName());
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
                      tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName());
                      spectrumLY->GetXaxis()->SetTitle("[Ph/MeV]");
                      spectrumLY->GetYaxis()->SetTitle("N");
                      CurrentCrystal->SetLYSpectrum(spectrumLY);
                      sname.str("");
                      var.str("");

                      TSpectrum *s;
                      s = new TSpectrum(20);
                      // 		Input[i].SumSpectraCanvas->cd(j+1);
                      Int_t CrystalPeaksN = s->Search(spectrumLY,2,"goff",0.5);
                      Float_t *CrystalPeaks = s->GetPositionX();
                      Float_t *CrystalPeaksY = s->GetPositionY();
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
                      Int_t CrystalPeaksN = s->Search(spectrumCharge,2,"goff",0.5);
                      Float_t *CrystalPeaks = s->GetPositionX();
                      Float_t *CrystalPeaksY = s->GetPositionY();
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
                      sname << "gaussCharge - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
                      TF1 *gauss = new TF1(sname.str().c_str(),  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitmin,fitmax);
                      gauss->SetParameter(0,par0);
                      gauss->SetParameter(1,par1);
                      gauss->SetParameter(2,par2);
                      spectrumCharge->Fit(sname.str().c_str(),"Q","",fitmin,fitmax);
                      //store the mean and sigma in the crystal
                      if(gauss->GetParameter(1) > 0) // otherwise the fit was very wrong..)
                        CurrentCrystal->SetPhotopeak(gauss->GetParameter(1),std::abs(gauss->GetParameter(2)));
                      CurrentCrystal->SetFit(gauss);
                      // 		std::cout << "Photopeak Mean for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakPosition() << std::endl;
                      // 		std::cout << "Photopeak Sigma for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakSigma() << std::endl;
                      // 		std::cout << "Photopeak Energy Resolution FWHM for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakEnergyResolution() << std::endl;
                      //Compute the energy Tcut
                      std::stringstream streamEnergyCut;
                      EnergyCutMin = gauss->GetParameter(1) - 2.0*std::abs(gauss->GetParameter(2));
                      EnergyCutMax = gauss->GetParameter(1) + 4.0*std::abs(gauss->GetParameter(2));
                      streamEnergyCut << SumChannels << " > " << EnergyCutMin << " && " << SumChannels << " < " << EnergyCutMax;
                      // 		  TCut PhotopeakEnergyCutCorrected;
                      PhotopeakEnergyCut  = streamEnergyCut.str().c_str();
                      sname.str("");

                      // then prepare the highlighted spectrum and store it in the crystal
                      sname << "Hg Charge Spectrum - Crystal " << CurrentCrystal->GetID();
                      var << SumChannels << " >> " << sname.str();
                      TH1F* spectrumChargeHighlighted = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                      tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCut);
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
                    TH1F* spectrumHistoW = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,0,1);
                    tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCut+triggerPhotopeakCut);
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

                    if(usingTaggingBench) //done on the standard W because there's no ADCvsW correction in tagging bench
                    {
                      sname << "gaussW histogram - Crystal " << CurrentCrystal->GetID();
                      TF1 *gaussW = new TF1(sname.str().c_str(),  "gaus",0,1);
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
                        for(int kAltDoi = 0; kAltDoi < inputDoi.size(); kAltDoi++)
                        {
                          int ik = inputDoi[kAltDoi].i;
                          int jk = inputDoi[kAltDoi].j;
                          if( ((CurrentCrystal->GetI() + doiColumnOffset ) == ik) && CurrentCrystal->GetJ() == jk)
                          {
                            sname << "Alternate DOI res - Crystal " << CurrentCrystal->GetID();
                            var << "(FloodZ * " << inputDoi[kAltDoi].m << "+" << inputDoi[kAltDoi].q << ") - (TreeZPosition  ) >> " << sname.str();
                            TH1F* spectrumAltDoiRes = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,-10,10); //FIXME boundaries hardcoded
                            tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCut+triggerPhotopeakCut);
                            spectrumAltDoiRes->GetXaxis()->SetTitle("Doi_calculated - DoiTagPosition [mm]");
                            spectrumAltDoiRes->GetYaxis()->SetTitle("N");
                            var.str("");
                            sname.str("");

                            sname << "gaussFit Alternate - Crystal " << CurrentCrystal->GetID();
                            TF1 *gaussFitAlt = new TF1(sname.str().c_str(),  "gaus",-10,10); //FIXME boundaries hardcoded
                            Int_t fitStatusAlt = spectrumAltDoiRes->Fit(sname.str().c_str(),"QR");
                            if(!fitStatus)
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
                      //it will be useful fot doi correction
                      //long long int nPoints;
                      sname << "ADC channels vs. W - Crystal " << CurrentCrystal->GetID();
                      var << SumChannels << ":FloodZ >> " << sname.str();
                      TH2F* spectrum2dADCversusW = new TH2F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,0,1,histo1Dbins,0,histo1Dmax);
                      tree->Draw(var.str().c_str(),CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCut+w20percCut,"COLZ");
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
                        tree->Draw(var.str().c_str(),CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName(),"COLZ");
                        spectrumChargeCorrected->GetXaxis()->SetTitle("ADC channels");
                        spectrumChargeCorrected->GetYaxis()->SetTitle("Counts");
                        CurrentCrystal->SetCorrectedSpectrum(spectrumChargeCorrected);
                        var.str("");
                        sname.str("");
                        //find peaks in each crystal spectrum, with TSpectrum
                        TSpectrum *s_corr;
                        s_corr = new TSpectrum(20);
                        // 		Input[i].SumSpectraCanvas->cd(j+1);
                        Int_t CrystalPeaksN_corr = s_corr->Search(spectrumChargeCorrected,2,"goff",0.5);
                        Float_t *CrystalPeaks_corr = s_corr->GetPositionX();
                        Float_t *CrystalPeaksY_corr = s_corr->GetPositionY();
                        // 		  delete s_corr;
                        float maxPeak_corr = 0.0;
                        int peakID_corr = 0;
                        for (int peakCounter = 0 ; peakCounter < CrystalPeaksN_corr ; peakCounter++ )
                        {
                          if(CrystalPeaks_corr[peakCounter] > maxPeak_corr)
                          {
                            maxPeak_corr = CrystalPeaks_corr[peakCounter];
                            peakID_corr = peakCounter;
                          }
                        }
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
                        streamEnergyCutCorrected << baseVar.str() << " > " << gauss_corr->GetParameter(1) - 2.5*std::abs(gauss_corr->GetParameter(2)) << " && " << baseVar.str() << " < " << gauss_corr->GetParameter(1) + 3.0*std::abs(gauss_corr->GetParameter(2));
                        PhotopeakEnergyCutCorrected = streamEnergyCutCorrected.str().c_str();
                        // 		CurrentCrystal->SetSpectrum(*spectrum);
                        sname.str("");

                        //then prepare the highlighted spectrum and store it in the crystal
                        sname << "Hg Charge Spectrum Correctd - Crystal " << CurrentCrystal->GetID();
                        // 			var << "("  <<  SumChannels<< " ) - ( ( FloodZ - " <<  meanW20 << " ) * ( " << parM << ") ) >> " << sname.str();
                        var << baseVar.str() << " >> " << sname.str();
                        TH1F* spectrumChargeCorrectedHighlighted = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);
                        tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCutCorrected);
                        spectrumChargeCorrectedHighlighted->GetXaxis()->SetTitle("ADC Channels");
                        spectrumChargeCorrectedHighlighted->GetYaxis()->SetTitle("N");
                        CurrentCrystal->SetHighlightedSpectrumCorrected(spectrumChargeCorrectedHighlighted);
                        var.str("");
                        sname.str("");

                      }
                    }
                    // a 3d historgram for this crystal, mainly to check the 3d cut
                    sname << "Flood Histogram 3D - Crystal " << CurrentCrystal->GetID();
                    var << "FloodZ:FloodY:FloodX >> " << sname.str();
                    // NB histo3DchannelBin is taken from the MPPC element earlier in the code
                    TH3I* spectrum3dCrystal = new TH3I(sname.str().c_str(),sname.str().c_str(),histo3DchannelBin,minX3Dplot,maxX3Dplot,histo3DchannelBin,minY3Dplot,maxY3Dplot,histo3DchannelBin,0,1);
                    tree->Draw(var.str().c_str(),CutXYZ + CutTrigger + CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName());
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
                    tree->Draw(var.str().c_str(),CutXYZ + CutTrigger + CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName());
                    spectrum2dCrystal->GetXaxis()->SetTitle("U");
                    spectrum2dCrystal->GetYaxis()->SetTitle("V");
                    CurrentCrystal->SetFloodMap2D(spectrum2dCrystal);
                    sname.str("");
                    var.str("");

                    //density histogram - 1d histo of entries per bin in the cutg volume
                    // 		    sname << "Density Histogram - Crystal " << CurrentCrystal->GetID();
                    // 		    Int_t u,v,w;
                    // 		    spectrum3dCrystal->GetMaximumBin(u,v,w); // get the maximum bin of the 3d histo
                    // 		    //get ax bin content
                    // 		    TH1F* densityHisto = new TH1F(sname.str().c_str(),sname.str().c_str(),spectrum3dCrystal->GetBinContent(u,v,w)-1,1,spectrum3dCrystal->GetBinContent(u,v,w));
                    // 		    densityHisto->GetXaxis()->SetTitle("Numb of Entries");
                    // 		    int NbinX = spectrum3dCrystal->GetXaxis()->GetNbins();
                    // 		    int NbinY = spectrum3dCrystal->GetYaxis()->GetNbins();
                    // 		    int NbinZ = spectrum3dCrystal->GetZaxis()->GetNbins();
                    // 		    for(int iContent = 1 ; iContent < NbinX+1 ; iContent++)
                    // 		    {
                    // 		      for(int jContent = 1 ; jContent < NbinY+1 ; jContent++)
                    // 		      {
                    // 			for(int kContent = 1 ; kContent < NbinZ+1 ; kContent++)
                    // 			{
                    // 			  densityHisto->Fill(spectrum3dCrystal->GetBinContent(iContent,jContent,kContent));
                    // 			}
                    // 		      }
                    // 		    }
                    // 		    CurrentCrystal->SetDensityHisto(densityHisto);
                    // 		    sname.str("");

                    // Histogram 2d of time evolution
                    sname << "ADC channels vs. Time - Crystal " << CurrentCrystal->GetID();
                    var << SumChannels << ":ExtendedTimeTag >> " << sname.str();
                    TH2F* spectrum2dVersusTime = new TH2F(sname.str().c_str(),sname.str().c_str(),250,0,tree->GetMaximum("ExtendedTimeTag"),histo1Dbins,0,histo1Dmax);
                    tree->Draw(var.str().c_str(),CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName(),"COLZ");
                    spectrum2dVersusTime->GetXaxis()->SetTitle("ExtendedTimeTag");
                    spectrum2dVersusTime->GetYaxis()->SetTitle("ADC channels");
                    CurrentCrystal->SetVersusTime(spectrum2dVersusTime);
                    var.str("");
                    sname.str("");

                    //histogram of w versus adc channels - this time without the photopeak cut (so it looks nicer in the paper...)
                    sname << "Complete ADC channels vs. W - Crystal " << CurrentCrystal->GetID();
                    var << SumChannels << ":FloodZ >> " << sname.str() ;
                    TH2F* spectrum2dADCversusWComplete = new TH2F(sname.str().c_str(),sname.str().c_str(),100,0,1,histo1Dbins,0,histo1Dmax);
                    tree->Draw(var.str().c_str(),CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName(),"COLZ");
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
                    TH1F* spectrumHistoWCorrected = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,0,1);
                    tree->Draw(var.str().c_str(),CutXYZ + CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCutCorrected+triggerPhotopeakCut);
                    spectrumHistoWCorrected->GetXaxis()->SetTitle("W");
                    spectrumHistoWCorrected->GetYaxis()->SetTitle("N");

                    //clone this histogram, to smooth and find consistent relevant points. fitting will be performed on the non-smoothed version
                    TH1F* spectrumHistoWCorrectedClone = (TH1F*) spectrumHistoWCorrected->Clone();
                    spectrumHistoWCorrectedClone->SetName("Smooth W corrected");
                    spectrumHistoWCorrectedClone->Smooth(10); //smooth the w histo
                    //relevant points in the w corrected plot, to be used later

                    Int_t  FirstBinAbove20perc = spectrumHistoWCorrectedClone->FindFirstBinAbove(wThreshold*spectrumHistoWCorrectedClone->GetMaximum());
                    Int_t  LastBinAbove20perc  = spectrumHistoWCorrectedClone->FindLastBinAbove(wThreshold*spectrumHistoWCorrectedClone->GetMaximum());
                    double FirstWAbove20perc   = spectrumHistoWCorrectedClone->GetBinCenter(FirstBinAbove20perc);
                    double LastWAbove20perc    = spectrumHistoWCorrectedClone->GetBinCenter(LastBinAbove20perc);

                    double AverageW            = (FirstWAbove20perc + LastWAbove20perc) /2.0;                                                                                  // average w between first 20% and last 20%
                    Int_t  AverageBin          = (int) (FirstBinAbove20perc + LastBinAbove20perc) /2.0;                                                                        // average bin between first 20% and last 20%
                    // look for the first "peak" in w
                    spectrumHistoWCorrectedClone->GetXaxis()->SetRange(FirstBinAbove20perc,AverageBin);                               // restric the range where to look for the max to first above 20% and average
                    double FirstWpeak        = spectrumHistoWCorrectedClone->GetBinCenter(spectrumHistoWCorrectedClone->GetMaximumBin());  // get the w where the first max is
                    double FirstWpeakValue   = spectrumHistoWCorrectedClone->GetBinContent(spectrumHistoWCorrectedClone->GetMaximumBin());                                         // get value of max in this range

                    //FIXME quick fix for second peak, i.e. take the mx in a much shorter w range
                    spectrumHistoWCorrectedClone->GetXaxis()->SetRange(AverageBin +  ((int) (0.7*(LastBinAbove20perc - AverageBin))),LastBinAbove20perc);
                    double LastWpeak        = spectrumHistoWCorrectedClone->GetBinCenter(spectrumHistoWCorrectedClone->GetMaximumBin());  // get the w where the last max is
                    double LastWpeakValue   = spectrumHistoWCorrectedClone->GetBinContent(spectrumHistoWCorrectedClone->GetMaximumBin());
                    // 		spectrumHistoWCorrected->GetXaxis()->SetRange(1,250); //reset the w plots limits
                    // fit w with theta function
                    // 		    TF1 *w_fit_func = new TF1("fa1",thetaFunction,0,1,3);
                    // 		    w_fit_func->SetParameter( 0, FirstWAbove20perc); // on this w histo, the first bin above 20% max
                    // 		    w_fit_func->SetParameter( 1, LastWAbove20perc);  // on this w histo, the last bin above 20% max
                    // 		    w_fit_func->SetParameter( 2, FirstWpeakValue);
                    // 		    spectrumHistoWCorrected->Fit(w_fit_func,"QNR");
                    //save the spectrum, fit and values in the crystal
                    CurrentCrystal->SetHistoWCorrected(spectrumHistoWCorrected);
                    // 		    CurrentCrystal->SetWbegin(w_fit_func->GetParameter(0));
                    // 		    CurrentCrystal->SetWend(w_fit_func->GetParameter(1));


                    // 		    CurrentCrystal->SetThetaFit(w_fit_func);
                    //fit the left part of the w plot with a gaussian, to get delta w //
                    // 		    TF1 *gaussDeltaW = new TF1("gaussDeltaW","[0]*exp(-0.5*((x-[1])/[2])**2)",0,FirstWpeak); //fitting function defined only in the fitting range (otherwise somehow i cannot make it work)
                    // 		    gaussDeltaW->SetParameter( 0, FirstWpeakValue);  // starting point as the maximum value
                    // 		    gaussDeltaW->SetParameter( 1, FirstWpeak);  // fix center to the peak value
                    // // 		    gaussDeltaW->FixParameter( 1, FirstWpeak);  // fix center to the peak value
                    // 		    gaussDeltaW->SetParameter( 2, 0.015);
                    // 		    gaussDeltaW->SetLineColor(3);
                    // 		    spectrumHistoWCorrected->Fit(gaussDeltaW,"QNR");

                    // 		    TF1 *gaussDeltaW_2 = new TF1("gaussDeltaW_2","[0]*exp(-0.5*((x-[1])/[2])**2)",LastWpeak,1); //fitting function defined only in the fitting range (otherwise somehow i cannot make it work)
                    // 		    gaussDeltaW_2->SetParameter( 0, LastWpeakValue);  // starting point as the maximum value
                    // 		    gaussDeltaW_2->SetParameter( 1, LastWpeak);  // fix center to the peak value
                    // 		    gaussDeltaW_2->SetParameter( 2, 0.015);
                    // 		    gaussDeltaW_2->SetLineColor(4);
                    // 		    spectrumHistoWCorrected->Fit(gaussDeltaW_2,"QNR");



                    // 		    CurrentCrystal->SetDeltaWfit(gaussDeltaW);
                    // 		    CurrentCrystal->SetDeltaWfit_2(gaussDeltaW_2);

                    //FIXME which one do we choose?
                    //CurrentCrystal->SetDeltaW( (gaussDeltaW->GetParameter(2) + gaussDeltaW_2->GetParameter(2))/2.0 ); // average
                    // 		    CurrentCrystal->SetDeltaW( gaussDeltaW->GetParameter(2) ); // just left one

                    // min among the two
                    // 		    double minDeltaW = gaussDeltaW->GetParameter(2);
                    // 		    if(gaussDeltaW_2->GetParameter(2) < minDeltaW) minDeltaW = gaussDeltaW_2->GetParameter(2);
                    // 		    CurrentCrystal->SetDeltaW( minDeltaW ); // min among the two found

                    CurrentCrystal->SetHistoWCorrectedSmooth(spectrumHistoWCorrectedClone);

                    //DEBUG
                    // 		    std::ofstream wfile;
                    // 		    wfile.open("wfile.dat", std::ofstream::out | std::ofstream::app);
                    // 		    wfile <<  wThreshold << "\t" << FirstWAbove20perc << "\t" << LastWAbove20perc <<"\t"<< gaussDeltaW->GetParameter(1) - gaussDeltaW->GetParameter(2) * TMath::Sqrt(2.0*TMath::Log(2)) << "\t"<< gaussDeltaW_2->GetParameter(1) + gaussDeltaW_2->GetParameter(2) * TMath::Sqrt(2.0*TMath::Log(2)) << "\t" << w_fit_func->GetParameter(0) << "\t" << w_fit_func->GetParameter(1) << "\t" << gaussDeltaW->GetParameter(0) << "\t" << gaussDeltaW->GetParameter(1) << "\t" << gaussDeltaW->GetParameter(2) << std::endl;
                    // 		    wfile.close();

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
                    TH1F* cumulativeW = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,0,1);
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
                        calibrationZ.push_back(-(sumPdf*crystalz) + crystalz); // if it's a background run, interaction probability is constant everywhere, so z is just a rescale of the cumulative
                        else
                          // if it's a far source run, interaction probability is exponential, and relation between z and w is give by
                          // z = L + l * ln( 1 - (1 - exp(-L/l))*integral_0^w(PDF(w)dw) )
                          // where
                          // L = crystal length
                          // l = interaction length of 511 gammas in lyso (from literature, 12.195 mm)
                          // integral_0^w(PDF(w)dw) = the cumulative of PDF(w) from 0 to w
                          calibrationZ.push_back( crystalz +  lambda511 * TMath::Log( 1.0 - (1.0 - TMath::Exp(-(crystalz/lambda511)) )* sumPdf )  );
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

                    //w(z) for compton analysis
                    TGraph *wzgraph = new TGraph(calibrationW.size(),&calibrationZ[0],&calibrationW[0]);
                    sname.str("");
                    sname << "w(z) Plot - Crystal " << CurrentCrystal->GetID();
                    wzgraph->SetTitle(sname.str().c_str());
                    wzgraph->SetName(sname.str().c_str());
                    wzgraph->GetXaxis()->SetTitle("Z [mm]");
                    wzgraph->GetYaxis()->SetTitle("W");
                    sname.str("");



                    // 		    wfile.close();

                    CurrentCrystal->SetWZgraph(wzgraph);
                    CurrentCrystal->SetCumulativeW(cumulativeW);

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

                      // TGraphDelaunay*** interpolationGraph;
                      // interpolationGraph = new TGraphDelaunay** [nmppcx];
                      // for(int iCal = 0; iCal < nmppcx ; iCal++) interpolationGraph[iCal] = new TGraphDelaunay*[nmppcx];

                      // TH3I***  ComptonCalibationHistogram;
                      // ComptonCalibationHistogram = new TH3I** [4];
                      // for(int iCal = 0; iCal < 4 ; iCal++) ComptonCalibationHistogram[iCal] = new TH3I*[4];



                      //DEBUG
                      // just one
                      // TH3I* oneHisto = new TH3I(
                      //   "one",
                      //   "one",
                      //   wHistogramsBins,
                      //   0,
                      //   1,
                      //   histo1Dbins,
                      //   0,
                      //   histo1Dmax,
                      //   histo1Dbins,
                      //   0,
                      //   histo1Dmax
                      //   );
                      //
                      // var << "(ch" << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel() << "):(" << SumChannels << "):(FloodZ)" << " >> one" ;
                      // tree->Draw(var.str().c_str(),CutXYZ + CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName());
                      // var.str("");
                      // CurrentCrystal->SetOne(oneHisto);


                      // fill all the TGraphDelaunay
                      // for(int iNeighbour = -1 ; iNeighbour < 2 ; iNeighbour++)
                      // {
                      //   for(int jNeighbour = -1 ; jNeighbour < 2 ; jNeighbour++)
                      //   {
                      //     if( (iMppc + iNeighbour >= 0) && (iMppc + iNeighbour <= nmppcx) )
                      //     {
                      //       if( (jMppc + jNeighbour >= 0) && (jMppc + jNeighbour <= nmppcy) )
                      //       {
                      for(int iComptMppc = 0; iComptMppc < nmppcx ; iComptMppc++)
                      {
                        for(int jComptMppc = 0; jComptMppc < nmppcy ; jComptMppc++)
                        {
                              sname.str("");
                              var.str("");
                              //FIXME temporary BAD BAD BAD!!!
                              //TODO find a way to fix the max values
                              // double maxZhisto = histo1Dmax;
                              // if(iNeighbour == -1 | iNeighbour == 1){
                              //   if(jNeighbour == -1 | jNeighbour == 1){
                              //      maxZhisto = histo1Dmax / 20.0;
                              //   }
                              //   else{
                              //     maxZhisto = histo1Dmax / 6.0;
                              //   }
                              // }
                              // else{
                              //   if(jNeighbour == -1 | jNeighbour == 1) maxZhisto = histo1Dmax / 6.0;
                              // }
                              // if(iNeighbour == 0) && (jNeighbour ==0) maxZhisto = histo1Dmax;
                              // if(iNeighbour == -1) && (jNeighbour ==0) maxZhisto = histo1Dmax;
                              // if(CurrentCrystal->GetID() == 28)
                              // {
                              // if( (iNeighbour == 0) && (jNeighbour ==0) )
                              //   maxZhisto = histo1Dmax;
                              // if( (iNeighbour == 0) && (jNeighbour == -1) )
                              //   maxZhisto = 800.0;
                              // if( (iNeighbour == 0) && (jNeighbour == 1) )
                              //   maxZhisto = 200.0;
                              // if( (iNeighbour == -1) && (jNeighbour ==0) )
                              //   maxZhisto = 150.0;
                              // if( (iNeighbour == -1) && (jNeighbour == -1) )
                              //   maxZhisto = 100.0;
                              // if( (iNeighbour == -1) && (jNeighbour == 1) )
                              //   maxZhisto = 100.0;
                              // if( (iNeighbour == 1) && (jNeighbour ==0) )
                              //   maxZhisto = 600.0;
                              // if( (iNeighbour == 1) && (jNeighbour == -1) )
                              //   maxZhisto = 300.0;
                              // if( (iNeighbour == 1) && (jNeighbour == 1) )
                              //   maxZhisto = 100.0;
                              // }
                              // if(CurrentCrystal->GetID() == 29)
                              // {
                              // if( (iNeighbour == 0) && (jNeighbour ==0) )
                              //   maxZhisto = histo1Dmax;
                              // if( (iNeighbour == 0) && (jNeighbour == -1) )
                              //   maxZhisto = 200.0;
                              // if( (iNeighbour == 0) && (jNeighbour == 1) )
                              //   maxZhisto = 800.0;
                              // if( (iNeighbour == -1) && (jNeighbour ==0) )
                              //   maxZhisto = 150.0;
                              // if( (iNeighbour == -1) && (jNeighbour == -1) )
                              //   maxZhisto = 100.0;
                              // if( (iNeighbour == -1) && (jNeighbour == 1) )
                              //   maxZhisto = 100.0;
                              // if( (iNeighbour == 1) && (jNeighbour ==0) )
                              //   maxZhisto = 600.0;
                              // if( (iNeighbour == 1) && (jNeighbour == -1) )
                              //   maxZhisto = 100.0;
                              // if( (iNeighbour == 1) && (jNeighbour == 1) )
                              //   maxZhisto = 300.0;
                              // }
                              int tempChannel = mppc[iComptMppc][jComptMppc]->GetDigitizerChannel();
                              var << "(ch"<< tempChannel << "):(" << SumChannels << "):(FloodZ)" ;
                              tree->Draw(var.str().c_str(),CutXYZ + CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName());
                              TGraph2D* tempGraph2d = new TGraph2D(tree->GetSelectedRows(),tree->GetV3(),tree->GetV2(), tree->GetV1());



                              //let's get random...
                              // TH3I* temp3dHisto = new TH3I("temp3dHisto","temp3dHisto",
                              // 20,
                              // tempGraph2d->GetXmin(),
                              // tempGraph2d->GetXmax(),
                              // 20,
                              // tempGraph2d->GetYmin(),
                              // tempGraph2d->GetYmax(),
                              // 100,
                              // tempGraph2d->GetZmin(),
                              // tempGraph2d->GetZmax());


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

                                // std::cout << sname.str().c_str()    << " "
                                          // << tempGraph2d->GetXmin() << " "
                                          // << tempGraph2d->GetXmax() << " "
                                          // << tempGraph2d->GetYmin() << " "
                                          // << tempGraph2d->GetYmax() << " "
                                          // << tempGraph2d->GetZmin() << " "
                                          // << tempGraph2d->GetZmax() << std::endl;
                              var.str("");

                              // ComptonCalibrationHistogram->SetTitle(sname.str().c_str());
                              // ComptonCalibrationHistogram->SetName(sname.str().c_str());
                              ComptonCalibrationHistogram->GetXaxis()->SetTitle("W");
                              ComptonCalibrationHistogram->GetYaxis()->SetTitle("Sum Charge [ADC ch]");
                              ComptonCalibrationHistogram->GetZaxis()->SetTitle("pi [ADC ch]");
                              var.str("");
                              var << "(ch"<< tempChannel << "):(" << SumChannels << "):(FloodZ)" << " >> " << sname.str();
                              tree->Draw(var.str().c_str(),CutXYZ + CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName());
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
                                ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetPoint(i,calibrationGraph->Eval(x[i]),y[i],z[i]);
                                // ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetPointError(i,2.0,sqrt(y[i]),sqrt(z[i]));
                              }
                              ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetTitle(sname.str().c_str());
                              ConvertedComptonCalibration[iComptMppc][jComptMppc]->SetName(sname.str().c_str());
                              ConvertedComptonCalibration[iComptMppc][jComptMppc]->GetXaxis()->SetTitle("W");
                              ConvertedComptonCalibration[iComptMppc][jComptMppc]->GetYaxis()->SetTitle("Sum Charge [ADC ch]");
                              ConvertedComptonCalibration[iComptMppc][jComptMppc]->GetZaxis()->SetTitle("pi [ADC ch]");
                              sname.str("");

                              // sname << "quadFit[" << iComptMppc <<  "][" << jComptMppc <<  "]_" << CurrentCrystal->GetID();
                              // std::cout << std::endl;
                              // std::cout << sname.str().c_str() << std::endl;
                              // TF2 *quadFit = new TF2(sname.str().c_str(),"[0] + [1]*x + [2]*y + [3]*x*x + [4]*x*y + [5]*y*y",0,15,0,6000);
                              // quadFit->SetParameter(0,2311);
                              // quadFit->SetParameter(1,-136.7);
                              // quadFit->SetParameter(2,854);
                              // quadFit->SetParameter(3,-34.1);
                              // quadFit->SetParameter(4,-43.3);
                              // quadFit->SetParameter(5,1.061);
                              // ConvertedComptonCalibration[iComptMppc][jComptMppc]->Fit(quadFit);
                              // sname.str("");
                              // std::cout << std::endl;

                              // interpolation
                              // std::cout << "Running on detector " << iMppc + iNeighbour << "," << jMppc + jNeighbour << std::endl;
                              // sname << "Interpolation["<< iMppc + iNeighbour <<  "][" << jMppc + jNeighbour <<  "]_" << CurrentCrystal->GetID();
                              // interpolationGraph[iMppc + iNeighbour][jMppc + jNeighbour] = new TGraphDelaunay(ComptonCalibation[iMppc + iNeighbour][jMppc + jNeighbour]);
                              // interpolationGraph[iMppc + iNeighbour][jMppc + jNeighbour]->SetTitle(sname.str().c_str());
                              // interpolationGraph[iMppc + iNeighbour][jMppc + jNeighbour]->SetName(sname.str().c_str());
                              // interpolationGraph[iMppc + iNeighbour][jMppc + jNeighbour]->SetMaxIter(100000);
                              // interpolationGraph[iMppc + iNeighbour][jMppc + jNeighbour]->SetMarginBinsContent(0);
                              // interpolationGraph[iMppc + iNeighbour][jMppc + jNeighbour]->ComputeZ(0,0);
                              // interpolationGraph[iMppc + iNeighbour][jMppc + jNeighbour]->FindAllTriangles();
                              // sname.str("");
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
                    // CurrentCrystal->SetInterpolationGraph(interpolationGraph);
                    // if(CurrentCrystal->GetID() == 28)
                    // {
                    //   // TGraph2D *testGraph = new TGraph2D();
                    //   // TH3F *testPi = new TH3F("testPi","testPi",100,0,1,100,0,15000,100,0,6000);
                    //   sname << "Pi for channel 10 - Scint in Crystal " << CurrentCrystal->GetID();
                    //   var.str("");
                    //   int tempChannel = 10;
                    //   // int channel = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
                    //   var << "(ch"<< tempChannel << "):(" << SumChannels << "):(FloodZ)" << " >> " << sname.str();
                    //   tree->Draw(var.str().c_str(),CutXYZ + CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName());
                    //   TGraph2D *testGraph = new TGraph2D(tree->GetSelectedRows(),tree->GetV3(),tree->GetV2(), tree->GetV1());
                    //   testGraph->SetTitle(sname.str().c_str());
                    //   testGraph->SetName(sname.str().c_str());
                    //   testGraph->GetXaxis()->SetTitle("W");
                    //   testGraph->GetYaxis()->SetTitle("Sum Charge [ADC ch]");
                    //   testGraph->GetZaxis()->SetTitle("pi [ADC ch]");
                    //   CurrentCrystal->SetComptonCalibration(testGraph);
                    //   var.str("");
                    //   sname.str("");
                    // }
                    // TH3F *Max_Energy_W = new TH3F("Max_Energy_W","Max_Energy_W",100,0,1,100,0,15000,100,0,6000);
                    // sname << "Compton Calibration Plot - Crystal " << CurrentCrystal->GetID();
                    // Max_Energy_W->SetTitle(sname.str().c_str());
                    // Max_Energy_W->SetName(sname.str().c_str());
                    // Max_Energy_W->GetXaxis()->SetTitle("W");
                    // Max_Energy_W->GetYaxis()->SetTitle("Sum Charge [ADC ch]");
                    // Max_Energy_W->GetZaxis()->SetTitle("Max Charge [ADC ch]");
                    // var.str("");
                    // var << "(ch"<< channel << "):(" << SumChannels << "):(FloodZ)" << " >> " << sname.str();
                    // tree->Draw(var.str().c_str(),CutXYZ + CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName());
                    // CurrentCrystal->SetComptonCalibration(Max_Energy_W);
                    // var.str("");
                    // sname.str("");
                    // //histogram of doi resolution
                    // 		    double sigmaWdistro = userSigmaW;
                    // 		    if(userSigmaW == -1)
                    // 		    {
                    // 		      sigmaWdistro = CurrentCrystal->GetDeltaW();
                    // 		    }
                    // 		    sname << "Doi Resolution Histogram - Crystal " << CurrentCrystal->GetID();
                    // 		    TH1F* derivativeDoiResolution = new TH1F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,0,1);
                    // 		    for(int iPdfHisto = 0 ; iPdfHisto < wHistogramsBins; iPdfHisto++)
                    // 		    {
                    // 		      if(backgroundRun | lateralRun)
                    // 			derivativeDoiResolution->Fill(pdfW->GetBinCenter(iPdfHisto+1),sigmaWdistro*crystalz*pdfW->GetBinContent(iPdfHisto+1));
                    // 		      else
                    // 			derivativeDoiResolution->Fill(pdfW->GetBinCenter(iPdfHisto+1),TMath::Abs(sigmaWdistro*( ( lambda511 * pdfW->GetBinContent(iPdfHisto+1) * (TMath::Exp( -(crystalz/lambda511) ) -1 ) ) / ( 1.0 - (1.0 - TMath::Exp(-(crystalz/lambda511)) )* cumulativeW->GetBinContent(iPdfHisto+1) ) ) ));
                    // 		    }
                    // 		    CurrentCrystal->SetDerivativeDoiResolution(derivativeDoiResolution);
                    // 		    sname.str("");

                    //graph of doi res as a function of z
                    // 		    sname << "Doi Res. vs. Z - Crystal " << CurrentCrystal->GetID();
                    // 		    TGraph* doiResZ = new TGraph();
                    // 		    int pCounter = 0;
                    // 		    for(int iPdfHisto = 0 ; iPdfHisto < wHistogramsBins; iPdfHisto++)
                    // 		    {
                    // 		      doiResZ->SetPoint(pCounter,calibrationGraph->Eval(pdfW->GetBinCenter(iPdfHisto+1)),2.355*derivativeDoiResolution->GetBinContent(iPdfHisto+1));
                    // 		      pCounter++;
                    // 		    }
                    // 		    doiResZ->SetTitle(sname.str().c_str());
                    // 		    doiResZ->SetName(sname.str().c_str());
                    // 		    doiResZ->GetXaxis()->SetTitle("Z [mm]");
                    // 		    doiResZ->GetYaxis()->SetTitle("DoiRes FWHM [mm]");
                    // 		    CurrentCrystal->SetDoiResZ(doiResZ);
                    // 		    sname.str("");



                    //histo of resolutions
                    // 		    sname << "Doi Res. All - " << CurrentCrystal->GetID();
                    // 		    TH1F* resolutions = new TH1F(sname.str().c_str(),sname.str().c_str(),100,0,crystalz);
                    // 		    resolutions->GetXaxis()->SetTitle("Doi Res FWHM [mm]");
                    // 		    for(int iPdfHisto = 0 ; iPdfHisto < wHistogramsBins; iPdfHisto++)
                    // 		    {
                    // 		      resolutions->Fill(2.355*derivativeDoiResolution->GetBinContent(iPdfHisto+1));
                    // 		    }
                    // 		    CurrentCrystal->SetDoiResolutions(resolutions);
                    // 		    //now a bit dummy, sample the curve and find an average
                    // 		    int samplingNumb = 20;
                    // 		    double sumResolutions = 0;
                    // 		    for(int iSample = 0 ; iSample < samplingNumb; iSample++)
                    // 		    {
                    // 		      sumResolutions += doiResZ->Eval(iSample*(crystalz/samplingNumb));
                    // 		    }
                    // 		    sname.str("");
                    // 		    if( sumResolutions/samplingNumb > 0 && sumResolutions/samplingNumb < crystalz )
                    // 		      CurrentCrystal->SetAverageDoiResolution(sumResolutions/samplingNumb);
                    // 		    else
                    // 		      CurrentCrystal->SetAverageDoiResolution(crystalz);

                    //histogram with the difference between the tag calibration and the analytical calibration
                    // 		sname << "Calibration Difference - Crystal " << CurrentCrystal->GetID();
                    // 		var << "((ch" << channel << "/(" << SumChannels << "))*"<<m_cal[CurrentCrystal->GetID()]<<"+"<<q_cal[CurrentCrystal->GetID()]<<")-((ch" << channel << "/(" << SumChannels << "))*"<<m_tag[CurrentCrystal->GetID()]<<"+"<<q_tag[CurrentCrystal->GetID()]<<")>> " << sname.str();
                    // 		TH1F* HistoCalibDiff = new TH1F(sname.str().c_str(),sname.str().c_str(),2500,-20,20);
                    // 		tree->Draw(var.str().c_str(),CutXYZ + CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCutCorrected+triggerPhotopeakCut);
                    // 		HistoCalibDiff->GetXaxis()->SetTitle("Delta Calibration");
                    // 		HistoCalibDiff->GetYaxis()->SetTitle("N");
                    // 		calib_accuracy_File<< CurrentCrystal->GetID() << " " << HistoCalibDiff->GetMean() << " " << HistoCalibDiff->GetRMS() << std::endl;;
                    // 		CurrentCrystal->SetHistoCalibDiff(HistoCalibDiff);
                    // 		var.str("");
                    // 		sname.str("");

                    if(usingRealSimData) // only if this is a sim dataset
                    {
                      long long int nPoints;
                      // a 2d plot of real vs. w, using the same cuts as before
                      if(correctingForDOI)
                        PhotopeakEnergyCut = PhotopeakEnergyCutCorrected;

                      int RealZBins = 100;

                      sname << "Real Z vs. W - Crystal " << CurrentCrystal->GetID();
                      var << "-(RealZ-" << CurrentCrystal->GetDimensionZ()/2.0 << "):FloodZ >> " << sname.str();
                      TH2F* spectrum2dSimDOIplot = new TH2F(sname.str().c_str(),sname.str().c_str(),wHistogramsBins,0,1,RealZBins,0,CurrentCrystal->GetDimensionZ());
                      nPoints = tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCut ,"COLZ");
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

    }
  }
  //----------------------------------------------------------//




  TCanvas* C_spectrum;
  TCanvas* C_multi;
  TCanvas* C_multi_2d;
  TCanvas* C_global;
  TCanvas* C_multi_2;


  //----------------------------------------------------------//
  // Produce some Canvases                                    //
  //----------------------------------------------------------//
  //multicanvases per module 0,0
  //FIXME extend to multiple modules please...
  TCanvas* RawCanvas = new TCanvas("RawSpectra","Rawspectra",1200,800);
  TCanvas* TriggerCanvas = new TCanvas("TriggerSpectra","TriggerSpectra",1200,800);
  TCanvas* FloodHistoCanvas = new TCanvas("FloodHisto","FloodHisto",800,800);
  TCanvas* FloodSeparatedCanvas = new TCanvas("FloodSeparatedCanvas","FloodSeparatedCanvas",800,800);
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
  // <<<<<<< HEAD
  // =======
  // TCanvas* GlobalFlood3DSeparation = new TCanvas("GlobalFlood3DSeparation","GlobalFlood3DSeparation",800,800);

  // >>>>>>> doiTag
  TCanvas* BigSpectraCanvas = new TCanvas("BigSpectra","",800,800);
  TCanvas* BigLYCanvas = new TCanvas("BigLY","",800,800);
  BigSpectraCanvas->Divide(ncrystalsx*nmppcx*nmodulex,ncrystalsy*nmppcy*nmoduley);
  BigLYCanvas->Divide(ncrystalsx*nmppcx*nmodulex,ncrystalsy*nmppcy*nmoduley);
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

  //   TCanvas* GlobalSpherical = new TCanvas("Spherical Plot","Spherical Plot",1200,800);
  //   TCanvas* GlobalCylindricalX = new TCanvas("Cylindrical Plot Theta:X","Cylindrical Plot Theta:X",1200,800);
  //   TCanvas* GlobalCylindricalY = new TCanvas("Cylindrical Plot Theta:Y","Cylindrical Plot Theta:Y",1200,800);
  //canvas for the tagging crystal
  TCanvas* TaggingCanvas = new TCanvas("Tagging Crystal","Tagging Crystal",1200,800);
  //draw canvases
  std::cout << "Saving data to " << outputFileName << " ..." << std::endl;
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      GlobalFlood2DClean->cd();
      module[iModule][jModule]->GetFloodMap2D()->SetTitle(""); // FIXME temporary removed title, for the poster...
      module[iModule][jModule]->GetFloodMap2D()->Draw("COLZ");
      GlobalFlood2D->cd();
      module[iModule][jModule]->GetFloodMap2D()->Draw("COLZ");
      GlobalFlood3D->cd();
      module[iModule][jModule]->GetFloodMap3D()->Draw();
      //      GlobalFlood3DSeparation->cd();
      //      module[iModule][jModule]->GetFloodMap3DSeparation()->Draw();

      //       GlobalSpherical->cd();
      //       module[iModule][jModule]->GetSphericalMap()->Draw("COLZ");
      //       GlobalCylindricalX->cd();
      //       module[iModule][jModule]->GetCylindricalXMap()->Draw("COLZ");
      //       GlobalCylindricalY->cd();
      //       module[iModule][jModule]->GetCylindricalYMap()->Draw("COLZ");
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
        for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
        {
          if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForModular())
          {
            RawCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetRawSpectrum()->Draw();
            TriggerCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrum()->Draw();
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrumHighlighted()->SetFillColor(3);
            mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrumHighlighted()->Draw("same");
            // 	  FloodHistoCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
            // 	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->Draw("COLZ");
            // 	  // 	  TSpectrum2 *peaks2D = new TSpectrum2(ncrystalsx*ncrystalsy,1);
            // 	  // 	  int nfound2D = peaks2D->Search(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D(),1,"col",0.6);
            //
            // 	  FloodHisto3DCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
            // 	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap3D()->Draw("COLZ");
            // 	  // 	  SphericalCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
            // 	  // 	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetSphericalMap()->Draw("COLZ");
            // 	  // 	  CylindricalXCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
            // 	  // 	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCylindricalXMap()->Draw("COLZ");
            // 	  // 	  CylindricalYCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition());
            // 	  // 	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCylindricalYMap()->Draw("COLZ");
            //
            // 	  // temp debug
            //
          }
        }
      }
    }
  }
  //----------------------------------------------------------//



  //----------------------------------------------------------//
  // Global distributions                                     //
  //----------------------------------------------------------//
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

  //   PeakPositionVsIJ->SetStats(1);
  //--Distribution of energy resolutions FHWM
  //histogram
  TH1F *PeakEnergyResolutionDistro = new TH1F("Energy res FWHM","Distribution photopeak energy resolutions FWHM",200,0,1);
  PeakEnergyResolutionDistro->GetXaxis()->SetTitle("Energy Resolution FWHM");
  PeakEnergyResolutionDistro->GetYaxis()->SetTitle("N");
  PeakEnergyResolutionDistro->SetStats(1);

  //   EnergyResolutionVsIJ->SetStats(1);
  //Distribution of FWHM of W plots
  //histogram of fwhm
  //   TH1F *WfwhmDistro = new TH1F("w_fwhm","Distribution of FWHM in W plots",200,0,0.5);
  //   WfwhmDistro->GetXaxis()->SetTitle("W");
  //   WfwhmDistro->GetYaxis()->SetTitle("N");
  //   WfwhmDistro->SetStats(1);
  //2d histogram of fwhm
  //   TH2F *WfwhmVsIJ = new TH2F("w_fwhm vs. i,j","Distribution of FWHM in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  //   WfwhmVsIJ->GetXaxis()->SetTitle("i");
  //   WfwhmVsIJ->GetYaxis()->SetTitle("j");
  //   WfwhmVsIJ->GetZaxis()->SetTitle("w FHWM");
  //   WfwhmVsIJ->GetZaxis()->SetRangeUser(0,0.25);
  //   WfwhmVsIJ->SetStats(1);
  //histogram of rms
  //   TH1F *WrmsDistro = new TH1F("w_rms","Distribution of RMS in W plots",200,0,0.5);
  //   WrmsDistro->GetXaxis()->SetTitle("W");
  //   WrmsDistro->GetYaxis()->SetTitle("N");
  //   WrmsDistro->SetStats(1);
  //2d histogram of rms
  //   TH2F *WrmsVsIJ = new TH2F("w_rms vs. i,j","Distribution of RMS in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  //   WrmsVsIJ->GetXaxis()->SetTitle("i");
  //   WrmsVsIJ->GetYaxis()->SetTitle("j");
  //   WrmsVsIJ->GetZaxis()->SetTitle("w RMS");
  //Distribution of FWHM of W plots
  //   TH1F *Wwidth20perc = new TH1F("w20","Distribution of width at 20% in W plots",200,0,0.5);
  //   Wwidth20perc->GetXaxis()->SetTitle("W");
  //   Wwidth20perc->GetYaxis()->SetTitle("N");
  //   Wwidth20perc->SetStats(1);

  //same plot but just for the "not lateral" channels
  //   TH1F *Wwidth20percCentral = new TH1F("Central w20 - Central Crystals","Distribution of width at 20% in W plots - Central Crystals",200,0,0.5);
  //   Wwidth20percCentral->GetXaxis()->SetTitle("W");
  //   Wwidth20percCentral->GetYaxis()->SetTitle("N");
  //   Wwidth20percCentral->SetStats(1);

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

  //   TH2F *DoiResolutionVsIJ = new TH2F("DOI res FWHM vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  //   DoiResolutionVsIJ->GetXaxis()->SetTitle("i (U axis)");
  //   DoiResolutionVsIJ->GetYaxis()->SetTitle("j (V axis)");
  //   DoiResolutionVsIJ->GetZaxis()->SetTitle("DOI Resolution FWHM [mm]");
  //   DoiResolutionVsIJ->GetXaxis()->SetTitleOffset(1.8);
  //   DoiResolutionVsIJ->GetYaxis()->SetTitleOffset(1.8);
  //   DoiResolutionVsIJ->GetZaxis()->SetTitleOffset(2.2);
  //   DoiResolutionVsIJ->GetZaxis()->SetRangeUser(0,DoiResolutionVsIJmax);

  //   TH2F *DeltaWvsIJ = new TH2F("Delta W vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  //   DeltaWvsIJ->GetXaxis()->SetTitle("i (U axis)");
  //   DeltaWvsIJ->GetYaxis()->SetTitle("j (V axis)");
  //   DeltaWvsIJ->GetZaxis()->SetTitle("Delta W");
  //   DeltaWvsIJ->GetXaxis()->SetTitleOffset(1.8);
  //   DeltaWvsIJ->GetYaxis()->SetTitleOffset(1.8);
  //   DeltaWvsIJ->GetZaxis()->SetTitleOffset(2.2);
  //   DeltaWvsIJ->GetZaxis()->SetRangeUser(0,0.05);

  //   TH2F *mCalVsIJ = new TH2F("Abs(mCal) vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  //   mCalVsIJ->GetXaxis()->SetTitle("i (U axis)");
  //   mCalVsIJ->GetYaxis()->SetTitle("j (V axis)");
  //   mCalVsIJ->GetZaxis()->SetTitle("abs(mCal)");
  //   mCalVsIJ->GetXaxis()->SetTitleOffset(1.8);
  //   mCalVsIJ->GetYaxis()->SetTitleOffset(1.8);
  //   mCalVsIJ->GetZaxis()->SetTitleOffset(2.2);
  //   mCalVsIJ->GetZaxis()->SetRangeUser(0,mCalVsIJmax);

  //   TH2F *RealmCalVsIJ = new TH2F("mCal vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  //   RealmCalVsIJ->GetXaxis()->SetTitle("i (U axis)");
  //   RealmCalVsIJ->GetYaxis()->SetTitle("j (V axis)");
  //   RealmCalVsIJ->GetZaxis()->SetTitle("mCal)");
  //   RealmCalVsIJ->GetXaxis()->SetTitleOffset(1.8);
  //   RealmCalVsIJ->GetYaxis()->SetTitleOffset(1.8);
  //   RealmCalVsIJ->GetZaxis()->SetTitleOffset(2.2);
  //   RealmCalVsIJ->GetZaxis()->SetRangeUser(RealmCalVsIJmin,RealmCalVsIJmax);

  //   TH2F *qCalVsIJ = new TH2F("qCal vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  //   qCalVsIJ->GetXaxis()->SetTitle("i (U axis)");
  //   qCalVsIJ->GetYaxis()->SetTitle("j (V axis)");
  //   qCalVsIJ->GetZaxis()->SetTitle("qCal)");
  //   qCalVsIJ->GetXaxis()->SetTitleOffset(1.8);
  //   qCalVsIJ->GetYaxis()->SetTitleOffset(1.8);
  //   qCalVsIJ->GetZaxis()->SetTitleOffset(2.2);
  //   qCalVsIJ->GetZaxis()->SetRangeUser(0,qCalVsIJmax);


  TH2F *EnergyResolutionVsIJ_corr = new TH2F("Corrected Energy res FWHM vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  EnergyResolutionVsIJ_corr->GetXaxis()->SetTitle("i (U axis)");
  EnergyResolutionVsIJ_corr->GetYaxis()->SetTitle("j (V axis)");
  EnergyResolutionVsIJ_corr->GetZaxis()->SetTitle("Energy Resolution FWHM");
  EnergyResolutionVsIJ_corr->GetXaxis()->SetTitleOffset(1.8);
  EnergyResolutionVsIJ_corr->GetYaxis()->SetTitleOffset(1.8);
  EnergyResolutionVsIJ_corr->GetZaxis()->SetTitleOffset(2.2);
  EnergyResolutionVsIJ_corr->GetZaxis()->SetRangeUser(0,0.3);


  //   TH2F *AverageDoiResolutionVsIJ = new TH2F("Average DOI res FWHM vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  //   AverageDoiResolutionVsIJ->GetXaxis()->SetTitle("i (U axis)");
  //   AverageDoiResolutionVsIJ->GetYaxis()->SetTitle("j (V axis)");
  //   AverageDoiResolutionVsIJ->GetZaxis()->SetTitle("DOI Resolution FWHM [mm]");
  //   AverageDoiResolutionVsIJ->GetXaxis()->SetTitleOffset(1.8);
  //   AverageDoiResolutionVsIJ->GetYaxis()->SetTitleOffset(1.8);
  //   AverageDoiResolutionVsIJ->GetZaxis()->SetTitleOffset(2.2);
  //   AverageDoiResolutionVsIJ->GetZaxis()->SetRangeUser(0,DoiResolutionVsIJmax);

  //   TH2F *Wwidht20percVsIJ = new TH2F("w20 vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  //   Wwidht20percVsIJ->GetXaxis()->SetTitle("i (U axis)");
  //   Wwidht20percVsIJ->GetYaxis()->SetTitle("i (V axis)");
  //   Wwidht20percVsIJ->GetZaxis()->SetRangeUser(0,0.25);
  //   Wwidht20percVsIJ->GetXaxis()->SetTitleOffset(1.8);
  //   Wwidht20percVsIJ->GetYaxis()->SetTitleOffset(1.8);
  //   Wwidht20percVsIJ->GetZaxis()->SetTitleOffset(2.2);
  //   Wwidht20percVsIJ->GetZaxis()->SetTitle("W width at 20%");


  //Distribution of DOI resolutions
  //   TH1F *WDoiDistro = new TH1F("Doi Res FWHM","Distribution of DOI res FWHM",50,0,10);
  //   WDoiDistro->GetXaxis()->SetTitle("DOI resolution FWHM [mm]");
  //   //   WDoiDistro->GetYaxis()->SetTitle("N");
  //   WDoiDistro->SetStats(1);

  //Distribution of DOI resolutions
  //   TH1F *WDoiDistroCentral = new TH1F("Central Doi Res FWHM","Central distribution of DOI res FWHM",50,0,10);
  //   //   WDoiDistroCentral->GetXaxis()->SetTitle("DOI resolution FWHM [mm]");
  //   WDoiDistroCentral->GetYaxis()->SetTitle("N");
  //   WDoiDistroCentral->SetStats(1);

  //Distribution of AVERAGE DOI resolutions
  //   TH1F *AverageDoiDistro = new TH1F("Average Doi Res FWHM","Average distribution of DOI res FWHM",50,0,10);
  //   AverageDoiDistro->GetXaxis()->SetTitle("DOI resolution FWHM [mm]");
  //   WDoiDistro->GetYaxis()->SetTitle("N");
  //   AverageDoiDistro->SetStats(1);

  //Distribution of AVERAGE DOI resolutions -CENTRAL
  //   TH1F *AverageDoiDistroCentral = new TH1F("Central - Average Doi Res FWHM","Central - Average distribution of DOI res FWHM",50,0,10);
  //   AverageDoiDistroCentral->GetXaxis()->SetTitle("DOI resolution FWHM [mm]");
  //   WDoiDistro->GetYaxis()->SetTitle("N");
  //   AverageDoiDistroCentral->SetStats(1);


  //Distribution of fit exp tau for W plots
  //   TH1F *WtauFit = new TH1F("Exp slope W","Distribution of exp slopes in W plots",1000,0,1);
  //   WtauFit->GetXaxis()->SetTitle("Tau");
  //   WtauFit->GetYaxis()->SetTitle("N");
  //   WtauFit->SetStats(1);


  //   TH2F *WtauFitVsIJ = new TH2F("Exp slopes W vs. i,j","Distribution of exp slopes in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  //   WtauFitVsIJ->GetXaxis()->SetTitle("i");
  //   WtauFitVsIJ->GetYaxis()->SetTitle("i");
  //   WtauFitVsIJ->GetZaxis()->SetTitle("tau");




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
      //       std::cout << ModuleDirStream.str() << std::endl;
      directory[iModule+jModule][0][0] = fPlots->mkdir(ModuleDirStream.str().c_str());
      directory[iModule+jModule][0][0]->cd();
      GlobalFlood2D->Write();
      GlobalFlood3D->Write();

      if(usingRealSimData)
      {
        C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
        C_spectrum->SetName(module[iModule][jModule]->GetFloodMap2DSingleCrystalHit()->GetName());
        C_spectrum->cd();
        module[iModule][jModule]->GetFloodMap2DSingleCrystalHit()->Draw("COLZ");
        C_spectrum->Write();
        delete C_spectrum;

      }

      if(usingTaggingBench)
      {
        C_TaggingCrystalSpectrum->cd();
        TaggingCrystalSpectrum->Draw();
        TriggerSpectrumHighlight->Draw("same");
        C_TaggingCrystalSpectrum->Write();
        // GlobalFlood3DSeparation->Write();

      }

      RawCanvas->Write();
      TriggerCanvas->Write();
      BigSpectraCanvas->Write();
      if(lightYieldComputation) BigLYCanvas->Write();


      //       C_global = new TCanvas("C_global","C_global",800,800);
      //       name = "3D Cuts - " + ModuleDirStream.str();
      //       C_global->SetName(name);
      //       C_global->SetTitle(name);
      //       int colorcounter = 0;

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
                Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
                if(CurrentCrystal->GetIsOnForModular())
                {
                  CrystalDirStream << "Crystal " <<  CurrentCrystal->GetID();
                  directory[iModule+jModule][(iMppc+jMppc)+1][(iCry+jCry)+1] = directory[iModule+jModule][(iMppc+jMppc)+1][0]->mkdir(CrystalDirStream.str().c_str());
                  directory[iModule+jModule][(iMppc+jMppc)+1][(iCry+jCry)+1]->cd();
                  if(CurrentCrystal->CrystalIsOn() /*| usingRealSimData*/) // save data only if the crystal was found
                  {
                    //create a pointer for the current crystal (mainly to make the code more readable)
                    Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
                    //fill the global distributions histograms
                    if(!backgroundRun)
                    {
                      PeakPositionDistro->Fill(CurrentCrystal->GetPhotopeakPosition());
                      PeakEnergyResolutionDistro->Fill(CurrentCrystal->GetPhotopeakEnergyResolution());

                      if(correctingForDOI)
                      {
                        PeakEnergyResolutionDistro_corr->Fill(CurrentCrystal->GetPhotopeakEnergyResolutionCorrected());
                      }
                      // 		WfwhmDistro->Fill(CurrentCrystal->GetWfwhm());
                      // 		WDoiDistro->Fill( (15.0/CurrentCrystal->GetWfwhm())*0.0158); // FIXME CAREFUL: here the 0.0158 value is hardcoded and taken from the sigma of W distros in DOI bench setup. 15.0 is the length of the crystals in mm.

                      PeakPositionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakPosition());
                      EnergyResolutionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakEnergyResolution());
                      EnergyResolutionVsIJ_corr->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakEnergyResolutionCorrected());
                    }

                    if(lightYieldComputation)
                    {
                      LYDistro->Fill(CurrentCrystal->GetLY());
                      LYVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetLY());
                    }

                    // 		    if(CurrentCrystal->GetDoiResolutionFWHM() > 0 && CurrentCrystal->GetDoiResolutionFWHM() < crystalz )
                    // 		      DoiResolutionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetDoiResolutionFWHM());
                    // 		    else
                    // 		      DoiResolutionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),0);
                    // 		    WDoiDistro->Fill(CurrentCrystal->GetDoiResolutionFWHM());


                    // 		    if(CurrentCrystal->GetAverageDoiResolution() > 0 && CurrentCrystal->GetAverageDoiResolution() < crystalz )
                    // 		    AverageDoiResolutionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetAverageDoiResolution());
                    // 		    AverageDoiDistro->Fill(CurrentCrystal->GetAverageDoiResolution());
                    // 		    else
                    // 		      AverageDoiResolutionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),0);

                    // 		    DeltaWvsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetDeltaW());
                    // 		    mCalVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),std::abs(CurrentCrystal->GetMcal()));
                    // 		    RealmCalVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetMcal());
                    // 		    qCalVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetQcal());

                    // 		WfwhmVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetWfwhm());
                    // 		WrmsDistro->Fill(CurrentCrystal->GetWrms());
                    // 		WrmsVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetWrms());
                    // 		Wwidht20percVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetWwidth20perc());
                    // 		Wwidth20perc->Fill(CurrentCrystal->GetWwidth20perc());


                    if( ((iModule*nmppcx)+iMppc) > 0 && (((iModule*nmppcx)+iMppc) < nmppcx -1) && ((jModule*nmppcy)+jMppc) > 0 && (((jModule*nmppcy)+jMppc) < nmppcy -1 ))
                    {
                      // 		  Wwidth20percCentral->Fill(CurrentCrystal->GetWwidth20perc());
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


                      // 		      AverageDoiDistroCentral->Fill(CurrentCrystal->GetAverageDoiResolution());
                      // 		      WDoiDistroCentral->Fill(CurrentCrystal->GetDoiResolutionFWHM());

                    }

                    //draw also the 3d cuts in the common canvas of this mppc
                    C_multi->cd();
                    CurrentCrystal->GetFloodMap3D()->SetMarkerColor(counter);
                    CurrentCrystal->GetFloodMap3D()->SetFillColor(counter);
                    legend->AddEntry(CurrentCrystal->GetFloodMap3D(),CurrentCrystal->GetFloodMap3D()->GetName(),"f");
                    CurrentCrystal->GetFloodMap3D()->Draw("same");
                    legend->Draw();
                    // 		counter++;

                    //same as above, but the 2D maxPeak
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

                    // 		C_global->cd();
                    // 		CurrentCrystal->GetFloodMap3D()->SetMarkerColor(colorcounter);
                    // 		CurrentCrystal->GetFloodMap3D()->SetFillColor(colorcounter);
                    // 		CurrentCrystal->GetFloodMap3D()->Draw("same");
                    // 		colorcounter++;

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

                    // spectrum without highligth
                    // 		    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    // 		    TString title = "mod_";
                    // 		    title += CurrentCrystal->GetSpectrum()->GetName() ;
                    // 		    C_spectrum->SetName(title);
                    // 		    C_spectrum->SetTitle("");
                    // 		    C_spectrum->cd();
                    // 		    CurrentCrystal->GetSpectrum()->SetFillStyle(3001);
                    // 		    CurrentCrystal->GetSpectrum()->SetFillColor(kBlue);
                    // 		    CurrentCrystal->GetSpectrum()->Draw();
                    // 		    C_spectrum->Write();
                    // 		    delete C_spectrum;

                    //w histo
                    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    C_spectrum->SetName(CurrentCrystal->GetHistoW()->GetName());
                    C_spectrum->cd();
                    CurrentCrystal->GetHistoW()->Draw();
                    if(usingTaggingBench) CurrentCrystal->GetHistoWfit()->Draw("same");
                    C_spectrum->Write();
                    delete C_spectrum;

                    // adc versus time
                    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    C_spectrum->SetName(CurrentCrystal->GetVersusTime()->GetName());
                    C_spectrum->cd();
                    CurrentCrystal->GetVersusTime()->Draw("COLZ");
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

                    // <<<<<<< HEAD
                    if(!backgroundRun)
                    {
                      if(correctingForDOI)
                      {
                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                        // 			C_spectrum->SetName(CurrentCrystal->GetSlicesMean()->GetName());
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
                    // 		    CurrentCrystal->GetThetaFit()->Draw("same");
                    //instead of drawing the real fitted function, draw another gaussian function with same parameters but in the 0-1 range
                    // 		    TF1 *gaussDraw = new TF1("gaussDraw","[0]*exp(-0.5*((x-[1])/[2])**2)",0,1);
                    // 		    gaussDraw->SetParameter(0,CurrentCrystal->GetDeltaWfit()->GetParameter(0));
                    // 		    gaussDraw->SetParameter(1,CurrentCrystal->GetDeltaWfit()->GetParameter(1));
                    // 		    gaussDraw->SetParameter(2,CurrentCrystal->GetDeltaWfit()->GetParameter(2));
                    // 		    gaussDraw->SetLineColor(3);
                    // 		CurrentCrystal->GetDeltaWfit()->Draw("same");
                    // 		    gaussDraw->Draw("same");
                    // 		    TF1 *gaussDraw_2 = new TF1("gaussDraw_2","[0]*exp(-0.5*((x-[1])/[2])**2)",0,1);
                    // 		    gaussDraw_2->SetParameter(0,CurrentCrystal->GetDeltaWfit_2()->GetParameter(0));
                    // 		    gaussDraw_2->SetParameter(1,CurrentCrystal->GetDeltaWfit_2()->GetParameter(1));
                    // 		    gaussDraw_2->SetParameter(2,CurrentCrystal->GetDeltaWfit_2()->GetParameter(2));
                    // 		    gaussDraw_2->SetLineColor(2);
                    // 		    gaussDraw_2->Draw("same");
                    C_spectrum->Write();
                    delete C_spectrum;

                    //w histo corrected smooth
                    // 		    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    // 		    C_spectrum->SetName(CurrentCrystal->GetHistoWCorrectedSmooth()->GetName());
                    // 		    C_spectrum->cd();
                    // 		    CurrentCrystal->GetHistoWCorrectedSmooth()->Draw();
                    // 		    C_spectrum->Write();
                    // 		    delete C_spectrum;

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

                    //calibration graph
                    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    C_spectrum->SetName(CurrentCrystal->GetCalibrationGraph()->GetName());
                    C_spectrum->cd();
                    CurrentCrystal->GetCalibrationGraph()->Draw("AL");
                    C_spectrum->Write();
                    delete C_spectrum;

                    //compton w(z)
                    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    C_spectrum->SetName(CurrentCrystal->GetWZgraph()->GetName());
                    C_spectrum->cd();
                    CurrentCrystal->GetWZgraph()->Draw("AL");
                    C_spectrum->Write();
                    delete C_spectrum;

                    CurrentCrystal->GetZYCut()->Write();
                    CurrentCrystal->GetZXCut()->Write();
                    // CurrentCrystal->GetOneHisto()->Write();


                    if(usingTaggingBench)
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




                    // compton calibration
                    if(comptonAnalysis)
                    {
                      for(int iComptMppc = 0; iComptMppc < nmppcx ; iComptMppc++)
                      {
                        for(int jComptMppc = 0; jComptMppc < nmppcy ; jComptMppc++)
                        {
                      // for(int iNeighbour = -1 ; iNeighbour < 2 ; iNeighbour++)
                      // {
                      //   for(int jNeighbour = -1 ; jNeighbour < 2 ; jNeighbour++)
                      //   {
                      //     if( (iMppc + iNeighbour >= 0) && (iMppc + iNeighbour <= nmppcx) )
                      //     {
                      //       if( (jMppc + jNeighbour >= 0) && (jMppc + jNeighbour <= nmppcy) )
                      //       {
                              // C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
                              TGraph2D ***tempGraph = CurrentCrystal->GetComptonCalibration();
                              // C_spectrum->SetName(tempGraph[iMppc + iNeighbour][jMppc + jNeighbour]->GetName());
                              // C_spectrum->cd();
                              // tempGraph[iMppc + iNeighbour][jMppc + jNeighbour]->Draw("AP");
                              tempGraph[iComptMppc][jComptMppc]->Write();
                              // C_spectrum->Write();
                              // delete C_spectrum;

                              // TGraphDelaunay ***tempDelaunay = CurrentCrystal->GetInterpolationGraph();
                              // tempDelaunay[iMppc + iNeighbour][jMppc + jNeighbour]->Write();
                              TGraph2D ***tempCorrGraph = CurrentCrystal->GetConvertedComptonCalibration();
                              tempCorrGraph[iComptMppc][jComptMppc]->Write();


                              // tempHisto[iMppc + iNeighbour][jMppc + jNeighbour]->Write();


                              // TH3I ***tempHisto = CurrentCrystal->GetComptonCalibrationHistogram();
                              // // tempHisto[iMppc + iNeighbour][jMppc + jNeighbour]->Write();
                              // C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
                              // C_spectrum->SetName(tempHisto[iMppc + iNeighbour][jMppc + jNeighbour]->GetName());
                              // C_spectrum->cd();
                              // tempHisto[iMppc + iNeighbour][jMppc + jNeighbour]->Draw();
                              // C_spectrum->Write();
                              // delete C_spectrum;
                            }
                          }
                      //       }
                      //     }
                      //   }
                      // }
                      for(int iCompt = 0 ; iCompt < CurrentCrystal->GetNumOfComptonHisto(); iCompt++)
                      {
                        C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
                        C_spectrum->SetName(CurrentCrystal->GetComptonHistogram(iCompt)->GetName());
                        C_spectrum->cd();
                        CurrentCrystal->GetComptonHistogram(iCompt)->Draw();
                        C_spectrum->Write();
                        delete C_spectrum;
                      }
                    }


                    // C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    // C_spectrum->SetName(CurrentCrystal->GetComptonCalibration()->GetName());
                    // C_spectrum->cd();
                    // CurrentCrystal->GetComptonCalibration()->Draw("AP");
                    // C_spectrum->Write();
                    // delete C_spectrum;



                    // 		    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    // 		    C_spectrum->SetName(CurrentCrystal->GetDerivativeDoiResolution()->GetName());
                    // 		    C_spectrum->cd();
                    // 		    CurrentCrystal->GetDerivativeDoiResolution()->Draw();
                    // 		    C_spectrum->Write();
                    // 		    delete C_spectrum;

                    // 		    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    // 		    C_spectrum->SetName(CurrentCrystal->GetDoiResZ()->GetName());
                    // 		    C_spectrum->cd();
                    // 		    CurrentCrystal->GetDoiResZ()->Draw("AL");
                    // 		    C_spectrum->Write();
                    // 		    delete C_spectrum;

                    // 		    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    // 		    C_spectrum->SetName(CurrentCrystal->GetDoiResolutions()->GetName());
                    // 		    C_spectrum->cd();
                    // 		    CurrentCrystal->GetDoiResolutions()->Draw();
                    // 		    C_spectrum->Write();
                    // 		    delete C_spectrum;
                    // =======
                    //
                    //
                    // 		if(correctingForDOI)
                    // 		{
                    // 		  C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    // 		  C_spectrum->SetName(CurrentCrystal->GetSlicesMean()->GetName());
                    // 		  C_spectrum->cd();
                    // 		  CurrentCrystal->GetSlicesMean()->Draw();
                    // 		  CurrentCrystal->GetSlicesMeanFit()->Draw("same");
                    // 		  C_spectrum->Write();
                    // 		  delete C_spectrum;
                    //
                    // 		  C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    // 		  C_spectrum->SetName(CurrentCrystal->GetCorrectedSpectrum()->GetName());
                    // 		  C_spectrum->cd();
                    // 		  CurrentCrystal->GetCorrectedSpectrum()->Draw();
                    // 		  CurrentCrystal->GetHighlightedSpectrumCorrected()->SetFillColor(3);
                    // 		  CurrentCrystal->GetHighlightedSpectrumCorrected()->Draw("same");
                    // 		  CurrentCrystal->GetFitCorrected()->Draw("same");
                    // 		  C_spectrum->Write();
                    // 		  delete C_spectrum;
                    //
                    //
                    //
                    // 		  //w histo corrected
                    // 		  C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    // 		  C_spectrum->SetName(CurrentCrystal->GetHistoWCorrected()->GetName());
                    // 		  C_spectrum->cd();
                    // 		  CurrentCrystal->GetHistoWCorrected()->Draw();
                    // 		  C_spectrum->Write();
                    // 		  delete C_spectrum;
                    //
                    // 		}
                    // >>>>>>> doiTag
                    //density histo
                    // 		    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
                    // 		    C_spectrum->SetName(CurrentCrystal->GetDensityHisto()->GetName());
                    // 		    C_spectrum->cd();
                    // 		    CurrentCrystal->GetDensityHisto()->Draw();
                    // 		    C_spectrum->Write();
                    // 		    delete C_spectrum;


                    if(usingRealSimData)
                    {
                      // 		      WtauFit->Fill(CurrentCrystal->GetSimFit()->GetParameter(1));
                      // 		      WtauFitVsIJ->Fill(CurrentCrystal->GetI() ,CurrentCrystal->GetJ() , CurrentCrystal->GetSimFit()->GetParameter(1));

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

            directory[iModule+jModule][(iMppc+jMppc)+1][0]->cd();
            //   	  pt->Draw();
            // 	  nameMppc = "3D Cuts - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
            // 	  C_multi->SetName(name);
            // 	  C_multi->SetTitle(name);
            // 	  C_multi->Update();
            C_multi->Write();
            C_multi_2d->Write();

            // 	  C_multi_2->Write();
            // 	  C_graph->Write();
            // 	  delete C_graph;
            delete C_multi_2d;
            delete C_multi;
          }
        }
      }


      directory[iModule+jModule][0][0]->cd(); // go back to main directory
      //       C_global->Write();
      //       delete C_global;
      //write the summary histos and canvases, that were filled during file saving



      //       TCanvas *C_WfwhmVsIJ = new TCanvas("C_WfwhmVsIJ","C_WfwhmVsIJ",800,800);
      //       C_WfwhmVsIJ->SetName(WfwhmVsIJ->GetName());
      //       C_WfwhmVsIJ->cd();
      //       WfwhmVsIJ->Draw("LEGO2");
      //       C_WfwhmVsIJ->Write();

      //       TCanvas *C_WrmsVsIJ = new TCanvas("C_WrmsVsIJ","C_WrmsVsIJ",800,800);
      //       C_WrmsVsIJ->SetName(WrmsVsIJ->GetName());
      //       C_WrmsVsIJ->cd();
      //       WrmsVsIJ->Draw("LEGO2");
      //       C_WrmsVsIJ->Write();

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

      //       TCanvas *C_DoiResolutionVsIJ = new TCanvas("C_DoiResolutionVsIJ","C_DoiResolutionVsIJ",800,800);
      //       C_DoiResolutionVsIJ->SetName(DoiResolutionVsIJ->GetName());
      //       C_DoiResolutionVsIJ->cd();
      //       DoiResolutionVsIJ->Draw("LEGO2");
      //       C_DoiResolutionVsIJ->SetLeftMargin(0.15);
      //       C_DoiResolutionVsIJ->Write();

      //       TCanvas *C_DeltaWvsIJ = new TCanvas("C_DeltaWvsIJ","C_DeltaWvsIJ",800,800);
      //       C_DeltaWvsIJ->SetName(DeltaWvsIJ->GetName());
      //       C_DeltaWvsIJ->cd();
      //       DeltaWvsIJ->Draw("LEGO2");
      //       C_DeltaWvsIJ->SetLeftMargin(0.15);
      //       C_DeltaWvsIJ->Write();

      //       TCanvas *C_mCalVsIJ = new TCanvas("C_mCalVsIJ","C_mCalVsIJ",800,800);
      //       C_mCalVsIJ->SetName(mCalVsIJ->GetName());
      //       C_mCalVsIJ->cd();
      //       mCalVsIJ->Draw("LEGO2");
      //       C_mCalVsIJ->SetLeftMargin(0.15);
      //       C_mCalVsIJ->Write();

      //       TCanvas *C_AverageDoiResolutionVsIJ = new TCanvas("C_AverageDoiResolutionVsIJ","C_AverageDoiResolutionVsIJ",800,800);
      //       C_AverageDoiResolutionVsIJ->SetName(AverageDoiResolutionVsIJ->GetName());
      //       C_AverageDoiResolutionVsIJ->cd();
      //       AverageDoiResolutionVsIJ->Draw("LEGO2");
      //       C_AverageDoiResolutionVsIJ->SetLeftMargin(0.15);
      //       C_AverageDoiResolutionVsIJ->Write();


      //       TCanvas *C_qCalVsIJ = new TCanvas("C_qCalVsIJ","C_qCalVsIJ",800,800);
      //       C_qCalVsIJ->SetName(qCalVsIJ->GetName());
      //       C_qCalVsIJ->cd();
      //       qCalVsIJ->Draw("LEGO2");
      //       C_qCalVsIJ->SetLeftMargin(0.15);
      //       C_qCalVsIJ->Write();

      //       TCanvas *C_RealmCalVsIJ = new TCanvas("C_RealmCalVsIJ","C_RealmCalVsIJ",800,800);
      //       C_RealmCalVsIJ->SetName(RealmCalVsIJ->GetName());
      //       C_RealmCalVsIJ->cd();
      //       RealmCalVsIJ->Draw("LEGO2");
      //       C_RealmCalVsIJ->SetLeftMargin(0.15);
      //       C_RealmCalVsIJ->Write();

      //       TCanvas *C_Wwidht20percVsIJ = new TCanvas("C_Wwidht20percVsIJ","C_Wwidht20percVsIJ",800,800);
      //       C_Wwidht20percVsIJ->SetName(Wwidht20percVsIJ->GetName());
      //       C_Wwidht20percVsIJ->cd();
      //       Wwidht20percVsIJ->Draw("LEGO2");
      //       C_Wwidht20percVsIJ->SetLeftMargin(0.15);
      //       C_Wwidht20percVsIJ->Write();

      //       gStyle->SetOptStat(1);

      //       if(usingRealSimData)
      //       {
      // 	WtauFit->Write();
      // 	TCanvas *C_WtauFitVsIJ = new TCanvas("C_WtauFitVsIJ","C_WtauFitVsIJ",800,800);
      // 	C_WtauFitVsIJ->SetName(WtauFitVsIJ->GetName());
      // 	C_WtauFitVsIJ->cd();
      // 	WtauFitVsIJ->Draw("LEGO2");
      // 	C_WtauFitVsIJ->Write();
      //       }

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
      //       AverageDoiDistro->Write();
      //       AverageDoiDistroCentral->Write();
      //       WDoiDistro->Write();
      //       WDoiDistroCentral->Write();
    }
  }
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

  delete crystal;
  delete mppc;
  delete module;

  return 0;
}

// Double_t thetaFunction(Double_t *x, Double_t *par)
// {
//   Double_t f;
//   Float_t xx =x[0];
//   if (xx>par[0] && xx<par[1])
//   {
//     f = par[2];
//   }
//   else
//   {
//     f = 0;
//   }
//   return f;
// }
