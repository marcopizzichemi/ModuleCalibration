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

// #include <omp.h>

#define ENERGY_RESOLUTION 0.12
#define ENERGY_RESOLUTION_SATURATION_CORRECTION 0.25





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
  input.FillTree();                // create the TTree that will be used in analysis
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
  int histo3DchannelBin         = config.read<int>("histo3DchannelBin",100);            // number of bins of the 3D flood histograms, for single channels
  int histo3DglobalBins         = config.read<int>("histo3DglobalBins");            // number of bins of the 3D flood histograms, for entire module
  int taggingPeakMin            = config.read<int>("taggingPeakMin",8000);          // min range of tagging crystal photopeak, in ADC channels - to help TSpectrum
  int taggingPeakMax            = config.read<int>("taggingPeakMax",12000);         // max range of tagging crystal photopeak, in ADC channels - to help TSpectrum
  int clusterLevelPrecision     = config.read<int>("clusterLevelPrecision",10);     // precision of the level search when separating the cluster of 3D points
  bool saveAnalysisTree         = config.read<bool>("saveAnalysisTree");            // choice to save or not the analysis TTree, in a file temp.root
  float taggingPosition         = config.read<float>("taggingPosition");            // position of the tagging bench in mm 
  bool usingTaggingBench        = config.read<bool>("usingTaggingBench");           // true if the input is using tagging bench, false if not
  int taggingCrystalChannel     = config.read<int>("taggingCrystalChannel");        // input channel where the tagging crystal information is stored
  bool correctingSaturation     = config.read<bool>("correctingSaturation");        // true if saturation correction is applied, false if it's not
  bool correctingForDOI         = config.read<bool>("correctingForDOI",0);            // true if the energy correction using DOI info is computed    
  float energyResolution        = config.read<float>("expectedEnergyResolution",0); // energy resolution input by the user, if any, otherwise 0
  bool usingRealSimData         = config.read<bool>("usingRealSimData",0);
  float moduleLateralSideX      = config.read<float>("moduleLateralSideX",7.0); //
  float moduleLateralSideY      = config.read<float>("moduleLateralSideY",7.0); // 
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
  //  Plots and spectra                                       //
  //----------------------------------------------------------//
  // get the TTree, to plot the spectra
  TTree* tree = input.GetTree();     
  // prepare spectra
  //TH1F* spectrum;
  //TH2F* spectrum2d;
  //TH3I* spectrum3d;
  //TGraph* graph; 
//   TGraph2D* graph2D;
//   TCanvas* C_graph;
  
  // doi bench specific part
  std::ofstream doiFile;
  TCut triggerPhotopeakCut = "" ;
  TH1F* TaggingCrystalSpectrum;
  TH1F *TriggerSpectrumHighlight;
  if(usingTaggingBench)
  {
    doiFile.open("doiData.txt", std::ofstream::out);
  }
  
  // angles for separating crystals, translations for having a nice 2d plot  
//   double lateralQ1;
//   double lateralQ2;
//   double lateralDeltaU;
//   double lateralDeltaV;
//   double lateralRescaleRL;
//   double lateralRescaleTB;
  //corners
//   double cornerQ1;
//   double cornerQ2;
//   double cornerDeltaU;
//   double cornerDeltaV;
//   double cornerRescale;
  
//   TString name;
  
  TCut SingleCrystalInteraction;
  TCut SingleEnergyDeposition;
  
  if(usingRealSimData) // only if this is a sim dataset
  {
    SingleCrystalInteraction = "CrystalsHit == 1";
    SingleEnergyDeposition = "NumbOfInteractions == 1";
  }
  // Loop on modules, mppcs and crystal
  
  
  
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      std::stringstream CutXYZstream;
      CutXYZstream << "FloodX > " << -moduleLateralSideX << " && FloodX < " << moduleLateralSideX  << "&& FloodY > " << -moduleLateralSideY <<   " && FloodY < " << moduleLateralSideY << "&& FloodZ > 0 && FloodZ < 1";
      TCut CutXYZ = CutXYZstream.str().c_str();
      std::cout << "Generating global spectra..." << std::endl;
      
      TString nameModule;
      std::stringstream varModule;
      
//       TH3I* spectrum3dModule;
      
      // GLOBAL SPECTRA
      // Flood histogram
      
      nameModule = "Flood Histogram 2D - " + module[iModule][jModule]->GetName();
      varModule << "FloodY:FloodX >> " << nameModule; 
//       std::cout << nameModule << " ... ";
      TH2F *spectrum2dModule = new TH2F(nameModule,nameModule,histo2DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo2DglobalBins,-moduleLateralSideY,moduleLateralSideY);
      tree->Draw(varModule.str().c_str(),"","COLZ");
      spectrum2dModule->SetName(nameModule); 
      spectrum2dModule->SetTitle(nameModule);
      spectrum2dModule->GetXaxis()->SetTitle("U");
      spectrum2dModule->GetYaxis()->SetTitle("V");
      module[iModule][jModule]->SetFloodMap2D(spectrum2dModule);
      varModule.str("");
//       std::cout << " done" << std::endl;
//       delete spectrum2dModule;
      
      //3D plot
      nameModule = "Flood Histogram 3D - Module " + module[iModule][jModule]->GetName();
      varModule << "FloodZ:FloodY:FloodX >> " << nameModule; 
//       std::cout << nameModule << " ... ";
      TH3I* spectrum3dModule = new TH3I(nameModule,nameModule,histo3DglobalBins,-moduleLateralSideX,moduleLateralSideX,histo3DglobalBins,-moduleLateralSideY,moduleLateralSideY,histo3DglobalBins,0,1);
      tree->Draw(varModule.str().c_str(),CutXYZ);
      spectrum3dModule->GetXaxis()->SetTitle("U");
      spectrum3dModule->GetYaxis()->SetTitle("V");
      spectrum3dModule->GetZaxis()->SetTitle("W");
      module[iModule][jModule]->SetFloodMap3D(spectrum3dModule);
      varModule.str("");
//       std::cout << " done" << std::endl;
//       delete spectrum3dModule;
      
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
//       #pragma omp parallel for
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  //these loop will run in parallel, so some variables need to exist only inside the loop 
	  TString name;
	  std::stringstream var,cut,sname; 
// 	  TH1F* spectrum;
          //TH2F* spectrum2d;
//           TH3I* spectrum3d;
          //TGraph* graph;
	  //FIXME careful, at the moment it works only because there's one module
	  // it should be fixed for more modules by using the same mppc[][] logic used for the crystals, below
	  int channel = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
	  std::cout << "Generating spectra for MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " ..." << std::endl;
	  cut << "TriggerChannel == " << channel  ;
	  TCut CutTrigger = cut.str().c_str();
	  cut.str("");
	  //same as the global ones, but selecting on TriggerChannel
	  // raw spectrum
	  name = "Raw Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() + " - Module " + module[iModule][jModule]->GetName();
// 	  channel = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
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
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTriggerSpectrum(spectrumTrigger);
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
// 	  spectrum3d = new TH3I("spectrum3d","spectrum3d",histo3DchannelBin,-moduleLateralSideX,moduleLateralSideX,histo3DchannelBin,-moduleLateralSideY,moduleLateralSideY,histo3DchannelBin,0,1);
	  //little trick to try and use less bins
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
	  TH3I* spectrum3dMPPC = new TH3I(name,name,histo3DchannelBin,minX3Dplot,maxX3Dplot,histo3DchannelBin,minY3Dplot,maxY3Dplot,histo3DchannelBin,0,1);//FIXME temp
	  tree->Draw(var.str().c_str(),CutXYZ+CutTrigger);
	  spectrum3dMPPC->GetXaxis()->SetTitle("U");
	  spectrum3dMPPC->GetYaxis()->SetTitle("V");
	  spectrum3dMPPC->GetZaxis()->SetTitle("W");
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetFloodMap3D(spectrum3dMPPC);
	  var.str("");
// 	  delete spectrum3d;
	  // automatic crystal finder on the mppc
	  // it runs always, unless the user has set onlyuserinput
	  TCutG**** cutg; // prepare the graphical cuts
	  //const int numbOfCrystals = 4;
	  cutg = new TCutG***[2]; // two planes of cuts, their intersection will create a 3d cut
	  for(int iCut =0 ; iCut < 2 ; iCut++)
	  {
	    cutg[iCut] = new TCutG**[ncrystalsx];
	    for(int iCry = 0; iCry < ncrystalsx ; iCry++)
	    {
	      cutg[iCut][iCry] = new TCutG*[ncrystalsy];
	    }
	  }
	  bool found = false;
	  if(!onlyuserinput)
	  {
	    if(usingTaggingBench) //if it's a tagging bench run, check first if the mppc is on for DOI bench measurement
	    {
	      if(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetIsOnForDoi())
	      {
		found = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->FindCrystalCuts(cutg,histo3DchannelBin,clusterLevelPrecision);
		// 		mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->Find2Dpeaks(ncrystalsx*ncrystalsy,mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D());
	      }
	    }
	    else
	    {
	      found = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->FindCrystalCuts(cutg,histo3DchannelBin,clusterLevelPrecision);
	    // 	      mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->Find2Dpeaks(ncrystalsx*ncrystalsy,mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D());
	    }
	  }
	  
	  
	  
	  // run on all the possible crystals (i.e. all the crystals coupled to this mppc)
	  
	  for(int iCry = 0; iCry < ncrystalsx ; iCry++)
	  {
	    for(int jCry = 0; jCry < ncrystalsy ; jCry++)
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
		//automatically look for the 511Kev peak to find the photopeak energy cut
		//find peaks in each crystal spectrum, with TSpectrum
		
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
		streamEnergyCut << SumChannels << " > " << gauss->GetParameter(1) - 2.0*std::abs(gauss->GetParameter(2)) << " && " << SumChannels << " < " << gauss->GetParameter(1) + 4.0*std::abs(gauss->GetParameter(2));
		TCut PhotopeakEnergyCutCorrected;
		TCut PhotopeakEnergyCut  = streamEnergyCut.str().c_str(); 
		CurrentCrystal->SetSpectrum(spectrumCharge);
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
		//-----------------------------------------------------------------------
		
		
		//w histogram with cut on crystal events, xyz and trigger channel and cut on photopeak
		sname << "W histogram - Crystal " << CurrentCrystal->GetID();
		var << "(ch" << channel << "/(" << SumChannels << ")) >> " << sname.str();
		TH1F* spectrumHistoW = new TH1F(sname.str().c_str(),sname.str().c_str(),250,0,1);	  
		tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCut+triggerPhotopeakCut);
		spectrumHistoW->GetXaxis()->SetTitle("W");
		spectrumHistoW->GetYaxis()->SetTitle("N");
		int bin1 = spectrumHistoW->FindFirstBinAbove(spectrumHistoW->GetMaximum()/2.0);
		int bin2 = spectrumHistoW->FindLastBinAbove(spectrumHistoW->GetMaximum()/2.0);
		int bin3 = spectrumHistoW->FindFirstBinAbove(spectrumHistoW->GetMaximum()/5.0);
		int bin4 = spectrumHistoW->FindLastBinAbove(spectrumHistoW->GetMaximum()/5.0);
		std::stringstream ssCut20w;
		ssCut20w << "(ch" << channel << "/(" << SumChannels << ")) > " << spectrumHistoW->GetBinCenter(bin3) << " && " << "(ch" << channel << "/(" << SumChannels << ")) < "<<  spectrumHistoW->GetBinCenter(bin4);
		TCut w20percCut = ssCut20w.str().c_str();  //cut for w to get only the "relevant" part - TODO find a reasonable way to define this
		double meanW20 = (spectrumHistoW->GetBinCenter(bin4) + spectrumHistoW->GetBinCenter(bin3)) / 2.0;
		CurrentCrystal->SetW20percCut(w20percCut);
		double width20perc =spectrumHistoW->GetBinCenter(bin4) - spectrumHistoW->GetBinCenter(bin3);
		double fwhm = spectrumHistoW->GetBinCenter(bin2) - spectrumHistoW->GetBinCenter(bin1);
		double rms = spectrumHistoW->GetRMS();
		CurrentCrystal->SetHistoW(spectrumHistoW);
		CurrentCrystal->SetHistoWfwhm(fwhm);
		CurrentCrystal->SetHistoWrms(rms);
		CurrentCrystal->SetHistoWwidth20perc(width20perc);
		if(usingTaggingBench)
		{
		  sname << "gaussW histogram - Crystal " << CurrentCrystal->GetID();
		  TF1 *gaussW = new TF1(sname.str().c_str(),  "gaus",0,1);
		  int binmax = spectrumHistoW->GetMaximumBin();
		  double maximum = spectrumHistoW->GetXaxis()->GetBinCenter(binmax);
		  gaussW->SetParameter(0,maximum);
		  gaussW->SetParameter(1,spectrumHistoW->GetMean());
		  gaussW->SetParameter(2,rms);
		  spectrumHistoW->Fit(sname.str().c_str(),"QR");
		  doiFile << CurrentCrystal->GetX() << " " << CurrentCrystal->GetY() << " " << gaussW->GetParameter(1) << " " <<  gaussW->GetParameter(2) << std::endl;;
		  CurrentCrystal->SetHistoWfit(gaussW);
		}
		var.str("");
		sname.str("");
		
		
		
		//histogram of w versus adc channels
		//it will be useful fot doi correction
		long long int nPoints;
		sname << "ADC channels vs. W - Crystal " << CurrentCrystal->GetID();
		var << SumChannels << ":FloodZ >> " << sname.str();
		TH2F* spectrum2dADCversusW = new TH2F(sname.str().c_str(),sname.str().c_str(),100,0,1,histo1Dbins,0,histo1Dmax);
		nPoints = tree->Draw(var.str().c_str(),CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCut+w20percCut,"COLZ");
		spectrum2dADCversusW->GetXaxis()->SetTitle("W");
		spectrum2dADCversusW->GetYaxis()->SetTitle("ADC channels");
		double parM;
		var.str("");
		sname.str("");
		
		if(correctingForDOI)
		{
		  //mayhem
		  //i.e. use FitSlicesX to get the gaussian fit of each slices of the TH2F. Sliced in x (bin by bin)
		  spectrum2dADCversusW->FitSlicesY(0, 0, -1, 0, "QNR");
		  sname << spectrum2dADCversusW->GetName() << "_1";
		  TH1D *spectrum2d_1 = (TH1D*)gDirectory->Get(sname.str().c_str()); // _1 is the TH1D automatically created by ROOT when FitSlicesX is called, holding the TH1F of the mean values
		  
		  
		  sname << "linearCrystal - Crystal " << CurrentCrystal->GetID();
		  TF1 *linearCrystal = new TF1(sname.str().c_str(),  "[0]*x + [1]",0,1);
		  spectrum2d_1->Fit(sname.str().c_str(),"Q");
		  
		  parM = linearCrystal->GetParameter(0); // m parameter for the linear fit to correct energy res for DOI
		  // 		  double parQ = linearCrystal->GetParameter(1); // q parameter for the linear fit to correct energy res for DOI
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
		  baseVar << "(("  <<  SumChannels<< " ) - ( ( FloodZ - " <<  meanW20 << " ) * ( " << parM << ") ))";
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
		  float fitmin_corr = par1_corr-1.5*par2_corr;
		  float fitmax_corr = par1_corr+1.8*par2_corr;
		  
		  sname << "gauss_corr - Crystal " << CurrentCrystal->GetID() << " - MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
		  TF1 *gauss_corr = new TF1(sname.str().c_str(),  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitmin_corr,fitmax_corr);
		  gauss_corr->SetParameter(0,par0);
		  gauss_corr->SetParameter(1,par1);
		  gauss_corr->SetParameter(2,par2); 
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
		  streamEnergyCutCorrected << baseVar.str() << " > " << gauss_corr->GetParameter(1) - 2.5*std::abs(gauss_corr->GetParameter(2)) << " && " << baseVar.str() << " < " << gauss_corr->GetParameter(1) + 4.0*std::abs(gauss_corr->GetParameter(2));
		  PhotopeakEnergyCutCorrected = streamEnergyCutCorrected.str().c_str();
		  // 		CurrentCrystal->SetSpectrum(*spectrum);
		  sname.str("");
		  
		  //then prepare the highlighted spectrum and store it in the crystal
		  sname << "Hg Charge Spectrum Correctd - Crystal " << CurrentCrystal->GetID();
		  var << "("  <<  SumChannels<< " ) - ( ( FloodZ - " <<  meanW20 << " ) * ( " << parM << ") ) >> " << sname.str();
		  TH1F* spectrumChargeCorrectedHighlighted = new TH1F(sname.str().c_str(),sname.str().c_str(),histo1Dbins,1,histo1Dmax);	  
		  tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCutCorrected);
		  spectrumChargeCorrectedHighlighted->GetXaxis()->SetTitle("ADC Channels");
		  spectrumChargeCorrectedHighlighted->GetYaxis()->SetTitle("N");
		  CurrentCrystal->SetHighlightedSpectrumCorrected(spectrumChargeCorrectedHighlighted);
		  var.str("");
		  sname.str("");
		  
		}
		
		// a 3d historgram for this crystal, mainly to check the 3d cut
		sname << "Flood Histogram 3D - Crystal " << CurrentCrystal->GetID();
		var << "FloodZ:FloodY:FloodX >> " << sname.str();
		TH3I* spectrum3dCrystal = new TH3I(sname.str().c_str(),sname.str().c_str(),histo3DchannelBin,minX3Dplot,maxX3Dplot,histo3DchannelBin,minY3Dplot,maxY3Dplot,histo3DchannelBin,0,1);
		tree->Draw(var.str().c_str(),CutXYZ + CutTrigger + CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName());
		spectrum3dCrystal->GetXaxis()->SetTitle("U");
		spectrum3dCrystal->GetYaxis()->SetTitle("V");
		spectrum3dCrystal->GetZaxis()->SetTitle("W");
		CurrentCrystal->SetFloodMap3D(spectrum3dCrystal);
		sname.str("");
		var.str("");
		
		// a 3d historgram for this crystal, mainly to check the 3d cut
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
		sname << "Density Histogram - Crystal " << CurrentCrystal->GetID();
		Int_t u,v,w;  
		spectrum3dCrystal->GetMaximumBin(u,v,w); // get the maximum bin of the 3d histo
		 //get ax bin content
		TH1F* densityHisto = new TH1F(sname.str().c_str(),sname.str().c_str(),spectrum3dCrystal->GetBinContent(u,v,w)-1,1,spectrum3dCrystal->GetBinContent(u,v,w));
		densityHisto->GetXaxis()->SetTitle("Numb of Entries");
		int NbinX = spectrum3dCrystal->GetXaxis()->GetNbins();
                int NbinY = spectrum3dCrystal->GetYaxis()->GetNbins();
                int NbinZ = spectrum3dCrystal->GetZaxis()->GetNbins();
		for(int iContent = 1 ; iContent < NbinX+1 ; iContent++) 
		{
		  for(int jContent = 1 ; jContent < NbinY+1 ; jContent++) 
		  {
		    for(int kContent = 1 ; kContent < NbinZ+1 ; kContent++)
		    {
		      densityHisto->Fill(spectrum3dCrystal->GetBinContent(iContent,jContent,kContent));
		    }
		  }
		}
		CurrentCrystal->SetDensityHisto(densityHisto);
		sname.str("");
		
		// Histogram 2d of the photopeak time evolution
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
		nPoints = tree->Draw(var.str().c_str(),CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName(),"COLZ");
		spectrum2dADCversusWComplete->GetXaxis()->SetTitle("W");
		spectrum2dADCversusWComplete->GetYaxis()->SetTitle("ADC channels");
		spectrum2dADCversusWComplete->SetName(sname.str().c_str());
		CurrentCrystal->SetADCversusWComplete(spectrum2dADCversusWComplete);
		var.str("");
		sname.str("");
		
		//w histogram with energy cut on the corrected spectrum
		sname << "W histogram Corrected - Crystal " << CurrentCrystal->GetID();
		var << "(ch" << channel << "/(" << SumChannels << ")) >> " << sname.str();
		TH1F* spectrumHistoWCorrected = new TH1F(sname.str().c_str(),sname.str().c_str(),250,0,1);	  
		tree->Draw(var.str().c_str(),CutXYZ + CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCutCorrected+triggerPhotopeakCut);
		spectrumHistoWCorrected->GetXaxis()->SetTitle("W");
		spectrumHistoWCorrected->GetYaxis()->SetTitle("N");
		CurrentCrystal->SetHistoWCorrected(spectrumHistoWCorrected);
		var.str("");
		sname.str("");
		
		
		if(usingRealSimData) // only if this is a sim dataset
		{
		  // 		  long long int nPoints;
		  // a 2d plot of real vs. w, using the same cuts as before
		  if(correctingForDOI)
		    PhotopeakEnergyCut = PhotopeakEnergyCutCorrected;
		  
		  sname << "Real Z vs. W - Crystal " << CurrentCrystal->GetID();
		  var << "-(RealZ-" << CurrentCrystal->GetDimensionZ()/2.0 << "):FloodZ >> " << sname.str(); 
		  TH2F* spectrum2dSimDOIplot = new TH2F(sname.str().c_str(),sname.str().c_str(),100,0,1,100,0,CurrentCrystal->GetDimensionZ());
		  nPoints = tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CurrentCrystal->GetZXCut()->GetName() + CurrentCrystal->GetZYCut()->GetName()+PhotopeakEnergyCut,"COLZ"); 
		  spectrum2dSimDOIplot->GetXaxis()->SetTitle("W");
		  spectrum2dSimDOIplot->GetYaxis()->SetTitle("Z");
		  CurrentCrystal->SetSimDOIplot(spectrum2dSimDOIplot);
		  sname.str("");
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
		  sname << "expfit - Crystal " << CurrentCrystal->GetID();
		  TF1 *expfit = new TF1(sname.str().c_str(),  "[0]*exp(-x/[1])",0,1);
// 		  linear->SetParameter(0,-100);
// 		  linear->SetParameter(1,50);
		  expfit->SetParameter(0,50);
		  expfit->SetParameter(1,0.1);
		  // 		  graph->SetStats(1);
		  graph->Fit(sname.str().c_str(),"Q","",0.1,0.7);
		  CurrentCrystal->SetSimFit(expfit);
		  CurrentCrystal->SetSimGraph(graph);
		  sname.str("");
		  var.str("");
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
  
  //   TCanvas* SphericalCanvas = new TCanvas("Spherical","Spherical",1200,800);
  //   TCanvas* CylindricalXCanvas = new TCanvas("CylindricalX","CylindricalX",1200,800);
  //   TCanvas* CylindricalYCanvas = new TCanvas("CylindricalY","CylindricalY",1200,800);
  RawCanvas->Divide(4,4);
  TriggerCanvas->Divide(4,4);
  FloodHistoCanvas->Divide(4,4);
  FloodHisto3DCanvas->Divide(4,4);
  //   SphericalCanvas->Divide(4,4);
  //   CylindricalXCanvas->Divide(4,4);
  //   CylindricalYCanvas->Divide(4,4);
  //canvases for the global plots
  
  TCanvas* C_TaggingCrystalSpectrum= new TCanvas("TaggingCrystalSpectrum","TaggingCrystalSpectrum",1200,800);
  
  
  TCanvas* GlobalFlood2D = new TCanvas("Flood Histogram 2D","Flood Histogram 2D",800,800);
  TCanvas* GlobalFlood2DClean = new TCanvas("Flood Histogram 2D Clean","",800,800);
  TCanvas* GlobalFlood3D = new TCanvas("Flood Histogram 3D","Flood Histogram 3D",800,800);
  
  TCanvas* BigSpectraCanvas = new TCanvas("BigSpectra","",800,800);
  BigSpectraCanvas->Divide(ncrystalsx*nmppcx*nmodulex,ncrystalsy*nmppcy*nmoduley);
  int canvascounter = 1;
  for(int jCrystal = ncrystalsy*nmppcy*nmoduley -1 ; jCrystal >=0  ; jCrystal--)
  {
    for(int iCrystal = 0 ; iCrystal < ncrystalsx*nmppcx*nmodulex   ; iCrystal++)
    {
      BigSpectraCanvas->cd(canvascounter);
      if(crystal[iCrystal][jCrystal]->CrystalIsOn())
      {
        crystal[iCrystal][jCrystal]->GetCorrectedSpectrum()->SetFillStyle(3001);
        crystal[iCrystal][jCrystal]->GetCorrectedSpectrum()->SetFillColor(kBlue);
        crystal[iCrystal][jCrystal]->GetCorrectedSpectrum()->Draw();
      }
      canvascounter++;
      //std::cout << crystal[iCrystal][jCrystal]->GetID() << "\t";
    } 
    //std::cout << std::endl;
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
	  RawCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()); 
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetRawSpectrum()->Draw(); 
	  TriggerCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()); 
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrum()->Draw();
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
  TH1F *WfwhmDistro = new TH1F("w_fwhm","Distribution of FWHM in W plots",200,0,0.5);
  WfwhmDistro->GetXaxis()->SetTitle("W");
  WfwhmDistro->GetYaxis()->SetTitle("N");
  WfwhmDistro->SetStats(1);
  //2d histogram of fwhm
  TH2F *WfwhmVsIJ = new TH2F("w_fwhm vs. i,j","Distribution of FWHM in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  WfwhmVsIJ->GetXaxis()->SetTitle("i");
  WfwhmVsIJ->GetYaxis()->SetTitle("j");
  WfwhmVsIJ->GetZaxis()->SetTitle("w FHWM");
  //   WfwhmVsIJ->GetZaxis()->SetRangeUser(0,0.25);
  //   WfwhmVsIJ->SetStats(1);
  //histogram of rms
  TH1F *WrmsDistro = new TH1F("w_rms","Distribution of RMS in W plots",200,0,0.5);
  WrmsDistro->GetXaxis()->SetTitle("W");
  WrmsDistro->GetYaxis()->SetTitle("N");
  WrmsDistro->SetStats(1);
  //2d histogram of rms
  TH2F *WrmsVsIJ = new TH2F("w_rms vs. i,j","Distribution of RMS in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  WrmsVsIJ->GetXaxis()->SetTitle("i");
  WrmsVsIJ->GetYaxis()->SetTitle("j");
  WrmsVsIJ->GetZaxis()->SetTitle("w RMS");
  //Distribution of FWHM of W plots
  TH1F *Wwidth20perc = new TH1F("w20","Distribution of width at 20% in W plots",200,0,0.5);
  Wwidth20perc->GetXaxis()->SetTitle("W");
  Wwidth20perc->GetYaxis()->SetTitle("N");
  Wwidth20perc->SetStats(1);
  
  //same plot but just for the "not lateral" channels
  TH1F *Wwidth20percCentral = new TH1F("Central w20 - Central Crystals","Distribution of width at 20% in W plots - Central Crystals",200,0,0.5);
  Wwidth20percCentral->GetXaxis()->SetTitle("W");
  Wwidth20percCentral->GetYaxis()->SetTitle("N");
  Wwidth20percCentral->SetStats(1);
  
  TH1F *PeakPositionDistroCentral = new TH1F("Central LY","Distribution photopeak positions - Central Crystals",200,0,histo1Dmax);
  PeakPositionDistroCentral->GetXaxis()->SetTitle("ADC Channels");
  PeakPositionDistroCentral->GetYaxis()->SetTitle("N");
  PeakPositionDistroCentral->SetStats(1);
  
  TH1F *PeakEnergyResolutionDistroCentral = new TH1F("Central En Res","Distribution photopeak energy resolutions FWHM - Central Crystals",200,0,1);
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
  PeakPositionVsIJ->GetZaxis()->SetRangeUser(0,12000);
  //2d histogram
  TH2F *EnergyResolutionVsIJ = new TH2F("Energy res FWHM vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  EnergyResolutionVsIJ->GetXaxis()->SetTitle("i (U axis)");
  EnergyResolutionVsIJ->GetYaxis()->SetTitle("j (V axis)");
  EnergyResolutionVsIJ->GetZaxis()->SetTitle("Energy Resolution FWHM");
  EnergyResolutionVsIJ->GetXaxis()->SetTitleOffset(1.8);
  EnergyResolutionVsIJ->GetYaxis()->SetTitleOffset(1.8);
  EnergyResolutionVsIJ->GetZaxis()->SetTitleOffset(2.2);
  EnergyResolutionVsIJ->GetZaxis()->SetRangeUser(0,0.3);
  
  TH2F *EnergyResolutionVsIJ_corr = new TH2F("Corrected Energy res FWHM vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  EnergyResolutionVsIJ_corr->GetXaxis()->SetTitle("i (U axis)");
  EnergyResolutionVsIJ_corr->GetYaxis()->SetTitle("j (V axis)");
  EnergyResolutionVsIJ_corr->GetZaxis()->SetTitle("Energy Resolution FWHM");
  EnergyResolutionVsIJ_corr->GetXaxis()->SetTitleOffset(1.8);
  EnergyResolutionVsIJ_corr->GetYaxis()->SetTitleOffset(1.8);
  EnergyResolutionVsIJ_corr->GetZaxis()->SetTitleOffset(2.2);
  EnergyResolutionVsIJ_corr->GetZaxis()->SetRangeUser(0,0.3);
  
  
  TH2F *Wwidht20percVsIJ = new TH2F("w20 vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  Wwidht20percVsIJ->GetXaxis()->SetTitle("i (U axis)");
  Wwidht20percVsIJ->GetYaxis()->SetTitle("i (V axis)");
  Wwidht20percVsIJ->GetZaxis()->SetRangeUser(0,0.25);
  Wwidht20percVsIJ->GetXaxis()->SetTitleOffset(1.8);
  Wwidht20percVsIJ->GetYaxis()->SetTitleOffset(1.8);
  Wwidht20percVsIJ->GetZaxis()->SetTitleOffset(2.2);
  Wwidht20percVsIJ->GetZaxis()->SetTitle("W width at 20%");
  
  
  //Distribution of DOI resolutions - not very nice since one parameter in the calculation is assumed (from the DOI bench results)
  TH1F *WDoiDistro = new TH1F("Doi Res","Distribution of doi res",20,0,6);
  WDoiDistro->GetXaxis()->SetTitle("doi");
  WDoiDistro->GetYaxis()->SetTitle("N");
  WDoiDistro->SetStats(1);
  
  //Distribution of fit exp tau for W plots
  TH1F *WtauFit = new TH1F("Exp slope W","Distribution of exp slopes in W plots",1000,0,1);
  WtauFit->GetXaxis()->SetTitle("Tau");
  WtauFit->GetYaxis()->SetTitle("N");
  WtauFit->SetStats(1);
  
  
  TH2F *WtauFitVsIJ = new TH2F("Exp slopes W vs. i,j","Distribution of exp slopes in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  WtauFitVsIJ->GetXaxis()->SetTitle("i");
  WtauFitVsIJ->GetYaxis()->SetTitle("i");
  WtauFitVsIJ->GetZaxis()->SetTitle("tau");
  
  
  
  
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
      }

      RawCanvas->Write();
      TriggerCanvas->Write();
      BigSpectraCanvas->Write();
      
//       C_global = new TCanvas("C_global","C_global",800,800);
//       name = "3D Cuts - " + ModuleDirStream.str();
//       C_global->SetName(name);
//       C_global->SetTitle(name);
//       int colorcounter = 0;
      
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
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
	      CrystalDirStream << "Crystal " <<  CurrentCrystal->GetID();
	      directory[iModule+jModule][(iMppc+jMppc)+1][(iCry+jCry)+1] = directory[iModule+jModule][(iMppc+jMppc)+1][0]->mkdir(CrystalDirStream.str().c_str());
	      directory[iModule+jModule][(iMppc+jMppc)+1][(iCry+jCry)+1]->cd(); 
	      if(CurrentCrystal->CrystalIsOn() /*| usingRealSimData*/) // save data only if the crystal was found
	      {
		//create a pointer for the current crystal (mainly to make the code more readable)
		Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
		//fill the global distributions histograms
		PeakPositionDistro->Fill(CurrentCrystal->GetPhotopeakPosition());
		PeakEnergyResolutionDistro->Fill(CurrentCrystal->GetPhotopeakEnergyResolution());
		
		if(correctingForDOI)
		{
		  PeakEnergyResolutionDistro_corr->Fill(CurrentCrystal->GetPhotopeakEnergyResolutionCorrected());
		}
		WfwhmDistro->Fill(CurrentCrystal->GetWfwhm());
		WDoiDistro->Fill( (15.0/CurrentCrystal->GetWfwhm())*0.0158); // FIXME CAREFUL: here the 0.0158 value is hardcoded and taken from the sigma of W distros in DOI bench setup. 15.0 is the length of the crystals in mm.
		PeakPositionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakPosition());
		EnergyResolutionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakEnergyResolution());
		EnergyResolutionVsIJ_corr->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakEnergyResolutionCorrected());
		WfwhmVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetWfwhm());
		WrmsDistro->Fill(CurrentCrystal->GetWrms());
		WrmsVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetWrms());
		Wwidht20percVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetWwidth20perc());
		Wwidth20perc->Fill(CurrentCrystal->GetWwidth20perc());
		
		
		if( ((iModule*nmppcx)+iMppc) > 0 && (((iModule*nmppcx)+iMppc) < nmppcx -1) && ((jModule*nmppcy)+jMppc) > 0 && (((jModule*nmppcy)+jMppc) < nmppcy -1 ))
		{
		  Wwidth20percCentral->Fill(CurrentCrystal->GetWwidth20perc());
		  PeakPositionDistroCentral->Fill(CurrentCrystal->GetPhotopeakPosition());
		  PeakEnergyResolutionDistroCentral->Fill(CurrentCrystal->GetPhotopeakEnergyResolution());
		  if(correctingForDOI)
		  {
		    PeakEnergyResolutionDistroCentral_corr->Fill(CurrentCrystal->GetPhotopeakEnergyResolutionCorrected());
		  }
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
		CurrentCrystal->GetHighlightedSpectrum()->SetFillColor(3);
		CurrentCrystal->GetHighlightedSpectrum()->Draw("same");
		CurrentCrystal->GetFit()->Draw("same");
		C_spectrum->Write();
		delete C_spectrum;
		
		// spectrum without highligth
		C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
		TString title = "mod_";
		title += CurrentCrystal->GetSpectrum()->GetName() ;
		C_spectrum->SetName(title);
		C_spectrum->SetTitle("");
		C_spectrum->cd();
		CurrentCrystal->GetSpectrum()->SetFillStyle(3001);
		CurrentCrystal->GetSpectrum()->SetFillColor(kBlue);
		CurrentCrystal->GetSpectrum()->Draw();
		C_spectrum->Write();
		delete C_spectrum;
		
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
		C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
		C_spectrum->SetName(CurrentCrystal->GetADCversusW()->GetName());
		C_spectrum->cd();
		CurrentCrystal->GetADCversusW()->Draw("COLZ");
		C_spectrum->Write();
		delete C_spectrum;
		
		// adc versus w complete
		C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
		C_spectrum->SetName(CurrentCrystal->GetADCversusWComplete()->GetName());
		C_spectrum->cd();
		CurrentCrystal->GetADCversusWComplete()->Draw("COLZ");
		C_spectrum->Write();
		delete C_spectrum;
		
		
		if(correctingForDOI)
		{
		  C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
		  C_spectrum->SetName(CurrentCrystal->GetSlicesMean()->GetName());
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
		
		//w histo corrected
		C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
		C_spectrum->SetName(CurrentCrystal->GetHistoWCorrected()->GetName());
		C_spectrum->cd();
		CurrentCrystal->GetHistoWCorrected()->Draw();
		C_spectrum->Write();
		delete C_spectrum;
		
		//density histo
		C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
		C_spectrum->SetName(CurrentCrystal->GetDensityHisto()->GetName());
		C_spectrum->cd();
		CurrentCrystal->GetDensityHisto()->Draw();
		C_spectrum->Write();
		delete C_spectrum;
		
		
		if(usingRealSimData)
		{
		  WtauFit->Fill(CurrentCrystal->GetSimFit()->GetParameter(1));
		  WtauFitVsIJ->Fill(CurrentCrystal->GetI() ,CurrentCrystal->GetJ() , CurrentCrystal->GetSimFit()->GetParameter(1));
		  
		  C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
		  C_spectrum->SetName(CurrentCrystal->GetSimDOIplot()->GetName());
		  C_spectrum->cd();
		  CurrentCrystal->GetSimDOIplot()->Draw("COLZ");
		  C_spectrum->Write();
		  delete C_spectrum;
		  
		  C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
		  C_spectrum->SetName(CurrentCrystal->GetSimGraph()->GetName());
		  C_spectrum->cd();
		  CurrentCrystal->GetSimGraph()->Draw("ap");
		  CurrentCrystal->GetSimFit()->Draw("same");
		  C_spectrum->Write();		
		  delete C_spectrum;
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
      
      
      directory[iModule+jModule][0][0]->cd(); // go back to main directory
//       C_global->Write();
//       delete C_global;
      //write the summary histos and canvases, that were filled during file saving
      PeakPositionDistro->Write();
      PeakEnergyResolutionDistro->Write();
      if(correctingForDOI)
      {
	PeakEnergyResolutionDistro_corr->Write();
      }
      WfwhmDistro->Write();
      WDoiDistro->Write();
      WrmsDistro->Write();
      Wwidth20perc->Write();
      Wwidth20percCentral->Write();
      PeakPositionDistroCentral->Write();
      PeakEnergyResolutionDistroCentral->Write();
      if(correctingForDOI)
      {
	PeakEnergyResolutionDistroCentral_corr->Write();
      }
      
      
      TCanvas *C_WfwhmVsIJ = new TCanvas("C_WfwhmVsIJ","C_WfwhmVsIJ",800,800);
      C_WfwhmVsIJ->SetName(WfwhmVsIJ->GetName());
      C_WfwhmVsIJ->cd();
      WfwhmVsIJ->Draw("LEGO2");
      C_WfwhmVsIJ->Write();
      
      TCanvas *C_WrmsVsIJ = new TCanvas("C_WrmsVsIJ","C_WrmsVsIJ",800,800);
      C_WrmsVsIJ->SetName(WrmsVsIJ->GetName());
      C_WrmsVsIJ->cd();
      WrmsVsIJ->Draw("LEGO2");
      C_WrmsVsIJ->Write();
      
      TCanvas *C_PeakPositionVsIJ = new TCanvas("C_PeakPositionVsIJ","C_PeakPositionVsIJ",800,800);
      C_PeakPositionVsIJ->SetName(PeakPositionVsIJ->GetName());
      C_PeakPositionVsIJ->cd();
      PeakPositionVsIJ->Draw("LEGO2");
      C_PeakPositionVsIJ->SetLeftMargin(0.15);
      C_PeakPositionVsIJ->Write();
      
      TCanvas *C_EnergyResolutionVsIJ = new TCanvas("C_EnergyResolutionVsIJ","C_EnergyResolutionVsIJ",800,800);
      C_EnergyResolutionVsIJ->SetName(EnergyResolutionVsIJ->GetName());
      C_EnergyResolutionVsIJ->cd();
      EnergyResolutionVsIJ->Draw("LEGO2");
      C_EnergyResolutionVsIJ->SetLeftMargin(0.15);
      C_EnergyResolutionVsIJ->Write();
      
      if(correctingForDOI)
      {
	TCanvas *C_EnergyResolutionVsIJ_corr = new TCanvas("C_EnergyResolutionVsIJ_corr","C_EnergyResolutionVsIJ_corr",800,800);
	C_EnergyResolutionVsIJ_corr->SetName(EnergyResolutionVsIJ_corr->GetName());
	C_EnergyResolutionVsIJ_corr->cd();
	EnergyResolutionVsIJ_corr->Draw("LEGO2");
	C_EnergyResolutionVsIJ_corr->SetLeftMargin(0.15);
	C_EnergyResolutionVsIJ_corr->Write();
      }
      
      TCanvas *C_Wwidht20percVsIJ = new TCanvas("C_Wwidht20percVsIJ","C_Wwidht20percVsIJ",800,800);
      C_Wwidht20percVsIJ->SetName(Wwidht20percVsIJ->GetName());
      C_Wwidht20percVsIJ->cd();
      Wwidht20percVsIJ->Draw("LEGO2");
      C_Wwidht20percVsIJ->SetLeftMargin(0.15);
      C_Wwidht20percVsIJ->Write();
      
      //       gStyle->SetOptStat(1);
      
      if(usingRealSimData)
      {
	WtauFit->Write();
	TCanvas *C_WtauFitVsIJ = new TCanvas("C_WtauFitVsIJ","C_WtauFitVsIJ",800,800);
	C_WtauFitVsIJ->SetName(WtauFitVsIJ->GetName());
	C_WtauFitVsIJ->cd();
	WtauFitVsIJ->Draw("LEGO2");
	C_WtauFitVsIJ->Write();
      }
      
      
      
    }
  }
  if(saveAnalysisTree) // save the TTree created for the analysis, if the user requires it in the config file
  {
    std::cout << "Saving analysis TTree to a file temp.root" << std::endl;
    TFile* fFile = new TFile("temp.root","recreate");
    fFile->cd();
    tree->Write();
    fFile->Close();
  }
  fPlots->Close();
  //----------------------------------------------------------//
  if(usingTaggingBench)
  {
    doiFile.close();
  }
  
  delete crystal;
  delete mppc;
  delete module;
  
  return 0;
}
