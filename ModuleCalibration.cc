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
#include "TH3F.h"
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
// #include "Classes.h"

#define ENERGY_RESOLUTION 0.12
#define ENERGY_RESOLUTION_SATURATION_CORRECTION 0.25

// void ReverseXAxis (TH1 *h)
// {
//    // Remove the current axis
//    h->GetXaxis()->SetLabelOffset(999);
//    h->GetXaxis()->SetTickLength(0);
//                                                                                 
//    // Redraw the new axis
//    gPad->Update();
//    TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
//                                 gPad->GetUymin(),
//                                 gPad->GetUxmin(),
//                                 gPad->GetUymin(),
//                                 h->GetXaxis()->GetXmin(),
//                                 h->GetXaxis()->GetXmax(),
//                                 510,"-");
//    newaxis->SetLabelOffset(-0.03);
//    newaxis->Draw();
// }
//                                                                                 
// void ReverseYAxis (TH1 *h)
// {
//    // Remove the current axis
//    h->GetYaxis()->SetLabelOffset(999);
//    h->GetYaxis()->SetTickLength(0);
//                                                                                 
//    // Redraw the new axis
//    gPad->Update();
//    TGaxis *newaxis = new TGaxis(gPad->GetUxmin(),
//                                 gPad->GetUymax(),
//                                 gPad->GetUxmin()-0.001,
//                                 gPad->GetUymin(),
//                                 h->GetYaxis()->GetXmin(),
//                                 h->GetYaxis()->GetXmax(),
//                                 510,"+");
//    newaxis->SetLabelOffset(-0.03);
//    newaxis->Draw();
// }



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
  std::cout<<"#                                                          #"<<std::endl;
  std::cout<<"#           New Clear PEM module calibration               #"<<std::endl;  
  std::cout<<"#                                                          #"<<std::endl;  
  std::cout<<"###########################################################"<<std::endl;
  std::cout<<"\n\n"<<std::endl;
  std::cout<<"=====>   C O N F I G U R A T I O N   <====\n"<<std::endl;
  
  
  //----------------------------------------------------------//
  //  Import input and parse the config file                  //
  //----------------------------------------------------------//
  // Set a default config file name
  std::string ConfigFileName = "config.cfg"; 
  // and assume as default that there is no config file name from command line 
  //then check
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
  input.CreateTree();                // create the TTree that will be used in analysis
  
  int ncrystalsx            = config.read<int>("ncrystalsx",2);             // number of crystals in x direction per mppc - default to 2 if the key is not found in the config file
  int ncrystalsy            = config.read<int>("ncrystalsy",2);             // number of crystals in y direction per mppc - default to 2 if the key is not found in the config file
  int nmppcx                = config.read<int>("nmppcx",2);                 // number of mppc in x direction per mppc - default to 2 if the key is not found in the config file
  int nmppcy                = config.read<int>("nmppcy",2);                 // number of mppc in y direction per mppc - default to 2 if the key is not found in the config file
  int nmodulex              = config.read<int>("nmodulex",1);               // number of modules in x direction per mppc - default to 1 if the key is not found in the config file
  int nmoduley              = config.read<int>("nmoduley",1);               // number of modules in y direction per mppc - default to 1 if the key is not found in the config file
  int histo1Dmax            = config.read<int>("histo1Dmax");               // max of the 1D charge histograms (in ADC channels)
  int histo1Dbins           = config.read<int>("histo1Dbins");              // number of bins of the 1D charge histograms
  int histo2DchannelBin     = config.read<int>("histo2DchannelBin");        // number of bins of the 2D flood histograms, for single channels
  int histo2DglobalBins     = config.read<int>("histo2DglobalBins");        // number of bins of the 2D flood histograms, for entire module
  int histo3DchannelBin     = config.read<int>("histo3DchannelBin");        // number of bins of the 3D flood histograms, for single channels
  int histo3DglobalBins     = config.read<int>("histo3DglobalBins");        // number of bins of the 3D flood histograms, for entire module
  bool saveAnalysisTree     = config.read<bool>("saveAnalysisTree");        // choice to save or not the analysis TTree, in a file temp.root
  float taggingPosition     = config.read<float>("taggingPosition");        // position of the tagging bench in mm 
  bool usingTaggingBench    = config.read<bool>("usingTaggingBench");       // true if the input is using tagging bench, false if not
  int taggingCrystalChannel = config.read<int>("taggingCrystalChannel");    // input channel where the tagging crystal information is stored
  bool correctingSaturation = config.read<bool>("correctingSaturation");;   // true if saturation correction is applied, false if it's not
  float energyResolution    = config.read<float>("expectedEnergyResolution",0); // energy resolution input by the user, if any, otherwise 0
  bool usingRealSimData     = config.read<bool>("usingRealSimData",0);
  
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
    sSumChannels << "+ch" << i; // the analysis ttree will always have channels in order, from 0 to input size
  SumChannels = sSumChannels.str();
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
  TH1F* spectrum;
  TH2F* spectrum2d;
  TH3F* spectrum3d;
  TGraph* simGraph; 
  TCanvas* C_spectrum;
  std::stringstream var,cut,sname; 
  
  TString name;
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      TCut CutXYZ = "FloodX > -7 && FloodX < 7 && FloodY > -7 && FloodY < 7 && FloodZ > 0 && FloodZ < 1";
      std::cout << "Generating global spectra..." << std::endl;
      
      // GLOBAL SPECTRA
      // Flood histogram
      spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DglobalBins,-7,7,histo2DglobalBins,-7,7);
      tree->Draw("FloodY:FloodX >> spectrum2d","","COLZ");
      name = "Flood Histogram 2D - " + module[iModule][jModule]->GetName();
      spectrum2d->SetName(name); 
      spectrum2d->SetTitle(name);
      spectrum2d->GetXaxis()->SetTitle("U");
      spectrum2d->GetYaxis()->SetTitle("V");
      module[iModule][jModule]->SetFloodMap2D(*spectrum2d);
      delete spectrum2d;
      //sherical coordinates plot
      spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DglobalBins,1.1,3.1415/2.0,histo2DglobalBins,-3.14/2.0,3.14/2.0); 
      tree->Draw("Phi:Theta >> spectrum2d",CutXYZ,"COLZ");
      name = "Spherical Plot - " + module[iModule][jModule]->GetName();
      spectrum2d->SetName(name);
      spectrum2d->SetTitle(name);
      spectrum2d->GetXaxis()->SetTitle("Theta");
      spectrum2d->GetYaxis()->SetTitle("Phi");
      module[iModule][jModule]->SetSphericalMap(*spectrum2d);
      delete spectrum2d;
      //cylindrical coordinates plot, x and theta
      spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DglobalBins,-7,7,histo2DglobalBins,1.1,3.1415/2.0); 
      tree->Draw("Theta:FloodX >> spectrum2d",CutXYZ,"COLZ");
      name = "Cylindrical Plot Theta:X - Module " + module[iModule][jModule]->GetName();
      spectrum2d->SetName(name);
      spectrum2d->SetTitle(name);
      spectrum2d->GetXaxis()->SetTitle("U");
      spectrum2d->GetYaxis()->SetTitle("Theta");
      module[iModule][jModule]->SetCylindricalXMap(*spectrum2d);
      delete spectrum2d;
      //cylindrical coordinates plot, y and theta
      spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DglobalBins,-7,7,histo2DglobalBins,1.1,3.1415/2.0); 
      tree->Draw("Theta:FloodY >> spectrum2d",CutXYZ,"COLZ");
      name = "Cylindrical Plot Theta:Y - Module " + module[iModule][jModule]->GetName();
      spectrum2d->SetName(name);
      spectrum2d->SetTitle(name);
      spectrum2d->GetXaxis()->SetTitle("V");
      spectrum2d->GetYaxis()->SetTitle("Theta");
      module[iModule][jModule]->SetCylindricalYMap(*spectrum2d);
      delete spectrum2d;
      //3D plot
      spectrum3d = new TH3F("spectrum3d","spectrum3d",histo3DglobalBins,-7,7,histo3DglobalBins,-7,7,histo3DglobalBins,0,1);
      tree->Draw("FloodZ:FloodY:FloodX >> spectrum3d",CutXYZ);
      name = "Flood Histogram 3D - Module " + module[iModule][jModule]->GetName();
      spectrum3d->SetName(name);
      spectrum3d->SetTitle(name);
      spectrum3d->GetXaxis()->SetTitle("U");
      spectrum3d->GetYaxis()->SetTitle("V");
      spectrum3d->GetZaxis()->SetTitle("W");
      module[iModule][jModule]->SetFloodMap3D(*spectrum3d);
      delete spectrum3d;
      //spectra for each mppc
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  //FIXME careful, at the moment it works only because there's one module
	  // it should be fixed for more modules by using the same mppc[][] logic used for the crystals, below
	  int channel = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
	  std::cout << "Generating spectra for MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " ..." << std::endl;
	  cut << "TriggerChannel == " << channel  ;
	  TCut CutTrigger = cut.str().c_str();
	  cut.str("");
	  //same as the global ones, but selecting on TriggerChannel
	  // raw spectrum
	  spectrum = new TH1F("spectrum","spectrum",histo1Dbins,1,histo1Dmax);
	  channel = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
	  var << "ch" << channel << " >> spectrum";
	  tree->Draw(var.str().c_str(),"");
	  name = "Raw Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() + " - Module " + module[iModule][jModule]->GetName();
	  spectrum->SetName(name);
	  spectrum->SetTitle(name);
	  spectrum->GetXaxis()->SetTitle("ADC Channels");
	  spectrum->GetYaxis()->SetTitle("N");
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetRawSpectrum(*spectrum);
	  var.str("");
	  delete spectrum;
	  
	  //trigger selected spectrum
	  spectrum = new TH1F("spectrum","spectrum",histo1Dbins,1,histo1Dmax);	  
	  var << "ch" << channel << " >> spectrum";
	  tree->Draw(var.str().c_str(),CutTrigger);
	  name = "Trigger Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
	  spectrum->SetName(name);
	  spectrum->SetTitle(name);
	  spectrum->GetXaxis()->SetTitle("ADC Channels");
	  spectrum->GetYaxis()->SetTitle("N");
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetTriggerSpectrum(*spectrum);
	  var.str("");
	  delete spectrum;
	  // Flood histogram
	  spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DchannelBin,-7,7,histo2DchannelBin,-7,7);
	  tree->Draw("FloodY:FloodX >> spectrum2d",CutTrigger,"COLZ");
	  name = "Flood Histogram 2D - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
	  spectrum2d->SetName(name); 
	  spectrum2d->SetTitle(name);
	  spectrum2d->GetXaxis()->SetTitle("U");
	  spectrum2d->GetYaxis()->SetTitle("V");
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetFloodMap2D(*spectrum2d);
	  delete spectrum2d;
	  //sherical coordinates plot
	  spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DchannelBin,1.1,3.1415/2.0,histo2DchannelBin,-3.14/2.0,3.14/2.0); 
	  tree->Draw("Phi:Theta >> spectrum2d",CutXYZ+CutTrigger,"COLZ");
	  name = "Spherical Plot - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
	  spectrum2d->SetName(name);
	  spectrum2d->SetTitle(name);
	  spectrum2d->GetXaxis()->SetTitle("Theta");
	  spectrum2d->GetYaxis()->SetTitle("Phi");
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetSphericalMap(*spectrum2d);
	  delete spectrum2d;
	  //cylindrical coordinates plot, x and theta
	  spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DchannelBin,-7,7,histo2DchannelBin,1.1,3.1415/2.0); 
	  tree->Draw("Theta:FloodX >> spectrum2d",CutXYZ+CutTrigger,"COLZ");
	  name = "Cylindrical Plot Theta:X - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
	  spectrum2d->SetName(name);
	  spectrum2d->SetTitle(name);
	  spectrum2d->GetXaxis()->SetTitle("U");
	  spectrum2d->GetYaxis()->SetTitle("Theta");
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetCylindricalXMap(*spectrum2d);
	  delete spectrum2d;
	  //cylindrical coordinates plot, y and theta
	  spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DchannelBin,-7,7,histo2DchannelBin,1.1,3.1415/2.0); 
	  tree->Draw("Theta:FloodY >> spectrum2d",CutXYZ+CutTrigger,"COLZ");
	  name = "Cylindrical Plot Theta:Y - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
	  spectrum2d->SetName(name);
	  spectrum2d->SetTitle(name);
	  spectrum2d->GetXaxis()->SetTitle("V");
	  spectrum2d->GetYaxis()->SetTitle("Theta");
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetCylindricalYMap(*spectrum2d);
	  delete spectrum2d;
	  spectrum3d = new TH3F("spectrum3d","spectrum3d",histo3DchannelBin,-7,7,histo3DchannelBin,-7,7,histo3DchannelBin,0,1);
	  tree->Draw("FloodZ:FloodY:FloodX >> spectrum3d",CutXYZ+CutTrigger);
	  name = "Flood Histogram 3D - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
	  spectrum3d->SetName(name);
	  spectrum3d->SetTitle(name);
	  spectrum3d->GetXaxis()->SetTitle("U");
	  spectrum3d->GetYaxis()->SetTitle("V");
	  spectrum3d->GetZaxis()->SetTitle("W");
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetFloodMap3D(*spectrum3d);
	  delete spectrum3d;
	  
	  //spectra for each crystal
	  for(int iCry = 0; iCry < ncrystalsx ; iCry++)
	  {
	    for(int jCry = 0; jCry < ncrystalsy ; jCry++)
	    {
	      
	      Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
	      if(CurrentCrystal->CrystalIsOn() /*| usingRealSimData*/)
	      {
		std::cout << "Generating spectra for crystal " << CurrentCrystal->GetID() << " ..." << std::endl;
		
		
		TCut CutCrystal;
// 		if(usingRealSimData)
// 		{
// 		  std::stringstream simCutCrystal;
// 		  simCutCrystal << "RealX > " << CurrentCrystal->GetX() - CurrentCrystal->GetDimensionX()/2.0 << " && RealX < " << CurrentCrystal->GetX() + CurrentCrystal->GetDimensionX()/2.0 << "  && RealY > " << CurrentCrystal->GetY() - CurrentCrystal->GetDimensionY()/2.0 << " && RealY < " << CurrentCrystal->GetY() + CurrentCrystal->GetDimensionY()/2.0;
// 		  CutCrystal = simCutCrystal.str().c_str();
// // 		  std::cout << CurrentCrystal->GetID() << " " << CutCrystal << std::endl;
// 		}
// 		else
// 		{
		CutCrystal = CurrentCrystal->GetCrystalCut();
// 		}
		
		//-------------------------------------------------------------------------
		//standard sum spectrum with cut on crystal events, xyz and trigger channel
		//------------------------------------------------------------------------
		//draw spectrum
		spectrum = new TH1F("spectrum","spectrum",histo1Dbins,1,histo1Dmax);	  
		var << SumChannels << " >> spectrum";
		tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CutCrystal);
		sname << "Charge Spectrum - Crystal " << CurrentCrystal->GetID() << " - " << CurrentCrystal->GetExtendedID();
		spectrum->SetName(sname.str().c_str());
		spectrum->SetTitle(sname.str().c_str());
		spectrum->GetXaxis()->SetTitle("ADC Channels");
		spectrum->GetYaxis()->SetTitle("N");
// 		std::cout << var.str().c_str() << std::endl;
// 		std::cout << CutXYZ << std::endl;
// 		std::cout << CutTrigger << std::endl;
// 		std::cout << CutCrystal << std::endl;
		//automatically look for the 511Kev peak to find the photopeak energy cut
		//find peaks in each crystal spectrum, with TSpectrum
		TSpectrum *s;
		s = new TSpectrum(20);
// 		Input[i].SumSpectraCanvas->cd(j+1);
		Int_t CrystalPeaksN = s->Search(spectrum,2,"goff",0.5); 
		Float_t *CrystalPeaks = s->GetPositionX();
		Float_t *CrystalPeaksY = s->GetPositionY();
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
		float par0 = CrystalPeaksY[peakID];
		float par1 = CrystalPeaks[peakID];
		float par2 = (CrystalPeaks[peakID]*energyResolution)/2.35;
		float fitmin = par1-1*par2;
		float fitmax = par1+1.8*par2;
		TF1 *gauss = new TF1("gauss",  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitmin,fitmax);
		gauss->SetParameter(0,par0);
		gauss->SetParameter(1,par1);
		gauss->SetParameter(2,par2); 
		spectrum->Fit("gauss","Q","",fitmin,fitmax);
		//store the mean and sigma in the crystal
		if(gauss->GetParameter(1) > 0) // otherwise the fit was very wrong..)
		  CurrentCrystal->SetPhotopeak(gauss->GetParameter(1),std::abs(gauss->GetParameter(2)));
		CurrentCrystal->SetFit(*gauss);
// 		std::cout << "Photopeak Mean for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakPosition() << std::endl;
// 		std::cout << "Photopeak Sigma for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakSigma() << std::endl;
// 		std::cout << "Photopeak Energy Resolution FWHM for crystal " << CurrentCrystal->GetID() << " = " << CurrentCrystal->GetPhotopeakEnergyResolution() << std::endl;
		//Compute the energy Tcut
		std::stringstream streamEnergyCut;
		streamEnergyCut << SumChannels << " > " << gauss->GetParameter(1) - 1.5*std::abs(gauss->GetParameter(2)) << " && " << SumChannels << " < " << gauss->GetParameter(1) + 3.0*std::abs(gauss->GetParameter(2));
		TCut PhotopeakEnergyCut = streamEnergyCut.str().c_str();
		CurrentCrystal->SetSpectrum(*spectrum);
		var.str("");
		sname.str("");
		delete gauss;
		delete spectrum;
		// then prepare the highlighted spectrum and store it in the crystal
		spectrum = new TH1F("spectrum","spectrum",histo1Dbins,1,histo1Dmax);	  
		var << SumChannels << " >> spectrum";
		tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CutCrystal+PhotopeakEnergyCut);
		sname << "Hg Charge Spectrum - Crystal " << CurrentCrystal->GetID();
		spectrum->SetName(sname.str().c_str());
		spectrum->SetTitle(sname.str().c_str());
		spectrum->GetXaxis()->SetTitle("ADC Channels");
		spectrum->GetYaxis()->SetTitle("N");
		CurrentCrystal->SetHighlightedSpectrum(*spectrum);
		var.str("");
		sname.str("");
		delete spectrum;
		//-----------------------------------------------------------------------
		
		//w histogram with cut on crystal events, xyz and trigger channel and cut on photopeak
		spectrum = new TH1F("spectrum","spectrum",250,0,1);	  
		var << "(ch" << channel << "/(" << SumChannels << ")) >> spectrum";
		tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CutCrystal+PhotopeakEnergyCut);
		sname << "W histogram - Crystal " << CurrentCrystal->GetID();
		spectrum->SetName(sname.str().c_str());
		spectrum->SetTitle(sname.str().c_str());
		spectrum->GetXaxis()->SetTitle("W");
		spectrum->GetYaxis()->SetTitle("N");
		int bin1 = spectrum->FindFirstBinAbove(spectrum->GetMaximum()/2.0);
		int bin2 = spectrum->FindLastBinAbove(spectrum->GetMaximum()/2.0);
		int bin3 = spectrum->FindFirstBinAbove(spectrum->GetMaximum()/5.0);
		int bin4 = spectrum->FindLastBinAbove(spectrum->GetMaximum()/5.0);
		double width20perc =spectrum->GetBinCenter(bin4) - spectrum->GetBinCenter(bin3);
		double fwhm = spectrum->GetBinCenter(bin2) - spectrum->GetBinCenter(bin1);
		double rms = spectrum->GetRMS();
		CurrentCrystal->SetHistoW(*spectrum);
		CurrentCrystal->SetHistoWfwhm(fwhm);
		CurrentCrystal->SetHistoWrms(rms);
		CurrentCrystal->SetHistoWwidth20perc(width20perc);
		var.str("");
		sname.str("");
		delete spectrum;
		
		// Flood histogram 2d for this crystal, to show the elliptic cut
		spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DglobalBins,-7,7,histo2DglobalBins,-7,7);
		tree->Draw("FloodY:FloodX >> spectrum2d",CutTrigger+CutCrystal,"COLZ");
		sname << "Flood Histogram 2D - Crystal " << CurrentCrystal->GetID();
		spectrum2d->SetName(sname.str().c_str()); 
		spectrum2d->SetTitle(sname.str().c_str());
		spectrum2d->GetXaxis()->SetTitle("U");
		spectrum2d->GetYaxis()->SetTitle("V");
		CurrentCrystal->SetFloodMap2D(*spectrum2d);
		sname.str("");
		delete spectrum2d;
		
		// Histogram 2d of the photopeak time evolution
		spectrum2d = new TH2F("spectrum2d","spectrum2d",250,0,tree->GetMaximum("ExtendedTimeTag"),histo1Dbins,0,histo1Dmax);
		var << SumChannels << ":ExtendedTimeTag >> spectrum2d";
		tree->Draw(var.str().c_str(),CutTrigger+CutCrystal,"COLZ");
		sname << "ADC channels vs. Time - Crystal " << CurrentCrystal->GetID();
		spectrum2d->SetName(sname.str().c_str()); 
		spectrum2d->SetTitle(sname.str().c_str());
		spectrum2d->GetXaxis()->SetTitle("ExtendedTimeTag");
		spectrum2d->GetYaxis()->SetTitle("ADC channels");
		CurrentCrystal->SetVersusTime(*spectrum2d);
		var.str("");
		sname.str("");
		delete spectrum2d;
		
		if(usingRealSimData) // only if this is a sim dataset
		{
		  long long int nPoints;
		  spectrum2d = new TH2F("spectrum2d","spectrum2d",100,0,1,100,0,15);
		  var << "-(RealZ-" << CurrentCrystal->GetDimensionZ()/2.0 << "):FloodZ >> spectrum2d"; 
		  nPoints = tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CutCrystal+PhotopeakEnergyCut,"COLZ"); // a 2d plot of real vs. w, using the same cuts as before
		  sname << "Real Z vs. W - Crystal " << CurrentCrystal->GetID();
		  spectrum2d->SetName(sname.str().c_str()); 
		  spectrum2d->SetTitle(sname.str().c_str());
		  spectrum2d->GetXaxis()->SetTitle("W");
		  spectrum2d->GetYaxis()->SetTitle("Z");
		  CurrentCrystal->SetSimDOIplot(*spectrum2d);
		  sname.str("");
		  sname << "Graph Z vs. W - Crystal " << CurrentCrystal->GetID();
		  simGraph = new TGraph(nPoints,tree->GetV2(),tree->GetV1()); // same but TGraph (so it can be fitted in 1D)
		  simGraph->SetName(sname.str().c_str()); 
		  simGraph->SetTitle(sname.str().c_str());
		  simGraph->GetXaxis()->SetTitle("W");
		  simGraph->GetYaxis()->SetTitle("Z");
		  simGraph->Draw("ap");
		  TF1 *linear = new TF1("linear",  "[0]*x + [1]",0,1);
		  TF1 *expfit = new TF1("expfit",  "[0]*exp(-x/[1])",0,1);
		  linear->SetParameter(0,-100);
		  linear->SetParameter(1,50);
		  expfit->SetParameter(0,50);
		  expfit->SetParameter(1,0.1);
// 		  simGraph->SetStats(1);
		  simGraph->Fit("expfit","Q","",0.1,0.7);
		  
		  CurrentCrystal->SetSimFit(*expfit);
		  CurrentCrystal->SetSimGraph(*simGraph);
		  sname.str("");
		  var.str("");
		  delete linear;
		  delete expfit;
		  delete spectrum2d;
		  delete simGraph;
		}
		
	      }
	      
	      // perform it always if this is a sim dataset
// 	      if(usingRealSimData)
// 	      {
// // 		std::cout << "Generating Z vs. W plots for crystal " << CurrentCrystal->GetID() << " ..." << std::endl;
// 		std::stringstream simCutCrystal;
// 		simCutCrystal << "RealX > " << CurrentCrystal->GetX() - CurrentCrystal->GetDimensionX()/2.0 << " && RealX < " << CurrentCrystal->GetX() + CurrentCrystal->GetDimensionX()/2.0 << "  && RealY > " << CurrentCrystal->GetY() - CurrentCrystal->GetDimensionY()/2.0 << " && RealY < " << CurrentCrystal->GetY() + CurrentCrystal->GetDimensionY()/2.0;
// 		
// 		TCut CutCrystal = simCutCrystal.str().c_str();
// // 		std::cout << CutCrystal << std::endl;
// 		long long int nPoints;
// 		spectrum2d = new TH2F("spectrum2d","spectrum2d",100,0,1,100,0,15);
// 		var << "-(RealZ-" << CurrentCrystal->GetDimensionZ()/2.0 << "):FloodZ >> spectrum2d"; 
// 		nPoints = tree->Draw(var.str().c_str(),CutCrystal,"COLZ");
// 		sname << "Real Z vs. W - Crystal " << CurrentCrystal->GetID();
// 		spectrum2d->SetName(sname.str().c_str()); 
// 		spectrum2d->SetTitle(sname.str().c_str());
// 		spectrum2d->GetXaxis()->SetTitle("W");
// 		spectrum2d->GetYaxis()->SetTitle("Z");
// 		CurrentCrystal->SetSimDOIplot(*spectrum2d);
// 		sname.str("");
// 		sname << "Graph Z vs. W - Crystal " << CurrentCrystal->GetID();
// 		simGraph = new TGraph(nPoints,tree->GetV2(),tree->GetV1());
// 		simGraph->SetName(sname.str().c_str()); 
// 		simGraph->SetTitle(sname.str().c_str());
// 		simGraph->GetXaxis()->SetTitle("W");
// 		simGraph->GetYaxis()->SetTitle("Z");
// 		simGraph->Draw("ap");
// 		TF1 *linear = new TF1("linear",  "[0]*x + [1]",0,1);
// 		TF1 *expfit = new TF1("expfit",  "[0]*exp(x/[1])",0,1);
// 		linear->SetParameter(0,-100);
// 		linear->SetParameter(1,50);
// 		simGraph->Fit("linear","Q","",0,1);
// 		
// 		CurrentCrystal->SetSimFit(*linear);
// 		CurrentCrystal->SetSimGraph(*simGraph);
// 		sname.str("");
// 		var.str("");
// // 		delete linear;
// 		delete spectrum2d;
// 		delete simGraph;
// 	      }
	      
	    }
	  }
	  
	}
      }
    }
  }
  //----------------------------------------------------------//
  
  
  
  //create the crystal 
  //   Crystal* 
  //   
  //   if(usingTaggingBench) 
  //   {
  //     
  //     // raw spectrum
  //     spectrum = new TH1F("spectrum","spectrum",histo1Dbins,1,histo1Dmax);
  //     channel = mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel();
  //     var << "ch" << channel << " >> spectrum";
  //     tree->Draw(var.str().c_str(),"");
  //     name = "Raw Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
  //     spectrum->SetName(name);
  //     spectrum->SetTitle(name);
  //     spectrum->GetXaxis()->SetTitle("ADC Channels");
  //     spectrum->GetYaxis()->SetTitle("N");
  //     mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetRawSpectrum(*spectrum);
  //     var.str("");
  //     delete spectrum;
  //     //plot the tagging bench spectrum
  //     TaggingCanvas = new TCanvas("RawSpectra","Rawspectra",1200,800);
  //     
  //   }
  
  
  //----------------------------------------------------------//
  // Produce some Canvases                                    //
  //----------------------------------------------------------//
  //multicanvases per module 0,0
  //FIXME extend to multiple modules please...
  TCanvas* RawCanvas = new TCanvas("RawSpectra","Rawspectra",1200,800);
  TCanvas* TriggerCanvas = new TCanvas("TriggerSpectra","TriggerSpectra",1200,800);
  TCanvas* FloodHistoCanvas = new TCanvas("FloodHisto","FloodHisto",800,800);
  TCanvas* FloodHisto3DCanvas = new TCanvas("FloodHisto3D","FloodHisto3D",800,800);
  TCanvas* SphericalCanvas = new TCanvas("Spherical","Spherical",1200,800);
  TCanvas* CylindricalXCanvas = new TCanvas("CylindricalX","CylindricalX",1200,800);
  TCanvas* CylindricalYCanvas = new TCanvas("CylindricalY","CylindricalY",1200,800);
  RawCanvas->Divide(4,4);
  TriggerCanvas->Divide(4,4);
  FloodHistoCanvas->Divide(4,4);
  FloodHisto3DCanvas->Divide(4,4);
  SphericalCanvas->Divide(4,4);
  CylindricalXCanvas->Divide(4,4);
  CylindricalYCanvas->Divide(4,4);
  //canvases for the global plots
  TCanvas* GlobalFlood2D = new TCanvas("Flood Histogram 2D","Flood Histogram 2D",800,800);
  TCanvas* GlobalFlood3D = new TCanvas("Flood Histogram 3D","Flood Histogram 3D",800,800);
  TCanvas* GlobalSpherical = new TCanvas("Spherical Plot","Spherical Plot",1200,800);
  TCanvas* GlobalCylindricalX = new TCanvas("Cylindrical Plot Theta:X","Cylindrical Plot Theta:X",1200,800);
  TCanvas* GlobalCylindricalY = new TCanvas("Cylindrical Plot Theta:Y","Cylindrical Plot Theta:Y",1200,800);
  //canvas for the tagging crystal
  TCanvas* TaggingCanvas = new TCanvas("Tagging Crystal","Tagging Crystal",1200,800);    
  //draw canvases
  std::cout << "Saving data to " << outputFileName << " ..." << std::endl;
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      GlobalFlood2D->cd();
      module[iModule][jModule]->GetFloodMap2D()->Draw("COLZ");
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)   
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  for(int iCry = 0; iCry < ncrystalsx ; iCry++)
	  {
	    for(int jCry = 0; jCry < ncrystalsy ; jCry++)
	    {
	      Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
	      if(CurrentCrystal->CrystalIsOn())
	      {
		CurrentCrystal->GetGraphicalCut()->SetFillStyle(4001);
		CurrentCrystal->GetGraphicalCut()->SetLineColor(kRed);
		CurrentCrystal->GetGraphicalCut()->SetLineWidth(2);
		CurrentCrystal->GetGraphicalCut()->Draw("same");
	      }
	    }
	  }
	}
      } 
      GlobalFlood3D->cd();
      module[iModule][jModule]->GetFloodMap3D()->Draw();
      GlobalSpherical->cd();
      module[iModule][jModule]->GetSphericalMap()->Draw("COLZ");
      GlobalCylindricalX->cd();
      module[iModule][jModule]->GetCylindricalXMap()->Draw("COLZ");
      GlobalCylindricalY->cd();
      module[iModule][jModule]->GetCylindricalYMap()->Draw("COLZ");
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)   
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  RawCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()); 
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetRawSpectrum()->Draw(); 
	  TriggerCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()); 
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetTriggerSpectrum()->Draw();
	  FloodHistoCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()); 
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap2D()->Draw("COLZ");
	  FloodHisto3DCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()); 
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap3D()->Draw("COLZ");	  
	  SphericalCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()); 
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetSphericalMap()->Draw("COLZ");
	  CylindricalXCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()); 
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCylindricalXMap()->Draw("COLZ");	  
	  CylindricalYCanvas->cd(mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCanvasPosition()); 
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetCylindricalYMap()->Draw("COLZ");
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
  TH1F *PeakPositionDistro = new TH1F("Distribution of photopeak positions","Distribution photopeak positions",100,0,12000);
  PeakPositionDistro->GetXaxis()->SetTitle("ADC Channels");
  PeakPositionDistro->GetYaxis()->SetTitle("N");
  PeakPositionDistro->SetStats(1);
  //2d histogram
  TH2F *PeakPositionVsIJ = new TH2F("Distribution of photopeak positions VS. crystal position i,j","Distribution photopeak positions VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  PeakPositionVsIJ->GetXaxis()->SetTitle("i");
  PeakPositionVsIJ->GetYaxis()->SetTitle("j");
  PeakPositionVsIJ->GetZaxis()->SetTitle("ADC Channels");
//   PeakPositionVsIJ->SetStats(1);
  //--Distribution of energy resolutions FHWM
  //histogram
  TH1F *PeakEnergyResolutionDistro = new TH1F("Distribution of photopeak energy resolutions FWHM","Distribution photopeak energy resolutions FWHM",100,0,1);
  PeakEnergyResolutionDistro->GetXaxis()->SetTitle("Energy Resolution FWHM");
  PeakEnergyResolutionDistro->GetYaxis()->SetTitle("N");
  PeakEnergyResolutionDistro->SetStats(1);
  //2d histogram
  TH2F *EnergyResolutionVsIJ = new TH2F("Distribution of photopeak energy resolutions FWHM VS. crystal position i,j","Distribution photopeak energy resolutions FWHM VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  EnergyResolutionVsIJ->GetXaxis()->SetTitle("i");
  EnergyResolutionVsIJ->GetYaxis()->SetTitle("j");
  EnergyResolutionVsIJ->GetZaxis()->SetTitle("En. Res.");
//   EnergyResolutionVsIJ->SetStats(1);
  //Distribution of FWHM of W plots
  //histogram of fwhm
  TH1F *WfwhmDistro = new TH1F("Distribution of FWHM in W plots","Distribution of FWHM in W plots",100,0,0.5);
  WfwhmDistro->GetXaxis()->SetTitle("W");
  WfwhmDistro->GetYaxis()->SetTitle("N");
  WfwhmDistro->SetStats(1);
  //2d histogram of fwhm
  TH2F *WfwhmVsIJ = new TH2F("Distribution of FWHM in W plots VS. crystal position i,j","Distribution of FWHM in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  WfwhmVsIJ->GetXaxis()->SetTitle("i");
  WfwhmVsIJ->GetYaxis()->SetTitle("j");
  WfwhmVsIJ->GetZaxis()->SetTitle("w FHWM");
//   WfwhmVsIJ->SetStats(1);
  //histogram of rms
  TH1F *WrmsDistro = new TH1F("Distribution of RMS in W plots","Distribution of RMS in W plots",100,0,0.5);
  WrmsDistro->GetXaxis()->SetTitle("W");
  WrmsDistro->GetYaxis()->SetTitle("N");
  WrmsDistro->SetStats(1);
  //2d histogram of rms
  TH2F *WrmsVsIJ = new TH2F("Distribution of RMS in W plots VS. crystal position i,j","Distribution of RMS in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  WrmsVsIJ->GetXaxis()->SetTitle("i");
  WrmsVsIJ->GetYaxis()->SetTitle("j");
  WrmsVsIJ->GetZaxis()->SetTitle("w RMS");
  //Distribution of FWHM of W plots
  TH1F *Wwidth20perc = new TH1F("Distribution of width at 20% in W plots","Distribution of width at 20% in W plots",100,0,0.5);
  Wwidth20perc->GetXaxis()->SetTitle("W");
  Wwidth20perc->GetYaxis()->SetTitle("N");
  Wwidth20perc->SetStats(1);
  
  
  TH2F *Wwidht20percVsIJ = new TH2F("Distribution of width at 20% in W plots VS. crystal position i,j","Distribution of width at 20% in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  Wwidht20percVsIJ->GetXaxis()->SetTitle("i");
  Wwidht20percVsIJ->GetYaxis()->SetTitle("i");
  Wwidht20percVsIJ->GetZaxis()->SetTitle("w width at 20%");
  
  
  //Distribution of DOI resolutions - not very nice since one parameter in the calculation is assumed (from the DOI bench results)
  TH1F *WDoiDistro = new TH1F("Distribution of doi res","Distribution of doi res",20,0,6);
  WDoiDistro->GetXaxis()->SetTitle("doi");
  WDoiDistro->GetYaxis()->SetTitle("N");
  WDoiDistro->SetStats(1);
  
  //Distribution of fit exp tau for W plots
  TH1F *WtauFit = new TH1F("Distribution of exp slopes in W plots","Distribution of exp slopes in W plots",1000,0,1);
  WtauFit->GetXaxis()->SetTitle("Tau");
  WtauFit->GetYaxis()->SetTitle("N");
  WtauFit->SetStats(1);
  
  
  TH2F *WtauFitVsIJ = new TH2F("Distribution of exp slopes in W plots VS. crystal position i,j","Distribution of exp slopes in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
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
      directory[iModule+jModule][0][0] = fPlots->mkdir(ModuleDirStream.str().c_str());
      directory[iModule+jModule][0][0]->cd();      
      GlobalFlood2D->Write();
      GlobalFlood3D->Write();
      GlobalSpherical->Write();
      GlobalCylindricalX->Write();
      GlobalCylindricalY->Write();
      RawCanvas->Write();
      TriggerCanvas->Write();
      FloodHistoCanvas->Write();
      FloodHisto3DCanvas->Write();
      SphericalCanvas->Write();
      CylindricalXCanvas->Write();
      CylindricalYCanvas->Write();
      //save the 3d flood maps separately for each channel
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  std::stringstream MppcDirStream;
	  MppcDirStream << "MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " - " <<  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetExtendedID();
	  directory[iModule+jModule][(iMppc+jMppc)+1][0] = directory[iModule+jModule][0][0]->mkdir(MppcDirStream.str().c_str());
	  directory[iModule+jModule][(iMppc+jMppc)+1][0]->cd();
	  mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetFloodMap3D()->Write();
	  for(int iCry = 0; iCry < ncrystalsx ; iCry++)
	  {
	    for(int jCry = 0; jCry < ncrystalsy ; jCry++)
	    {
	      std::stringstream CrystalDirStream;
	      Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
	      CrystalDirStream << "Crystal " <<  CurrentCrystal->GetID();
	      directory[iModule+jModule][(iMppc+jMppc)+1][(iCry+jCry)+1] = directory[iModule+jModule][(iMppc+jMppc)+1][0]->mkdir(CrystalDirStream.str().c_str());
	      directory[iModule+jModule][(iMppc+jMppc)+1][(iCry+jCry)+1]->cd(); 
	      if(CurrentCrystal->CrystalIsOn() /*| usingRealSimData*/) // save data only if the crystal was specified in the config file
	      {
		//create a pointer for the current crystal (mainly to make the code more readable)
		Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
		//fill the global distributions histograms
		PeakPositionDistro->Fill(CurrentCrystal->GetPhotopeakPosition());
		PeakEnergyResolutionDistro->Fill(CurrentCrystal->GetPhotopeakEnergyResolution());
		WfwhmDistro->Fill(CurrentCrystal->GetWfwhm());
		WDoiDistro->Fill( (15.0/CurrentCrystal->GetWfwhm())*0.0158); // CAREFUL: here the 0.0158 value is hardcoded and taken from the sigma of W distros in DOI bench setup. 15.0 is the length of the crystals in mm.
		PeakPositionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakPosition());
		EnergyResolutionVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetPhotopeakEnergyResolution());
		WfwhmVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetWfwhm());
		WrmsDistro->Fill(CurrentCrystal->GetWrms());
		WrmsVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetWrms());
		Wwidht20percVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetWwidth20perc());
		Wwidth20perc->Fill(CurrentCrystal->GetWwidth20perc());
		
		C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
		C_spectrum->SetName(CurrentCrystal->GetSpectrum()->GetName());
		C_spectrum->cd();
		CurrentCrystal->GetSpectrum()->Draw();
		CurrentCrystal->GetHighlightedSpectrum()->SetFillColor(3);
		CurrentCrystal->GetHighlightedSpectrum()->Draw("same");
		CurrentCrystal->GetFit()->Draw("same");
		C_spectrum->Write();
		delete C_spectrum;
		
		C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
		C_spectrum->SetName(CurrentCrystal->GetHistoW()->GetName());
		C_spectrum->cd();
		CurrentCrystal->GetHistoW()->Draw();
		C_spectrum->Write();
		delete C_spectrum;
		
		C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
		C_spectrum->SetName(CurrentCrystal->GetFloodMap2D()->GetName());
		C_spectrum->cd();
		CurrentCrystal->GetFloodMap2D()->Draw("COLZ");
		C_spectrum->Write();
		delete C_spectrum;
		
		C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
		C_spectrum->SetName(CurrentCrystal->GetVersusTime()->GetName());
		C_spectrum->cd();
		CurrentCrystal->GetVersusTime()->Draw("COLZ");
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
	}
      }
      directory[iModule+jModule][0][0]->cd(); // go back to main directory
      //write the summary histos and canvases, that were filled during file saving
      PeakPositionDistro->Write();
      PeakEnergyResolutionDistro->Write();
      WfwhmDistro->Write();
      WDoiDistro->Write();
      WrmsDistro->Write();
      Wwidth20perc->Write();
      
      TCanvas *C_Wwidht20percVsIJ = new TCanvas("C_Wwidht20percVsIJ","C_Wwidht20percVsIJ",800,800);
      C_Wwidht20percVsIJ->SetName(Wwidht20percVsIJ->GetName());
      C_Wwidht20percVsIJ->cd();
      Wwidht20percVsIJ->Draw("LEGO2");
      C_Wwidht20percVsIJ->Write();
      
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
      C_PeakPositionVsIJ->Write();
      
      TCanvas *C_EnergyResolutionVsIJ = new TCanvas("C_EnergyResolutionVsIJ","C_EnergyResolutionVsIJ",800,800);
      C_EnergyResolutionVsIJ->SetName(EnergyResolutionVsIJ->GetName());
      C_EnergyResolutionVsIJ->cd();
      EnergyResolutionVsIJ->Draw("LEGO2");
      C_EnergyResolutionVsIJ->Write();
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
  
  
  delete crystal;
  delete mppc;
  delete module;
  
  return 0;
}