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
#define ENERGY_RESOLUTION_SATURATION_CORRECTION 0.35


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
  
  ConfigFile config(ConfigFileName);
  InputFile input(argc,argv,config); // read the input chain of root files,
  input.CreateTree();                // create the TTree that will be used in analysis
  
  int ncrystalsx            = config.read<int>("ncrystalsx",2);
  int ncrystalsy            = config.read<int>("ncrystalsy");
  int nmppcx                = config.read<int>("nmppcx");
  int nmppcy                = config.read<int>("nmppcy");
  int nmodulex              = config.read<int>("nmodulex");
  int nmoduley              = config.read<int>("nmoduley");
  int histo1Dmax            = config.read<int>("histo1Dmax");       
  int histo1Dbins           = config.read<int>("histo1Dbins");      
  int histo2DchannelBin     = config.read<int>("histo2DchannelBin");
  int histo2DglobalBins     = config.read<int>("histo2DglobalBins");
  int histo3DchannelBin     = config.read<int>("histo3DchannelBin");
  int histo3DglobalBins     = config.read<int>("histo3DglobalBins");
  bool saveAnalysisTree     = config.read<bool>("saveAnalysisTree");
  float taggingPosition     = config.read<float>("taggingPosition");
  bool usingTaggingBench    = config.read<bool>("usingTaggingBench");
  int taggingCrystalChannel = config.read<int>("taggingCrystalChannel");
  bool correctingSaturation = config.read<bool>("correctingSaturation");;
  // set output file name
  std::string outputFileName = config.read<std::string>("output");
  outputFileName += ".root";
  // digitizer channels
  std::string digitizer_s   = config.read<std::string>("digitizer");
  std::vector<std::string> digitizer_f;
  config.split( digitizer_f, digitizer_s, "," );
  std::vector<int> digitizer;
  for(int i = 0 ; i < digitizer_f.size() ; i++)
  {
    config.trim(digitizer_f[i]);
    digitizer.push_back(atoi(digitizer_f[i].c_str()));
  }
  //create "sum channels" string
  std::stringstream sSumChannels;
  std::string SumChannels;
  sSumChannels << "ch" <<  digitizer[0];
  for(int i = 1 ; i < digitizer.size() ; i++)
    sSumChannels << "+ch" <<digitizer[i];
  SumChannels = sSumChannels.str();
  
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
	  name = "Raw Spectrum - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
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
	      if(CurrentCrystal->CrystalIsOn())
	      {
		std::cout << "Generating spectra for crystal " << CurrentCrystal->GetID() << " ..." << std::endl;
		TCut CutCrystal = CurrentCrystal->GetCrystalCut();
		// 	      std::cout << CutCrystal << std::endl;
		
		//-------------------------------------------------------------------------
		//standard sum spectrum with cut on crystal events, xyz and trigger channel
		//------------------------------------------------------------------------
		//draw spectrum
		spectrum = new TH1F("spectrum","spectrum",histo1Dbins,1,histo1Dmax);	  
		var << SumChannels << " >> spectrum";
		tree->Draw(var.str().c_str(),CutXYZ+CutTrigger+CutCrystal);
		sname << "Charge Spectrum - Crystal " << CurrentCrystal->GetID();
		spectrum->SetName(sname.str().c_str());
		spectrum->SetTitle(sname.str().c_str());
		spectrum->GetXaxis()->SetTitle("ADC Channels");
		spectrum->GetYaxis()->SetTitle("N");
		//automatically look for the 511Kev peak to find the photopeak energy cut
		//find peaks in each crystal spectrum, with TSpectrum
		TSpectrum *s;
		s = new TSpectrum(5);
// 		Input[i].SumSpectraCanvas->cd(j+1);
		Int_t CrystalPeaksN = s->Search(spectrum,1,"goff",0.5); 
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
		float energyResolution;
		if (correctingSaturation)
		  energyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION; 
		else
		  energyResolution = ENERGY_RESOLUTION;
		float par0 = CrystalPeaksY[peakID];
		float par1 = CrystalPeaks[peakID];
		float par2 = (CrystalPeaks[peakID]*energyResolution)/2.35;
		float fitmin = par1-1*par2;
		float fitmax = par1+1.8*par2;
		TF1 *gauss = new TF1("gauss",  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitmin,fitmax);
		gauss->SetParameter(0,par0);
		gauss->SetParameter(1,par1);
		gauss->SetParameter(2,par2); //expected FWHM en res = 12%
		spectrum->Fit("gauss","Q","",fitmin,fitmax);
		//store the mean and sigma in the crystal
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
		int bin1 = spectrum->FindFirstBinAbove(spectrum->GetMaximum()/2);
		int bin2 = spectrum->FindLastBinAbove(spectrum->GetMaximum()/2);
		double fwhm = spectrum->GetBinCenter(bin2) - spectrum->GetBinCenter(bin1);
		CurrentCrystal->SetHistoW(*spectrum);
		CurrentCrystal->SetHistoWfwhm(fwhm);
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
		
	      }
	      
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
  // Canvases                                                 //
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
  // Global distributions                                     //
  //----------------------------------------------------------//
  TH1F *PeakPositionDistro = new TH1F("Distribution photopeak positions","Distribution photopeak positions",50,0,12000);
  PeakPositionDistro->GetXaxis()->SetTitle("ADC Channels");
  PeakPositionDistro->GetYaxis()->SetTitle("N");
  TH1F *PeakEnergyResolutionDistro = new TH1F("Distribution photopeak energy resolutions FWHM","Distribution photopeak energy resolutions FWHM",50,0,1);
  PeakEnergyResolutionDistro->GetXaxis()->SetTitle("Energy Resolution FWHM");
  PeakEnergyResolutionDistro->GetYaxis()->SetTitle("N");
  TH1F *WfwhmDistro = new TH1F("Distribution of FWHM in W plots","Distribution of FWHM in W plots",50,0,0.5);
  WfwhmDistro->GetXaxis()->SetTitle("W");
  WfwhmDistro->GetYaxis()->SetTitle("N");
  TH1F *WDoiDistro = new TH1F("Distribution of doi res","Distribution of doi res",20,0,6);
  WDoiDistro->GetXaxis()->SetTitle("doi");
  WDoiDistro->GetYaxis()->SetTitle("N");
  TH2F *WfwhmVsIJ = new TH2F("Distribution of FWHM in W plots VS. crystal position i,j","Distribution of FWHM in W plots VS. crystal position i,j",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  WfwhmVsIJ->GetXaxis()->SetTitle("i");
  WfwhmVsIJ->GetYaxis()->SetTitle("j");
  WfwhmVsIJ->GetZaxis()->SetTitle("w FHWM");
  
  
  //----------------------------------------------------------//
  // Write output to root file(s)                             //
  //----------------------------------------------------------//
  //write output plots
  TFile* fPlots = new TFile(outputFileName.c_str(),"recreate");
  fPlots->cd();
  TDirectory ****directory; //TDirectory 
  directory = new TDirectory***[nmodulex*nmoduley]; 
  for(int i = 0; i < nmodulex*nmoduley+1 ; i++) 
  {
    directory[i] = new TDirectory** [nmppcx*nmppcy+1];
    for(int j = 0; j < nmppcx*nmppcy+1 ; j++) 
    {
      directory[i][j] = new TDirectory* [(nmppcx*ncrystalsx)*(nmppcy*ncrystalsy)+1];
    }
  }
  
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
	  MppcDirStream << "MPPC " << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
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
	      
	      if(CurrentCrystal->CrystalIsOn())
	      {
		
		Crystal *CurrentCrystal = crystal[(iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry)][(jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry)];
		//fill the global distributions histograms
		PeakPositionDistro->Fill(CurrentCrystal->GetPhotopeakPosition());
		PeakEnergyResolutionDistro->Fill(CurrentCrystal->GetPhotopeakEnergyResolution());
		WfwhmDistro->Fill(CurrentCrystal->GetWfwhm());
		WDoiDistro->Fill( (15.0/CurrentCrystal->GetWfwhm())*0.0158);
		WfwhmVsIJ->Fill(CurrentCrystal->GetI(),CurrentCrystal->GetJ(),CurrentCrystal->GetWfwhm());
		
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
		CurrentCrystal->GetFloodMap2D()->Draw();
		C_spectrum->Write();
		delete C_spectrum;
		
		
		
	      }
	    }
	  }
	  
	  
	  
	}
      }
      directory[iModule+jModule][0][0]->cd();
      PeakPositionDistro->Write();
      PeakEnergyResolutionDistro->Write();
      WfwhmDistro->Write();
      WDoiDistro->Write();
      TCanvas *C_WfwhmVsIJ = new TCanvas("C_WfwhmVsIJ","C_WfwhmVsIJ",800,800);
      C_WfwhmVsIJ->SetName(WfwhmVsIJ->GetName());
      C_WfwhmVsIJ->cd();
      WfwhmVsIJ->Draw("LEGO2");
      C_WfwhmVsIJ->Write();
    }
  }
  if(saveAnalysisTree)
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