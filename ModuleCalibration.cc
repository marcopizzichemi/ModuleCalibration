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




int main (int argc, char** argv)
{
  //----------------------------------------------------------//
  //  Check input                                             //
  //----------------------------------------------------------//
  if(argc<2) 
  {
    std::cout << " Usage: " 							<< std::endl;
    std::cout << " ModuleCalibration [-c config file ] <input-files> " 		<< std::endl;
    std::cout << "   note: -c option is optional, but if you use it, it has to be the first argument"	<< std::endl;
    std::cout << "         without it, a configFile.cfg will be assumed" 	<< std::endl;
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
  
  int ncrystalsx = config.read<int>("ncrystalsx");
  int ncrystalsy = config.read<int>("ncrystalsy");
  int nmppcx     = config.read<int>("nmppcx");
  int nmppcy     = config.read<int>("nmppcy");
  int nmodulex   = config.read<int>("nmodulex");
  int nmoduley   = config.read<int>("nmoduley");
  
  int histo1Dmax        = config.read<int>("histo1Dmax");       
  int histo1Dbins       = config.read<int>("histo1Dbins");      
  int histo2DchannelBin = config.read<int>("histo2DchannelBin");
  int histo2DglobalBins = config.read<int>("histo2DglobalBins");
  int histo3DchannelBin = config.read<int>("histo3DchannelBin");
  int histo3DglobalBins = config.read<int>("histo3DglobalBins");
  
  std::string outputFileName = config.read<std::string>("output");
  outputFileName += ".root";
  
  Module*** module;                           // create the elements
  Mppc*** mppc;
  Crystal*** crystal;
  module = new Module**[nmodulex]; // make an array of module pointers
  for(int j = 0; j < nmodulex ; j++)
  {
    module[j] = new Module* [nmoduley];
  }
  mppc = new Mppc**[nmodulex*nmppcx]; // make an array of mppc pointers
  for(int j = 0; j < nmodulex*nmppcx ; j++)
  {
    mppc[j] = new Mppc* [nmoduley*nmppcy];
  }
  crystal = new Crystal**[nmodulex*nmppcx*ncrystalsx]; // make an array of crystal pointers
  for(int j = 0; j < nmodulex*nmppcx*ncrystalsx ; j++)
  {
    crystal[j] = new Crystal* [nmoduley*nmppcy*ncrystalsy];
  }
  
  
  input.FillElements(module,mppc,crystal);    // fill the elements
  
  
  TTree* tree = input.GetTree();     // get the TTree, to plot the spectra
 
  TH1F* spectrum;
  TH2F* spectrum2d;
  TH3F* spectrum3d;
  std::stringstream var,cut;
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
      spectrum3d->GetYaxis()->SetTitle("W");
      module[iModule][jModule]->SetFloodMap3D(*spectrum3d);
      delete spectrum3d;
      
      //spectra for each mppc
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  
	  int channel = mppc[iMppc][jMppc]->GetDigitizerChannel();
	  std::cout << "Generating spectra for MPPC " << mppc[iMppc][jMppc]->GetLabel() << " ..." << std::endl;
	  
	  cut << "TriggerChannel == " << channel  ;
	  TCut CutTrigger = cut.str().c_str();
	  cut.str("");
	  //same as the global ones, but selecting on TriggerChannel
	  
	  
	  // raw spectrum
	  spectrum = new TH1F("spectrum","spectrum",histo1Dbins,1,histo1Dmax);
	  channel = mppc[iMppc][jMppc]->GetDigitizerChannel();
	  var << "ch" << channel << " >> spectrum";
	  tree->Draw(var.str().c_str(),"");
	  name = "Raw Spectrum - MPPC " + mppc[iMppc][jMppc]->GetLabel();
	  spectrum->SetName(name);
	  spectrum->SetTitle(name);
	  spectrum->GetXaxis()->SetTitle("ADC Channels");
          spectrum->GetYaxis()->SetTitle("N");
	  mppc[iMppc][jMppc]->SetRawSpectrum(*spectrum);
	  var.str("");
	  delete spectrum;
	  //trigger selected spectrum
	  spectrum = new TH1F("spectrum","spectrum",histo1Dbins,1,histo1Dmax);	  
	  var << "ch" << channel << " >> spectrum";
	  tree->Draw(var.str().c_str(),CutTrigger);
	  name = "Trigger Spectrum - MPPC " + mppc[iMppc][jMppc]->GetLabel();
	  spectrum->SetName(name);
	  spectrum->SetTitle(name);
	  spectrum->GetXaxis()->SetTitle("ADC Channels");
          spectrum->GetYaxis()->SetTitle("N");
	  mppc[iMppc][jMppc]->SetTriggerSpectrum(*spectrum);
	  var.str("");
	  delete spectrum;
	  
	  // Flood histogram
	  spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DchannelBin,-7,7,histo2DchannelBin,-7,7);
	  tree->Draw("FloodY:FloodX >> spectrum2d",CutTrigger,"COLZ");
	  name = "Flood Histogram 2D - MPPC " + mppc[iMppc][jMppc]->GetLabel();
	  spectrum2d->SetName(name); 
	  spectrum2d->SetTitle(name);
	  spectrum2d->GetXaxis()->SetTitle("U");
	  spectrum2d->GetYaxis()->SetTitle("V");
	  mppc[iMppc][jMppc]->SetFloodMap2D(*spectrum2d);
	  delete spectrum2d;
	  //sherical coordinates plot
	  spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DchannelBin,1.1,3.1415/2.0,histo2DchannelBin,-3.14/2.0,3.14/2.0); 
	  tree->Draw("Phi:Theta >> spectrum2d",CutXYZ+CutTrigger,"COLZ");
	  name = "Spherical Plot - MPPC " + mppc[iMppc][jMppc]->GetLabel();
	  spectrum2d->SetName(name);
	  spectrum2d->SetTitle(name);
	  spectrum2d->GetXaxis()->SetTitle("Theta");
	  spectrum2d->GetYaxis()->SetTitle("Phi");
	  mppc[iMppc][jMppc]->SetSphericalMap(*spectrum2d);
	  delete spectrum2d;
	  //cylindrical coordinates plot, x and theta
	  spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DchannelBin,-7,7,histo2DchannelBin,1.1,3.1415/2.0); 
	  tree->Draw("Theta:FloodX >> spectrum2d",CutXYZ+CutTrigger,"COLZ");
	  name = "Cylindrical Plot Theta:X - MPPC " + mppc[iMppc][jMppc]->GetLabel();
	  spectrum2d->SetName(name);
	  spectrum2d->SetTitle(name);
	  spectrum2d->GetXaxis()->SetTitle("U");
	  spectrum2d->GetYaxis()->SetTitle("Theta");
	  mppc[iMppc][jMppc]->SetCylindricalXMap(*spectrum2d);
	  delete spectrum2d;
	  //cylindrical coordinates plot, y and theta
	  spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DchannelBin,-7,7,histo2DchannelBin,1.1,3.1415/2.0); 
	  tree->Draw("Theta:FloodY >> spectrum2d",CutXYZ+CutTrigger,"COLZ");
	  name = "Cylindrical Plot Theta:Y - MPPC " + mppc[iMppc][jMppc]->GetLabel();
	  spectrum2d->SetName(name);
	  spectrum2d->SetTitle(name);
	  spectrum2d->GetXaxis()->SetTitle("V");
	  spectrum2d->GetYaxis()->SetTitle("Theta");
	  mppc[iMppc][jMppc]->SetCylindricalYMap(*spectrum2d);
	  delete spectrum2d;
	  
	  
	  
	  spectrum3d = new TH3F("spectrum3d","spectrum3d",histo3DchannelBin,-7,7,histo3DchannelBin,-7,7,histo3DchannelBin,0,1);
	  tree->Draw("FloodZ:FloodY:FloodX >> spectrum3d",CutXYZ+CutTrigger);
	  name = "Flood Histogram 3D - MPPC " + mppc[iMppc][jMppc]->GetLabel();
	  spectrum3d->SetName(name);
	  spectrum3d->SetTitle(name);
	  spectrum3d->GetXaxis()->SetTitle("U");
	  spectrum3d->GetYaxis()->SetTitle("V");
	  spectrum3d->GetYaxis()->SetTitle("W");
	  mppc[iMppc][jMppc]->SetFloodMap3D(*spectrum3d);
	  delete spectrum3d;
	  
	  
	}
      }
    }
  }
  
  //multicanvases
  TCanvas* RawCanvas = new TCanvas("RawSpectra","Rawspectra",1200,800);
  TCanvas* TriggerCanvas = new TCanvas("TriggerSpectra","TriggerSpectra",1200,800);
  TCanvas* FloodHistoCanvas = new TCanvas("FloodHisto","FloodHisto",1200,800);
  TCanvas* FloodHisto3DCanvas = new TCanvas("FloodHisto3D","FloodHisto3D",1200,800);
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
  TCanvas* GlobalFlood2D = new TCanvas("Flood Histogram 2D","Flood Histogram 2D",1200,800);
  GlobalFlood2D->cd();
  module[0][0]->GetFloodMap2D()->Draw("COLZ");
  TCanvas* GlobalFlood3D = new TCanvas("Flood Histogram 3D","Flood Histogram 3D",1200,800);
  GlobalFlood3D->cd();
  module[0][0]->GetFloodMap3D()->Draw();
  TCanvas* GlobalSpherical = new TCanvas("Spherical Plot","Spherical Plot",1200,800);
  GlobalSpherical->cd();
  module[0][0]->GetSphericalMap()->Draw("COLZ");
  TCanvas* GlobalCylindricalX = new TCanvas("Cylindrical Plot Theta:X","Cylindrical Plot Theta:X",1200,800);
  GlobalCylindricalX->cd();
  module[0][0]->GetCylindricalXMap()->Draw("COLZ");
  TCanvas* GlobalCylindricalY = new TCanvas("Cylindrical Plot Theta:Y","Cylindrical Plot Theta:Y",1200,800);
  GlobalCylindricalY->cd();
  module[0][0]->GetCylindricalYMap()->Draw("COLZ");
  
  
  
  
  std::cout << "Saving data to " << outputFileName << " ..." << std::endl;
  for(int iMppc = 0; iMppc < nmppcx ; iMppc++)   
  {
    for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
    {
      RawCanvas->cd(mppc[iMppc][jMppc]->GetCanvasPosition()); 
      mppc[iMppc][jMppc]->GetRawSpectrum()->Draw();
      
      TriggerCanvas->cd(mppc[iMppc][jMppc]->GetCanvasPosition()); 
      mppc[iMppc][jMppc]->GetTriggerSpectrum()->Draw();
      
      FloodHistoCanvas->cd(mppc[iMppc][jMppc]->GetCanvasPosition()); 
      mppc[iMppc][jMppc]->GetFloodMap2D()->Draw("COLZ");
      
      FloodHisto3DCanvas->cd(mppc[iMppc][jMppc]->GetCanvasPosition()); 
      mppc[iMppc][jMppc]->GetFloodMap3D()->Draw("COLZ");
      
      SphericalCanvas->cd(mppc[iMppc][jMppc]->GetCanvasPosition()); 
      mppc[iMppc][jMppc]->GetSphericalMap()->Draw("COLZ");
      
      CylindricalXCanvas->cd(mppc[iMppc][jMppc]->GetCanvasPosition()); 
      mppc[iMppc][jMppc]->GetCylindricalXMap()->Draw("COLZ");
      
      CylindricalYCanvas->cd(mppc[iMppc][jMppc]->GetCanvasPosition()); 
      mppc[iMppc][jMppc]->GetCylindricalYMap()->Draw("COLZ");
      
    }
  }
  
  

  TFile* fPlots = new TFile(outputFileName.c_str(),"recreate");
  fPlots->cd();
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      
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
	  mppc[iMppc][jMppc]->GetFloodMap3D()->Write();
	}
      }
    }
  }
  
  
  
  TFile* fFile = new TFile("temp.root","recreate");
  fFile->cd();
  tree->Write();
  fFile->Close();
  
  fPlots->Close();
  //----------------------------
  
  delete crystal;
  delete mppc;
  delete module;
  
  return 0;
}