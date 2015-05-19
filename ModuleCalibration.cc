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
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TF2.h"
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
  
  
  //open the config file
  //TODO modify this in such a way that when calling 
  //confif.read on a std::vector<std::string>, it automatically does
  //the split and trim part
  ConfigFile config(ConfigFileName);
  //--------------------------------------------------------------------------------------------------------------
  //read the strings that describe the input channels
  std::string digitizer_s    = config.read<std::string>("digitizer");
  std::string mppc_s         = config.read<std::string>("mppc");
  std::string plotPositions_s = config.read<std::string>("plotPositions");
  std::string xPositions_s   = config.read<std::string>("xPositions");
  std::string yPositions_s   = config.read<std::string>("yPositions");
  //prepare std::vectors to split the strings
  std::vector <std::string> digitizer_f;
  std::vector <std::string> mppc_f;
  std::vector <std::string> plotPositions_f;
  std::vector <std::string> xPositions_f;
  std::vector <std::string> yPositions_f;
  std::vector <int> digitizer;
  std::vector <std::string> mppc_label;
  std::vector <int> plotPositions;
  std::vector <float> xPositions;
  std::vector <float> yPositions;
  //split them using the config file class
  config.split( digitizer_f, digitizer_s, "," );
  config.split( mppc_f, mppc_s, "," );
  config.split( plotPositions_f, plotPositions_s, "," );
  config.split( xPositions_f, xPositions_s, "," );
  config.split( yPositions_f, yPositions_s, "," );
  //trim them using the config file class (i.e. remove spaces)
  //and at the same time put in vectors with numbers for the ones that are numbers
  for(int i = 0 ; i < digitizer_f.size() ; i++)
  {
    config.trim(digitizer_f[i]);
    digitizer.push_back(atoi(digitizer_f[i].c_str()));
  }
  for(int i = 0 ; i < mppc_f.size() ; i++)
  {
    config.trim(mppc_f[i]);
    mppc_label.push_back(mppc_f[i]);
  }
  for(int i = 0 ; i < plotPositions_f.size() ; i++)
  {
    config.trim(plotPositions_f[i]);
    plotPositions.push_back(atoi(plotPositions_f[i].c_str()));
  }
  for(int i = 0 ; i < xPositions_f.size() ; i++)
  {
    config.trim(xPositions_f[i]);
    xPositions.push_back(atof(xPositions_f[i].c_str()));
  }
  for(int i = 0 ; i < yPositions_f.size() ; i++)
  {
    config.trim(yPositions_f[i]);
    yPositions.push_back(atof(yPositions_f[i].c_str()));
  }
  //check if the vectors just built have the same size
  assert( (digitizer.size() == mppc_label.size() ) && (digitizer.size() == plotPositions.size()) && (digitizer.size() == xPositions.size()) && (digitizer.size() == yPositions.size()) );
  
  if(digitizer.size() > 16) 
  {
    std::cout << "ERROR: Only one module can be analyzed at a time! Set 16 or less input channels in the config file!" << std::endl;
    return 1;
  }
  //feedback to the user
  std::cout << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << " Channels configuration " << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << "ADC input\tMPPC ch\tCanvas\tx[mm]\ty[mm]" << std::endl;
  std::cout << "------------------------" << std::endl;
  for(int i = 0 ; i < digitizer.size() ; i++)
  {
    std::cout << "Channel[" << digitizer[i] << "] = \t" <<  mppc_label[i] << "\t" << plotPositions[i] << "\t" << xPositions[i] << "\t" << yPositions[i] << std::endl;
  }
  std::cout << "------------------------" << std::endl;
  std::cout << std::endl;
  //-------------------------------------------------------------------------------------------
  
  //temp
  int translateCh[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  
  std::string chainName = config.read<std::string>("chainName");
  InputFile input(argc,argv,chainName,digitizer.size()); // read the input chain of root files, prepares the ttree that will be used in the analysis

  
//   TChain* fchain = input.GetChain();
  input.CreateTree(digitizer,xPositions,yPositions);
  TTree* tree = input.GetTree();
  
  //temporary save ttree to check sanity
  //-----------------------------
//   TFile* fTree = new TFile("temp.root","recreate");
//   fTree->cd();
//   tree->Write();
//   fTree->Close();
//   std::cout << "DONE" << std::endl;
  //----------------------------
  
  //--------------------------------------------------------------------------------------------------------------
  int ncrystalsx = config.read<int>("ncrystalsx");
  int ncrystalsy = config.read<int>("ncrystalsy");
  int nmppcx     = config.read<int>("nmppcx");
  int nmppcy     = config.read<int>("nmppcy");
  int nmodulex   = config.read<int>("nmodulex");
  int nmoduley   = config.read<int>("nmoduley");
  
  //----------------------------------------------------------
  //TEST of pointer  
//   Crystal *crystal = new Crystal();
//   crystal->SetName("The Crystal");
//   Mppc *mppc = new Mppc();
//   mppc->SetName("The Mppc");
//   Module *module = new Module();
//   module->SetName("The Module");
//   crystal->SetMppc( (Element*) mppc );
//   mppc->SetModule( (Element*) module );
//   mppc->MakeCrystalPointers(1,1);
//   mppc->SetCrystal(0,0, (Element*) crystal);
//   module->MakeMppcPointers(1,1);
//   module->SetMppc(0,0, (Element*) mppc);
//   std::cout << std::endl;
//   std::cout << crystal->GetMppc() << " " << mppc << " " << crystal->GetMppc()->GetName() << " " << mppc->GetName() << std::endl;
//   std::cout << mppc->GetModule() << " " << module << " "  << mppc->GetModule()->GetName() << " " << module->GetName() << std::endl;
// //   std::cout << test->GetModule()->GetName()  << " " << parent->GetName() << std::endl;
//   delete crystal;
//   delete module;
//   delete mppc;
  //----------------------------------------------------------


  int moduleCounter = 0;
  Module*** module = new Module**[nmodulex]; // make an array of module pointers
  for(int j = 0; j < nmodulex ; j++)
  {
    module[j] = new Module* [nmoduley];
  }
  
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      std::stringstream sname;
      sname << "Module " << iModule << "." << jModule;
      module[iModule][jModule] = new Module(); // creates a default module   
      module[iModule][jModule]->SetName(sname.str().c_str());          // assign a name
      module[iModule][jModule]->SetID(moduleCounter);                  // assign an ID number
      module[iModule][jModule]->SetI(iModule); 
      module[iModule][jModule]->SetJ(jModule);
      module[iModule][jModule]->SetChildrenI(nmppcx);
      module[iModule][jModule]->SetChildrenJ(nmppcy);
      moduleCounter++;
    }
  }
  
  int mppcCounter = 0;
  Mppc*** mppc = new Mppc**[nmodulex*nmppcx]; // make an array of module pointers
  for(int j = 0; j < nmodulex*nmppcx ; j++)
  {
    mppc[j] = new Mppc* [nmoduley*nmppcy];
  }
  for(int iMppc = 0; iMppc < nmppcx*nmodulex ; iMppc++)
  {
    for(int jMppc = 0; jMppc < nmppcy*nmoduley ; jMppc++)
    {
      std::stringstream sname;
      sname << "Mppc " << mppc_label[mppcCounter];
      mppc[iMppc][jMppc] = new Mppc();
      mppc[iMppc][jMppc]->SetName(sname.str().c_str());   // assign a name
      mppc[iMppc][jMppc]->SetLabel(mppc_label[mppcCounter]);
      mppc[iMppc][jMppc]->SetID(mppcCounter);             // assign an ID
      mppc[iMppc][jMppc]->SetI(iMppc); 
      mppc[iMppc][jMppc]->SetJ(jMppc);
      mppc[iMppc][jMppc]->SetChildrenI(ncrystalsx);
      mppc[iMppc][jMppc]->SetChildrenJ(ncrystalsy);
      mppc[iMppc][jMppc]->SetPosition(xPositions[mppcCounter],yPositions[mppcCounter],0);
      mppc[iMppc][jMppc]->SetDigitizerChannel(translateCh[mppcCounter]);
      mppc[iMppc][jMppc]->SetCanvasPosition(plotPositions[mppcCounter]);
      mppc[iMppc][jMppc]->SetParentName(module[iMppc/nmppcx][jMppc/nmppcy]->GetName());
      mppcCounter++;
    }
  }
  
  int crystalCounter = 0;
  Crystal*** crystal = new Crystal**[nmodulex*nmppcx*ncrystalsx]; // make an array of module pointers
  for(int j = 0; j < nmodulex*nmppcx*ncrystalsx ; j++)
  {
    crystal[j] = new Crystal* [nmoduley*nmppcy*ncrystalsy];
  }
  for(int iCrystal = 0; iCrystal < ncrystalsx*nmppcx*nmodulex ; iCrystal++)
  {
    for(int jCrystal = 0; jCrystal < ncrystalsy*nmppcy*nmoduley ; jCrystal++)
    {
      std::stringstream stream;
      stream << "Crystal " << crystalCounter;
      crystal[iCrystal][jCrystal] = new Crystal();
      crystal[iCrystal][jCrystal]->SetName(stream.str().c_str());   // assign a name
      crystal[iCrystal][jCrystal]->SetID(crystalCounter);          // assign an ID
      crystal[iCrystal][jCrystal]->SetI(iCrystal); 
      crystal[iCrystal][jCrystal]->SetJ(jCrystal);
      crystal[iCrystal][jCrystal]->SetParentName(mppc[iCrystal/ncrystalsx][jCrystal/ncrystalsy]->GetName());
      crystalCounter++;
    }
  }
  
    
  
  
  
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  module[iModule][jModule]->AddChild( mppc[(iModule * nmppcx) + iMppc][(jModule * nmppcy) + jMppc]->GetName() );
	  for(int iCrystal = 0; iCrystal < ncrystalsx ; iCrystal++)
	  {
	    for(int jCrystal = 0; jCrystal < ncrystalsy ; jCrystal++)
	    {
	      mppc[iMppc][jMppc]->AddChild( crystal[(iMppc * ncrystalsx) + iCrystal][(jMppc * ncrystalsy) + jCrystal]->GetName() );
	    }
	  }
	}
      }
    } 
  }
  
  
  
//   std::cout << "--------------------------------------------" << std::endl;
//   module[0][0]->Print();
//   std::cout <<module[0][0]->GetMppcsNumber() << std::endl;
//   for(int iModule = 0; iModule < nmodulex ; iModule++)
//   {
//     for(int jModule = 0; jModule < nmoduley ; jModule++)
//     {
//       module[iModule][jModule]->Print();
//       for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
//       {
// 	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
// 	{
// 	  mppc[(iModule * nmppcx) + iMppc][(jModule * nmppcy) + jMppc]->Print();
// 	  for(int iCrystal = 0; iCrystal < ncrystalsx ; iCrystal++)
// 	  {
// 	    for(int jCrystal = 0; jCrystal < ncrystalsy ; jCrystal++)
// 	    {
// 	      crystal[(iMppc * ncrystalsx) + iCrystal][(jMppc * ncrystalsy) + jCrystal]->Print();
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
// //   std::cout <<module[0][0]->GetMppcsNumber() << std::endl;
// //   std::cout << "--------------------------------------------" << std::endl;
// //   std::cout << std::endl;
//   
//   // draw a map
//   
//   // MPPC labels
//   
//   
  for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
  {
//     std::cout << "|---"
    for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
    {
      
      std::cout << mppc[iMppc][jMppc]->GetLabel() << "\t";
    }  
    std::cout << std::endl;
  }
  
  for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
  {
    for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
    {
      std::cout << mppc[iMppc][jMppc]->GetDigitizerChannel() << "\t";
    }  
    std::cout << std::endl;
  }
  
  TCanvas* c2 = new TCanvas("c2","c2",1200,800);
  c2->Divide(4,4);
  
//   TFile* fPlots = new TFile("plots.root","recreate");
  
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      TH2F* spectrum2d = new TH2F("spectrum2d","spectrum2d",150,-7,7,150,-7,7);
      tree->Draw("FloodY:FloodX >> spectrum2d","","COLZ");
      spectrum2d->SetName("Flood Histogram");
      module[iModule][jModule]->SetFloodMap2D(*spectrum2d);
      delete spectrum2d;
      
      TH2F *spherical = new TH2F("spherical","spherical",200,1.37,1.5,200,0.6,1.1); 
      tree->Draw("Phi:Theta >> spherical","FloodX > -7 && FloodX < 7 && FloodY > -7 && FloodY < 7 && FloodZ > 0 && FloodZ < 1","COLZ");
      
      module[iModule][jModule]->SphericalMap = *spherical;
      delete spherical;
      
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  TH1F* spectrum = new TH1F("spectrum","spectrum",2048,1,5000);
	  int channel    = mppc[iMppc][jMppc]->GetDigitizerChannel();
	  std::stringstream var,cut;
	  var << "ch" << channel << " >> spectrum";
	  cut << "TriggerChannel == " << channel  ;
	  c2->cd(mppc[iMppc][jMppc]->GetCanvasPosition());
	  tree->Draw(var.str().c_str(),cut.str().c_str());
	  TString name = "Spectrum - " + mppc[iMppc][jMppc]->GetLabel();
	  spectrum->SetName(name);
	  spectrum->SetTitle(name);
	  mppc[iMppc][jMppc]->SetRawSpectrum(*spectrum);
// 	  delete spectrum;
	}
      }
    }
  }
  
  
  //temporary save 
  //-----------------------------
  TFile* fPlots = new TFile("plots.root","recreate");
  fPlots->cd();
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      module[iModule][jModule]->GetFloodMap2D()->Write();
      (&(module[iModule][jModule]->SphericalMap))->Write();
      c2->Write();
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  mppc[iMppc][jMppc]->GetRawSpectrum()->Write();
	}
      }
    }
  }
  
  fPlots->Close();
  //----------------------------
  
  delete crystal;
  delete mppc;
  delete module;
  
  return 0;
}