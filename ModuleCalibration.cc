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
  // and prepare the TChain
//   TChain *chain =  new TChain("adc");
  // and assume as default that there is no config file name from command line 
  //then check
  if(std::string(argv[1]) == std::string("-c")) // first argument is -c, then the config file name is passed by command line
  {
    ConfigFileName = argv[2];
    std::cout << "Configuration file: '" << argv[2] << "'"<< std::endl;
//     for (int i = 3; i < argc ; i++) // run on the remaining arguments to add all the input files
//     {
//       std::cout << "Adding file " << argv[i] << std::endl;
//       chain->Add(argv[i]);
//     }
  }
  else // the config file was indeed the default one
  {
    std::cout << "Configuration file set to default: config.cfg "<< std::endl;
//     for (int i = 1; i < argc ; i++) // run on the remaining arguments to add all the input files
//     {
//       std::cout << "Adding file " << argv[i] << std::endl;
//       chain->Add(argv[i]);
//     }
  }
  
  
  //open the config file
  //TODO modify this in sucha  way that when calling 
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
  //--------------------------------------------------------------------------------------------------------------
  
  
  
  std::string chainName = config.read<std::string>("chainName");
  InputFile input(argc,argv,chainName,digitizer.size()); // read the input chain of root files, produces the ttree that will be used in the analysis

  
  //--------------------------------------------------------------------------------------------------------------
  int ncrystalsx = config.read<int>("ncrystalsx");
  int ncrystalsy = config.read<int>("ncrystalsy");
  int nmppcx     = config.read<int>("nmppcx");
  int nmppcy     = config.read<int>("nmppcy");
  int nmodulex   = config.read<int>("nmodulex");
  int nmoduley   = config.read<int>("nmoduley");
  
  //TEST
  
  
  
  // Module(s)
  Module*** module = new Module** [nmodulex];
  for(int j = 0; j < nmodulex ; j++)
  {
    module[j] = new Module* [nmoduley];
  }
  // MPPCs
  Mppc*** mppc = new Mppc** [nmppcx];
  for(int j = 0; j < nmppcx ; j++)
  {
    mppc[j] = new Mppc* [nmppcy];
  }
  // Crystals
  Crystal*** crystal = new Crystal** [ncrystalsx];
  for(int j = 0; j < ncrystalsx ; j++)
  {
    crystal[j] = new Crystal* [ncrystalsy];
  }
  
  //counters
  int moduleCounter = 0;
  int mppcCounter = 0;
  int crystalCounter = 0;
  
  
//   for(int iModule = 0; iModule < nmodulex ; iModule++)
//   {
//     for(int jModule = 0; jModule < nmoduley ; jModule++)
//     {
//       std::stringstream sname;
//       sname << "Module " << iModule << "." << jModule;
//       module[iModule][jModule] = new Module();                 // creates a default module
//       module[iModule][jModule]->SetName(sname.str().c_str());  // assign a name
//       module[iModule][jModule]->SetID(moduleCounter);          // assign an ID
//       module[iModule][jModule]->SetI(iModule); 
//       module[iModule][jModule]->SetJ(jModule);
//       moduleCounter++;
//     }
//   }
//   
//   for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
//   {
//     for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
//     {
//       std::stringstream sname;
//       sname << "Mppc " << iMppc << "." << jMppc;
//       mppc[iMppc][jMppc] = new Mppc();                 // creates a default mppc
//       mppc[iMppc][jMppc]->SetName(sname.str().c_str());   // assign a name
//       mppc[iMppc][jMppc]->SetID(mppcCounter);          // assign an ID
//       mppc[iMppc][jMppc]->SetI(iMppc); 
//       mppc[iMppc][jMppc]->SetJ(jMppc); 
//       mppcCounter++;
//     }
//   }
//   
//   for(int iCrystal = 0; iCrystal < ncrystalsx ; iCrystal++)
//   {
//     for(int jCrystal = 0; jCrystal < ncrystalsy ; jCrystal++)
//     {
//       std::stringstream stream;
//       stream << "Crystal " << crystalCounter;
//       crystal[iCrystal][jCrystal] = new Crystal();                 // creates a default crystal
//       crystal[iCrystal][jCrystal]->SetName(stream.str().c_str());   // assign a name
//       crystal[iCrystal][jCrystal]->SetID(crystalCounter);          // assign an ID
//       crystal[iCrystal][jCrystal]->SetI(iCrystal); 
//       crystal[iCrystal][jCrystal]->SetJ(jCrystal);
//       crystalCounter++;
//     }
//   }
//   
//   for(int iModule = 0; iModule < (nmodulex*nmoduley) ; iModule++) // run on modules
//   {
//     //no assignment for the modules, they have no parent 
//     for(int iMppc = 0 ; iMppc < (nmppcx*nmppcy) ; iMppc++) // run on mppcs
//     {
//       mppc[int (iMppc / nmppcx)][iMppc % nmppcy]->SetModule(module[iModule / nmodulex][jModule % nmoduley]);
//       for(int iCrystal = 0; iCrystal < ncrystalsx*ncrystalsy ; iCrystal++) // run on crystals
//       {
// 	
//       }
//     }
//   }
  
  
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      std::stringstream sname;
      std::string       moduleName;
      sname << "Module " << iModule << "." << jModule;
      moduleName = sname.str();
      module[iModule][jModule] = new Module();                 // creates a default module
      module[iModule][jModule]->SetName(moduleName.c_str());  // assign a name
      module[iModule][jModule]->SetID(moduleCounter);          // assign an ID
      module[iModule][jModule]->SetI(iModule); 
      module[iModule][jModule]->SetJ(jModule);
      module[iModule][jModule]->MakeChildrenPointers(nmppcx,nmppcy); //prepare the pointers for the children
      
      mppcCounter = 0;
      
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  std::stringstream sname;
	  std::string       mppcName;
	  sname << moduleName << " - Mppc " << iMppc << "." << jMppc;
	  mppcName = sname.str();
	  mppc[iMppc][jMppc] = new Mppc();                 // creates a default mppc
	  mppc[iMppc][jMppc]->SetName(mppcName.c_str());   // assign a name
	  mppc[iMppc][jMppc]->SetID(mppcCounter);          // assign an ID
	  mppc[iMppc][jMppc]->SetI(iMppc); 
	  mppc[iMppc][jMppc]->SetJ(jMppc); 
	  mppc[iMppc][jMppc]->MakeChildrenPointers(ncrystalsx,ncrystalsy);
	  mppc[iMppc][jMppc]->SetModule(module[iModule][jModule]); //assign to this mppc its parent module
	  module[iModule][jModule]->SetChild(iMppc,jMppc,mppc[iMppc][jMppc]);  //assign this mppc to its parent module
	  
	  crystalCounter = 0;
	  for(int iCrystal = 0; iCrystal < ncrystalsx ; iCrystal++)
	  {
	    for(int jCrystal = 0; jCrystal < ncrystalsy ; jCrystal++)
	    {
	      std::stringstream sname;
	      std::string       crystalName;
	      sname << mppcName << " - Crystal " << crystalCounter;
	      crystalName = sname.str();
	      crystal[iCrystal][jCrystal] = new Crystal();                 // creates a default crystal
	      crystal[iCrystal][jCrystal]->SetName(crystalName.c_str());   // assign a name
	      crystal[iCrystal][jCrystal]->SetID(crystalCounter);          // assign an ID
	      crystal[iCrystal][jCrystal]->SetI(iCrystal); 
	      crystal[iCrystal][jCrystal]->SetJ(jCrystal); 
	      crystal[iCrystal][jCrystal]->SetMppc(mppc[iMppc][jMppc]);
	      mppc[iMppc][jMppc]->SetChild(iCrystal,jCrystal,crystal[iCrystal][jCrystal]);
	      
	      crystalCounter++;
	      sname.str("");
	    }
	  }
	  
	  mppcCounter++;
	}
      }
      
      moduleCounter++;
    }
  }
  
  
  
  
//   for(int i = 0; i < ncrystalsx ; i++)
//   {
//     for(int j = 0; j < ncrystalsy ; j++)
//     {
//       std::cout << "--------------------------------------------" << std::endl;
//       crystal[i][j]->Print();
//       std::cout << "--------------------------------------------" << std::endl;
//     }
//   }
  std::cout << std::endl;
  
  
  
  
//   for(int i = 0; i < nmppcx ; i++)
//   {
//     for(int j = 0; j < nmppcy ; j++)
//     {
//       std::cout << "--------------------------------------------" << std::endl;
//       mppc[i][j]->Print();
//       std::cout << "--------------------------------------------" << std::endl;
//     }
//   }
  
  std::cout << std::endl;
  
  
//   for(int i = 0; i < nmodulex ; i++)
//   {
//     for(int j = 0; j < nmoduley ; j++)
//     {
//       std::cout << "--------------------------------------------" << std::endl;
//       module[i][j]->Print();
//       std::cout << "--------------------------------------------" << std::endl;
//     }
//   }
  
  std::cout << std::endl;
  //--------------------------------------------------------------------------------------------------------------
  
  
  
  // build the structure, i.e. assign crystals to mppcs and mppcs to modules
  //TEST
  //crystal[0][0]->SetMppc(mppc[0][0]);
//   crystal[0][0]->GetMppc()->Print();
//   crystal[0][0]->Print();
  
  
  module[0][0]->GetChild(2,3)->GetChild(1,0)->Print();
  
  
  
  return 0;
}