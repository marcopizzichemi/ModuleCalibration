//------------------------------------------------------------//
//                                                            //
//  PROGRAM FOR matrix calibration USING DT5740               //
//                                                            //
//------------------------------------------------------------//

// compile with 
// g++ -o ModuleCalibration ModuleCalibration.cpp `root-config --cflags --glibs` -lSpectrum -lMLP -lTreePlayer

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

#include "ConfigFile.h"



int main (int argc, char** argv)
{
  
  //----------------------------------------------------------//
  //  Check input                                             //
  //----------------------------------------------------------//
  if(argc<2) 
  {
    std::cout << " Usage: " 							<< std::endl;
    std::cout << " ModuleCalibration [-c config file ] <input-files> " 		<< std::endl;
    std::cout << "   note: -c option is optional" 				<< std::endl;
    std::cout << "         without it, a configFile.txt will be assumed" 	<< std::endl;
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
  //  Parse the config file                                   //
  //----------------------------------------------------------//
  // Set a default config file name
  std::string ConfigFileName = "config.cfg"; 
  // and assume there is no config input from command line 
  bool configIsFromCommandLine = false;
  // if passed by command line, override the config file name
  for (int i = 0; i < argc ; i++)
  {
    if(std::string(argv[i]) == std::string("-c")) 
    {
      ConfigFileName = argv[i+1];
      std::cout << "Configuration file: '" << argv[i+1] << "'"<< std::endl;
      configIsFromCommandLine = true;
    }
  }
  if(!configIsFromCommandLine)
  {
    std::cout << "Configuration file set to default: config.cfg "<< std::endl;
  }
  //open the config file
  ConfigFile config(ConfigFileName);
  
  //read crystal dimensions
  std::string dummy = config.read<std::string>("dummy");
  int d1 = config.read<int>("ddd");
  int d2 = config.read<int>("ddd");
  std::cout << "Dummy: " << dummy << " " << d1 << " " << d2 << " "  << std::endl;
  
  
  
  return 0;
}