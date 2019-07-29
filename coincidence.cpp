// template program to import n calibration data and loop on m acquisition ttrees
// use it as a base for analysis programs

// compile with
// g++ -o ../build/coincidence coincidence.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas


#include "TROOT.h"
#include "TFile.h"
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
#include "TVector.h"
#include "TNamed.h"
#include "TPaveLabel.h"
#include "THStack.h"
#include "TFitResult.h"
#include "TMatrixD.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <getopt.h>
#include <algorithm>    // std::sort
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "./libraries/CrystalStructs.h"     // Crystal_t , detector_t
#include "./libraries/Calibration.h"        // readTaggingData , readCalibration , setWandZcuts
#include "./libraries/Utilities.h"          // read_directory , invert_a_matrix
#include "./libraries/Extract.h"            // extractWithEMG , extractCTR , FindSmallestInterval

// forward declaration of usage info output
void usage();


//----------------//
//  MAIN PROGRAM  //
//----------------//
int main (int argc, char** argv)
{
  if(argc < 2)
  {
    std::cout	<< "Usage: " << std::endl << argv[0] ;
    usage();
    return 1;
  }

  std::stringstream streamCommand;
  for(int i=0 ; i < argc; i++)
  {
    streamCommand << argv[i] << " ";
  }

  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  ROOT::v5::TFormula::SetMaxima(10000,10000,10000);

  // default args
  std::string calibrationFileNames = "";
  std::string inputFolderName = "./";
  std::string inputFilePrefix = "";
  std::string outputFileName = "outputBareboneFile.root";

  // parse command line arguments
  static struct option longOptions[] =
  {
      { "calibration", required_argument, 0, 0 },
      { "folder", required_argument, 0, 0 },
      { "prefix", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "c:f:p:o:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'c'){
			calibrationFileNames = (char *)optarg;
    }
		else if (c == 'f'){
      inputFolderName = (char *)optarg;
    }
    else if (c == 'p'){
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 'o'){
      outputFileName = (char *)optarg;
    }

		else if (c == 0 && optionIndex == 0){
      calibrationFileNames = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      inputFolderName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 3){
      outputFileName = (char *)optarg;
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}



  // check if required are given and files actually exists
  // first, input given and not empty
  if(inputFilePrefix == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide the prefix of input files!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  if(calibrationFileNames == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide calibration files!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }



  //----------------------------------//
  // GET CALIBRATION(S)               //
  //----------------------------------//
  //calibration files
  std::vector<std::string> listCalibrationFiles;
  split( listCalibrationFiles, calibrationFileNames, "," );  // split the entry
  bool calibrationFilesExist = true;
  for(unsigned int i = 0 ; i < listCalibrationFiles.size() ; i++)
  {
    if(!fileExists(listCalibrationFiles[i]))
    {
      calibrationFilesExist = false;
      std::cout << "ERROR! File " << listCalibrationFiles[i] << " does NOT exist!!!" << std::endl;
    }
  }

  if(calibrationFilesExist == false)
  {
    std::cout << std::endl;
    std::cout << "ERROR! Some input files do not exists! Aborting." << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  //open the calibration files
  std::vector<TFile*> calibrationFile;
  for(unsigned int i = 0 ; i < listCalibrationFiles.size() ; i++)
  {
    TFile* pCalibrationFile = new TFile(listCalibrationFiles[i].c_str());
    calibrationFile.push_back(pCalibrationFile);
  }





  //---------------------------------------//
  // FEEDBACK PARAMETERS                   //
  //---------------------------------------//
  std::cout << std::endl;
  std::cout << "//-------------------------------------//"  << std::endl;
  std::cout << "// INPUT PARAMETERS                    //"  << std::endl;
  std::cout << "//-------------------------------------//"  << std::endl;
  std::cout << "Input folder             = " << inputFolderName   << std::endl;
  std::cout << "Input file prefix        = " << inputFilePrefix   << std::endl;
  std::cout << "Calibration files        = " ;
  for(unsigned int i = 0 ; i < listCalibrationFiles.size() ; i++)
  {
    std::cout << listCalibrationFiles[i];
    if(i < (listCalibrationFiles.size() -1))
    {
      std::cout << ",";
    }
  }


  //----------------------------------//
  // GET INPUT FILES(S)               //
  //----------------------------------//
  // read file in dir
  std::cout << std::endl;
  std::cout << "|----------------------------------------|" << std::endl;
  std::cout << "|         ANALYSIS FILES                 |" << std::endl;
  std::cout << "|----------------------------------------|" << std::endl;
  std::cout << std::endl;
  // get input files list
  std::vector<std::string> v;
  read_directory(inputFolderName, v);
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;
  for(unsigned int i = 0 ; i < v.size() ; i++)
  {
    if(!v[i].compare(0,inputFilePrefix.size(),inputFilePrefix))
    {
      listInputFiles.push_back(inputFolderName + "/" + v[i]);
    }
  }
  // check if it's empty
  if(listInputFiles.size() == 0)
  {
    std::cout << std::endl;
    std::cout << "ERROR! Some input files do not exists! Aborting." << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }


  //----------------------------------------------------------//
  //  Get TChain from input TTree files                       //
  //----------------------------------------------------------//
  TChain* tree = new TChain("adc");  // create the input tchain and the analysis ttree
  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    std::cout << "Adding file " << listInputFiles[i] << std::endl;
    tree->Add(listInputFiles[i].c_str());
  }
  std::cout << "|----------------------------------------|" << std::endl;
  std::cout << std::endl;
  std::vector<int> detector_channels;
  TObjArray *leavescopy = tree->GetListOfLeaves();
  int nLeaves = leavescopy->GetEntries();
  std::vector<std::string> leavesName;
  // fill a vector with the leaves names
  for(int i = 0 ; i < nLeaves ; i++)
  {
    leavesName.push_back(leavescopy->At(i)->GetName());
  }
  // count the entries that start with "ch"
  int numOfCh = 0;
  std::string ch_prefix("ch");
  std::string t_prefix("t");
  for(int i = 0 ; i < nLeaves ; i++)
  {
    if (!leavesName[i].compare(0, ch_prefix.size(), ch_prefix))
    {
      numOfCh++;
      detector_channels.push_back(atoi( (leavesName[i].substr(ch_prefix.size(),leavesName[i].size()-ch_prefix.size())).c_str() )) ;
    }
  }
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;

  //set variables and branches
  ULong64_t     ChainExtendedTimeTag;                                // extended time tag
  ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous
  UShort_t      *charge;
  Float_t      *timeStamp;
  TBranch      *bChainExtendedTimeTag;                               // branches for above data
  TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  TBranch      **bCharge;
  TBranch      **btimeStamp;
  charge = new UShort_t[numOfCh];
  timeStamp = new Float_t[numOfCh];
  bCharge = new TBranch*[numOfCh];
  btimeStamp = new TBranch*[numOfCh];
  // set branches for reading the input files
  tree->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag, &bChainExtendedTimeTag);
  tree->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag, &bChainDeltaTimeTag);
  for (int i = 0 ; i < detector_channels.size() ; i++)
  {
    //empty the stringstreams
    std::stringstream sname;
    sname << "ch" << detector_channels[i];
    tree->SetBranchAddress(sname.str().c_str(),&charge[detector_channels[i]],&bCharge[detector_channels[i]]);
    sname.str("");
    sname << "t" << detector_channels[i];
    tree->SetBranchAddress(sname.str().c_str(),&timeStamp[detector_channels[i]],&btimeStamp[detector_channels[i]]);
    sname.str("");
  }


  TList *formulasAnalysis = new TList();
  std::vector<Crystal_t> crystal;

  for(unsigned int i = 0 ; i < calibrationFile.size() ; i++)
  {
    readCalibration(calibrationFile[i],       // this calib file
                    tree,                     // input TChain (same for everyone)
                    formulasAnalysis,         // TList of all TTreeFormula
                    crystal);                 // structure of all crystals found in all calib lifes


  }
  // optionally set w and z limits, and write values into crystal struct
  // setWandZcuts(crystal);


  // list the crystals with calibration data found
  std::cout << "Calibration data found for crystals: " << std::endl;
  for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  {
    if(crystal[i].accepted)
    {
      std::cout << crystal[i].number << std::endl;
    }
  }




  // MAIN LOOP
  long long int counter = 0;
  tree->SetNotify(formulasAnalysis);
  long long int neventAnalysis = tree->GetEntries();
  std::cout << "Total number of events in analysis file = " << neventAnalysis << std::endl;
  long int goodEventsAnalysis = 0;
  long int counterAnalysis = 0;
  for (long long int i=0;i<neventAnalysis;i++)
  {

    tree->GetEvent(i);              //read complete accepted event in memory

    // LOOP ON ALL ACCEPTED CRYSTALS
    // I.E. ALL CRYSTALS IN THE CALIBRATION FILES
    int c = 0;
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {

      // do whatever you want
      // in particular, some tips:
      // 1) if(crystal[iCry].FormulaTagAnalysis->EvalInstance()){} --> this will select events in the photopeak of the external reference crystal, if there is one. Notice that each crystal can have its own FormulaTagAnalysis, since they can come from different calibration runs
      // 2) if(crystal[iCry].FormulaAnalysis->EvalInstance()) --> this will select events in a crystal (geometrical position in u-v-w and photopeak condition)
      // 3) to calculate w, if needed
      // float FloodZ = calculateFloodZ(charge,crystal[iCry]);
      // notice that the detectorSaturation values are written directly into the crystal struct, because they can be acquired in different conditions (but it has to be the same photodetector!)

      //example here, selecting photopeak events in ext ref, photopeak events the crystal(s) found in calibration file(s), and calculating FloodZ for each
      if(crystal[iCry].accepted)
      {

        // if(crystal[iCry].FormulaAnalysis)
        // {
          if(crystal[iCry].FormulaAnalysis->EvalInstance())
          {
            // std::cout << i << " " << std::endl;
            // c++;

            // std::cout << crystal[iCry].number << " ";
            for(unsigned int jCry = 0 ;  jCry < crystal.size() ; jCry++)
            {
              if(crystal[jCry].number != crystal[iCry].number )
              {
                if(crystal[jCry].accepted)
                {
            //       // if(crystal[jCry].FormulaAnalysis)
            //       // {
                    if(crystal[jCry].FormulaAnalysis->EvalInstance())
                    {
                      // std::cout << crystal[iCry].timingChannel << " " << crystal[jCry].timingChannel << std::endl;
                      if(fabs( timeStamp[crystal[iCry].detectorChannel] - timeStamp[crystal[jCry].detectorChannel] ) < 1e-9) // FIXME this works only because channels are the same N for charge and t, and it is used here because the calibration files were obtained wihtout timing analisys
                      {
                        goodEventsAnalysis++;
                      }
                      // std::cout << crystal[iCry].number << " "
                                // << crystal[jCry].number << std::endl;
            //

            //           // std::cout << "ciao" << std::endl;
                    }
            //       // }
            //
                }
              }
            }
          }
        // }

      }
    }
    // if(c > 1)
    // {
      // std::cout << i << std::endl;
    // }



    // LOOP COUNTER
    counterAnalysis++;
    int perc = ((100*counterAnalysis)/neventAnalysis);
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }

  std::cout << "Good events = " << goodEventsAnalysis << std::endl;

  // sort crystals struct (can be useful)
  std::sort(crystal.begin(), crystal.end(), compareByNumber);

  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  // write whatever you want to save

  //save the full command line
  TNamed CommandNameD("Command",streamCommand.str().c_str());
  CommandNameD.Write();
  outputFile->Close();
  std::cout << "Results saved in file " << outputFileName << std::endl;

  return 0;
}
// end of main program




// feedback to user
void usage()
{
  std::cout << "\t"
            << "[-c|--calibration] <list> "
            << "\t"
            << "[-f|--folder] <folder> "
            << "\t"
            << "[-p|--prefix] <prefix> "
            << "\t"
            << "[-o|--output] <output> "
            << std::endl
            << "\t\t"
            << "<list>           - csv list of calibration.root files (outputs of ModuleCalibration) "
            << std::endl
            << "\t\t"
            << "<folder>         - path to folder were input files are located - default = ./ "
            << std::endl
            << "\t\t"
            << "<prefix>         - prefix of input TTree files "
            << std::endl
            << "\t\t"
            << "<output>         - output file name - default = outputBareboneFile.root "
            << std::endl
            << std::endl;
}
