// compile with
// g++ -o ../build/shiftTree shiftTree.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

// small program to extract timing calibration and data

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

#include "./libraries/CrystalStructs.h"
// #include "./include/ConfigFile.h"


// typedef std::vector<std::string> stringvec;
// list files in directory
// taken from
// http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
// void read_directory(const std::string& name, std::vector<std::string> &v)
// {
//     DIR* dirp = opendir(name.c_str());
//     struct dirent * dp;
//     while ((dp = readdir(dirp)) != NULL) {
//         v.push_back(dp->d_name);
//     }
//     closedir(dirp);
// }

void usage()
{
  std::cout << "\t\t" << "[-i|--input] <input_file> [-o|--output] <output_file> [-c|--config] <config_file>" << std::endl
            << "\t\t\t" << "<input_file>                                      - name of input file " << std::endl
            << "\t\t\t" << "<output_file>                                     - name of output file "   << std::endl
            << "\t\t\t" << "<config_file>                                     - name of config file "   << std::endl
            << std::endl;
}

//----------------//
//  MAIN PROGRAM  //
//----------------//
int main (int argc, char** argv)
{

  // check if there are args, otherwise print the usage info
  if(argc < 2)
  {
    std::cout << argv[0];
    usage();
    return 1;
  }


  // save the entire command line
  std::stringstream streamCommand;
  for(int i=0 ; i < argc; i++)
  {
    streamCommand << argv[i] << " ";
  }

  // set stat and fit information level in root files
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);


  //---------------------------------------//
  // PARSE COMMAND LINE INPUTs             //
  //---------------------------------------//

  std::string inputFileName = "";
  std::string outputFileName = "";
  std::string configFileName = "";

  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "config", required_argument, 0, 0 },
      { NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			inputFileName = (char *)optarg;
    }
    if (c == 'o'){
			outputFileName = (char *)optarg;
    }
    if (c == 'c'){
			configFileName = (char *)optarg;
    }
		else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      configFileName = (char *)optarg;
    }
    else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
  }


  // check if required are given and files actually exists
  // first, input given and not empty
  if(inputFileName == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide the input file name!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  if(outputFileName == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide the output file name!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  if(configFileName == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide the config output file name!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }

  std::cout << "input file = "  << inputFileName << std::endl
            << "ouptut file = " << outputFileName << std::endl
            << "config file = " << configFileName << std::endl;
  //READ SHIFT CONFIG file
  std::ifstream configFile;
  configFile.open(configFileName.c_str());
  std::vector<int> original;
  std::vector<int> target;
  if (configFile.is_open())
  {
    while(!configFile.eof())
    {
      int a,b;
      configFile >> a >> b;
      if(!configFile.eof())
      {
        original.push_back(a);
        target.push_back(b);
      }
    }
  }
  else
  {
    std::cout << std::endl;
    std::cout << "ERROR! cannot open " << configFileName << " file! Aborting..." << std::endl;
    return 1;
  }
  configFile.close();



  //INPUT TChain
  //find detector channels
  TChain* tree = new TChain("adc");  // create the input tchain and the analysis ttree
  tree->Add(inputFileName.c_str());
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
  // int numOfCry = 0;
  std::string ch_prefix("ch");
  std::string t_prefix("t");
  // std::string cry_prefix("cry");
  for(int i = 0 ; i < nLeaves ; i++)
  {
    if (!leavesName[i].compare(0, ch_prefix.size(), ch_prefix))
    {
      numOfCh++;
      detector_channels.push_back(atoi( (leavesName[i].substr(ch_prefix.size(),leavesName[i].size()-ch_prefix.size())).c_str() )) ;
    }
  }
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;

  //Prepare branches
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

  //OUTPUT TTree
  //declare ROOT ouput TTree and file
  ULong64_t out_DeltaTimeTag = 0;
  ULong64_t out_ExtendedTimeTag = 0;
  ULong64_t out_startTimeTag = 0;
  UShort_t *out_charge;
  Float_t *out_timestamp;
  //the ttree variable
  TTree *t1 ;
  //strings for the names
  std::stringstream snames;
  std::stringstream stypes;
  std::string names;
  std::string types;

  out_charge = new UShort_t[numOfCh];
  out_timestamp = new Float_t[numOfCh];

  t1 = new TTree("adc","adc");
  t1->Branch("ExtendedTimeTag",&out_ExtendedTimeTag,"ExtendedTimeTag/l");   //absolute time tag of the event
  t1->Branch("DeltaTimeTag",&out_DeltaTimeTag,"DeltaTimeTag/l");                    //delta time from previous event
  for (int i = 0 ; i < detector_channels.size() ; i++)
  {
    //empty the stringstreams
    snames.str(std::string());
    stypes.str(std::string());
    out_charge[i] = 0;
    snames << "ch" << i;
    stypes << "ch" << i << "/s";
    names = snames.str();
    types = stypes.str();
    t1->Branch(names.c_str(),&out_charge[i],types.c_str());
  }
  for (int i = 0 ; i < detector_channels.size() ; i++)
  {
    //empty the stringstreams
    snames.str(std::string());
    stypes.str(std::string());
    out_timestamp[i] = 0;
    snames << "t" << i;
    stypes << "t" << i << "/F";
    names = snames.str();
    types = stypes.str();
    t1->Branch(names.c_str(),&out_timestamp[i],types.c_str());
  }


  if((numOfCh != original.size()) || (numOfCh != target.size()))
  {
    std::cout << "ERROR! Mismatch in TTree and config file length!" << std::endl;
    return 1;
  }
  else
  {
    std::cout << "numOfCh         = " << numOfCh         << std::endl
              << "original.size() = " << original.size() << std::endl
              << "target.size()   = " << target.size()   << std::endl;
  }

  //shift map
  for(int iTar = 0 ; iTar < target.size() ; iTar++)
  {
    std::cout << "original " << original[iTar] << " to " << target[iTar] << std::endl;
  }


  //LOOP
  long long int nevent = tree->GetEntries();
  long long int counter = 0;
  std::cout << "Total number of events in analysis file = " << nevent << std::endl;

  for (long long int i=0;i<nevent;i++)
  {
    // std::cout << "Event " << i << std::endl;
    tree->GetEvent(i);              //read complete accepted event in memory
    // copy timetags
    out_ExtendedTimeTag = ChainExtendedTimeTag;
    out_DeltaTimeTag = ChainDeltaTimeTag;
    // shift data
    for(int iTar = 0 ; iTar < target.size() ; iTar++)
    {
      out_charge[target[iTar]] = charge[original[iTar]];
      out_timestamp[target[iTar]]= timeStamp[original[iTar]];
    }

    t1->Fill();

    counter++;
    int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
  }
  std::cout << std::endl;

  std::cout << "Writing output TTree to " << outputFileName.c_str() << std::endl;
  TFile* fTree = new TFile(outputFileName.c_str(),"recreate");
  fTree->cd();
  t1->Write();
  fTree->Close();
  return 0;
}
