// compile with
// g++ -o ../../build/filterForDeepLearning_MiniPET_old filterForDeepLearning_MiniPET_old.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

// small program to take doi bench data and filter events in one cyrstal, producing an output for deeplearning study
// this version is for input from MINI pet electronics

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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <getopt.h>
#include <algorithm>    // std::sort
#include <sys/types.h>

#include "../libraries/CrystalStructs.h"        // Crystal_t , detector_t
#include "../libraries/Calibration.h"           // readTaggingData , readCalibration , setWandZcuts
#include "../libraries/Utilities.h"             // read_directory , invert_a_matrix





void usage()
{
  std::cout << "\t\t" << "[-i|--input] <file_prefix>  [-o|--output] <output.txt> " << std::endl
            << "\t\t" << "<file_prefix>                                     - prefix of TTree files to analyze"   << std::endl
            << "\t\t" << "<output.txt>                                      - output file name"   << std::endl
            << "\t\t" << "<output.txt>                                      - output file name"   << std::endl
            << "\t\t" << std::endl;
}





int main (int argc, char** argv)
{
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }
  std::string inputFileName = "";
  std::string outputFileName = "";
  std::string calibrationFileName = "";
  int selectedCrystal = -1;
  bool verbose = false;
  bool UseAllCuts             = true;
  bool UseTaggingPhotopeakCut = false;
  bool UseTriggerChannelCut   = false;
  bool UseBroadCut            = false;
  bool UseCutNoise            = false;
  bool UsePhotopeakEnergyCut  = false;
  bool UseCutgs               = false;

  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "calibration", required_argument, 0, 0 },
      { "crystal", required_argument, 0, 0 },
      { "verbose", no_argument, 0, 0 },
      { "tagCut", no_argument, 0, 0 },
      { "triggerCut", no_argument, 0, 0 },
      { "broadCut", no_argument, 0, 0 },
      { "noiseCut", no_argument, 0, 0 },
      { "photopeakCut", no_argument, 0, 0 },
      { "cutgCut", no_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:c:n:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			inputFileName = (char *)optarg;
    }
		else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
    else if (c == 'c'){
      calibrationFileName = (char *)optarg;
    }
    else if (c == 'n'){
      selectedCrystal = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      calibrationFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 3){
      selectedCrystal = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      verbose = true;
    }
    else if (c == 0 && optionIndex == 5){
      UseAllCuts = false;
      UseTaggingPhotopeakCut = true;
    }
    else if (c == 0 && optionIndex == 6){
      UseAllCuts = false;
      UseTriggerChannelCut = true;
    }
    else if (c == 0 && optionIndex == 7){
      UseAllCuts = false;
      UseBroadCut = true;
    }
    else if (c == 0 && optionIndex == 8){
      UseAllCuts = false;
      UseCutNoise = true;
    }
    else if (c == 0 && optionIndex == 9){
      UseAllCuts = false;
      UsePhotopeakEnergyCut = true;
    }
    else if (c == 0 && optionIndex == 10){
      UseAllCuts = false;
      UseCutgs = true;
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}


  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");

  if(UseAllCuts)
  {
    UseTaggingPhotopeakCut = true;
    UseTriggerChannelCut   = true;
    UseBroadCut            = true;
    UseCutNoise            = true;
    UsePhotopeakEnergyCut  = true;
    UseCutgs               = true;
  }

  // read file in dir
  std::vector<std::string> v;
  read_directory(".", v);
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;
  for(unsigned int i = 0 ; i < v.size() ; i++)
  {
    if(!v[i].compare(0,inputFileName.size(),inputFileName))
    {
      listInputFiles.push_back(v[i]);
    }
  }



  //----------------------------------------------------------//
  //  Get TChain of input files                               //
  //----------------------------------------------------------//
  TChain* tree = new TChain("adc");  // create the input tchain and the analysis ttree
  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    std::cout << "Adding file " << listInputFiles[i] << std::endl;
    tree->Add(listInputFiles[i].c_str());
  }

  std::vector<int> detector_channels;
  std::vector<int> t_channels;
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
  int numOfT = 0;
  std::string t_prefix("t");

  for(int i = 0 ; i < nLeaves ; i++)
  {
    if (!leavesName[i].compare(0, ch_prefix.size(), ch_prefix))
    {
      numOfCh++;
      detector_channels.push_back(atoi( (leavesName[i].substr(ch_prefix.size(),leavesName[i].size()-ch_prefix.size())).c_str() )) ;
    }
  }

  for(int i = 0 ; i < nLeaves ; i++)
  {
    if (!leavesName[i].compare(0, t_prefix.size(), t_prefix))
    {
      numOfT++;
      t_channels.push_back(atoi( (leavesName[i].substr(t_prefix.size(),leavesName[i].size()-t_prefix.size())).c_str() )) ;
    }
  }

  std::cout << "Charge Channels \t= " << numOfCh << std::endl;
  std::cout << "Time Channels \t= " << numOfT << std::endl;

  ULong64_t     ChainExtendedTimeTag;                                // extended time tag
  ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous
  UShort_t     *charge;
  Float_t      *timeStamp;
  TBranch      *bChainExtendedTimeTag;                               // branches for above data
  TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  TBranch      **bCharge;
  TBranch      **btimeStamp;
  charge = new UShort_t[numOfCh];
  timeStamp = new Float_t[numOfT];
  bCharge = new TBranch*[numOfCh];
  btimeStamp = new TBranch*[numOfT];

  // set branches for reading the input files
  tree->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag, &bChainExtendedTimeTag);
  tree->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag, &bChainDeltaTimeTag);
  for (int i = 0 ; i < detector_channels.size() ; i++)
  {
    std::stringstream sname;
    sname << "ch" << detector_channels[i];
    tree->SetBranchAddress(sname.str().c_str(),&charge[detector_channels[i]],&bCharge[detector_channels[i]]);
    sname.str("");
  }

  for (int i = 0 ; i < t_channels.size() ; i++)
  {
    std::stringstream sname;
    sname << "t" << t_channels[i];
    tree->SetBranchAddress(sname.str().c_str(),&timeStamp[t_channels[i]],&btimeStamp[t_channels[i]]);
    sname.str("");
  }
  //----------------------------------------------------------//


  //----------------------------------------------------------//
  //  Set TTree for output                                    //
  //----------------------------------------------------------//
  //declare ROOT ouput TTree and file
  Float_t   out_doiFromTag = 0;
  ULong64_t out_DeltaTimeTag = 0;
  ULong64_t out_ExtendedTimeTag = 0;
  Float_t   out_charge[64];
  Float_t   out_timestamp[64];
  //the ttree variable
  TTree *t1 ;
  //strings for the names
  std::stringstream snames;
  std::stringstream stypes;
  std::string names;
  std::string types;

  t1 = new TTree("adc","adc");

  t1->Branch("ExtendedTimeTag",&out_ExtendedTimeTag,"ExtendedTimeTag/l");   //absolute time tag of the event
  t1->Branch("DeltaTimeTag",   &out_DeltaTimeTag,"DeltaTimeTag/l");                    //delta time from previouevent
  t1->Branch("zFromTag",&out_doiFromTag,"zFromTag/F");
  for (int i = 0 ; i < 64 ; i++)
  {
    //empty the stringstreams
    snames.str(std::string());
    stypes.str(std::string());
    out_charge[i] = 0.;
    snames << "ch" << i;
    stypes << "ch" << i << "/F";
    names = snames.str();
    types = stypes.str();
    t1->Branch(names.c_str(),&out_charge[i],types.c_str());
  }
  for (int i = 0 ; i < 64 ; i++)
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


  //----------------------------------------------------------//
  //  Read and import from Calibration file                   //
  //----------------------------------------------------------//
  // open calibration file
  TFile *calibrationFile = new TFile(calibrationFileName.c_str());
  // prepare a list of TTreeFormula
  TList* formulasAnalysis = new TList();
  // prepare a vector of crystals (we will use only one)
  std::vector<Crystal_t> crystal;

  readCalibration(calibrationFile,          // calib file
                  tree,                     // input TChain (same for everyone)
                  formulasAnalysis,         // TList of all TTreeFormula
                  crystal,                  // structure of all crystals found in all calib lifes
                  UseTriggerChannelCut,                      // include TriggerChannelCuthannel cut in crystalCut
                  UseBroadCut,                               // include broadCut in crystalCut
                  UseCutNoise,                               // include CutNoise in crystalCut
                  UsePhotopeakEnergyCut,                     // include PhotopeakEnergyCut in crystalCut
                  UseCutgs                                   // include CutGs in crystalCut
                 );

  // optionally set w and z limits, and write values into crystal struct
  // setWandZcuts(crystal);

  // list the crystals with calibration data found
  std::cout << "Calibration data found for crystals: " << std::endl;
  for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  {
    if(crystal[i].accepted)
    {
      std::cout << crystal[i].number << " "
                // << crystal[i].detectorChannel << " "
                << std::endl;
      // for(int iDet = 0 ; iDet < crystal[i].detectorSaturation.size(); iDet++)
      // {
      //   std::cout << crystal[i].detectorSaturation[iDet].digitizerChannel << " ";
      // }
      // std::cout << std::endl;
    }
  }


  // switch off crystals if they were not chosen (only one chosen allowed)
  // std::cout << "Calibration data found for crystals: " << std::endl;
  int crystalsFound = 0;
  int id = -1;
  for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  {
    if(crystal[i].accepted)
    {
      if(crystal[i].number != selectedCrystal)
      {
        crystal[i].accepted = false;
      }
      else
      {
        crystalsFound++;
        id = i;
      }
    }
  }

  if(crystalsFound == 0)
  {
    std::cout << "ERROR: No chosen crystal found? Aborting." << std::endl;
    return 1;
  }
  if(crystalsFound > 1)
  {
    std::cout << "ERROR: More than one chosen crystal found?? Aborting."<< std::endl << "crystalsFound = " << crystalsFound << std::endl;
    return 1;
  }
  if(crystalsFound == 1 )
  {
    if(id == -1 )
    {
      std::cout << "ERROR: This cannot happen. Aborting." << std::endl
                << "crystalsFound = " << crystalsFound << std::endl
                << "id            = " << id << std::endl;
      return 1;
    }
  }

  std::cout << "Chosen crystal = " << crystal[id].number << std::endl;


  //----------------------------------------------------------//
  //  MAIN LOOP                                               //
  //----------------------------------------------------------//

  std::cout << "-----" << std::endl;
  for(int iDet = 0 ; iDet < crystal[id].detectorSaturation.size(); iDet++)
  {
    std::cout << crystal[id].detectorSaturation[iDet].digitizerChannel << " " << crystal[id].detectorSaturation[iDet].saturation << " " << crystal[id].detectorSaturation[iDet].pedestal << " ";
  }
  std::cout << "-----" << std::endl;
  //MAIN LOOP
  long long int nevent = tree->GetEntries();
  std::cout << "Total number of events = " << nevent << std::endl;
  long int counter = 0;
  long int passCounter = 0;
  tree->SetNotify(formulasAnalysis);
  for (long long int i=0;i<nevent;i++)
  {
    // std::cout << "Event " << i << std::endl;
    tree->GetEvent(i);              //read complete accepted event in memory

    // std::cout << ChainExtendedTimeTag << "\t"
              // << ChainDeltaTimeTag    << "\t";
    // for (int j = 0 ; j < detector_channels.size() ; j++)
    // {
      // std::cout << charge[j] << "\t";
    // }
    // std::cout << std::endl;
    //
    counter++;

    if(crystal[id].FormulaTagAnalysis->EvalInstance()) // tag photopeak
    {
      if(crystal[id].FormulaAnalysis->EvalInstance())  //cut of crystal
      {
        passCounter++;
        out_doiFromTag      = crystal[id].taggingPosition;

        out_ExtendedTimeTag = ChainExtendedTimeTag;
        out_DeltaTimeTag    = ChainDeltaTimeTag;
        for(int iCh = 0; iCh < numOfCh ; iCh++)
        {
          if(charge[iCh] > 0) // old readout could give negative integrals (why?) anyway filter them out
          {
            if(iCh == crystal[id].taggingCrystalChannel)
            {
              out_charge[iCh] = charge[iCh];
            }
            for(int iDet = 0 ; iDet < crystal[id].detectorSaturation.size(); iDet++)
            {
              // std::cout << crystal[id].detectorSaturation[iDet].digitizerChannel << " " << crystal[id].detectorSaturation[iDet].saturation << " " << crystal[id].detectorSaturation[iDet].pedestal << " ";
              // sat correction is
              // - q_max * ln( 1 - (ch - ped)/(q_max))
              // where
              // q_max    = saturation
              // ch       = charge not corrected
              // ped      = pedestal (not corrected)

              if(iCh == crystal[id].detectorSaturation[iDet].digitizerChannel)
              {
                // std::cout << std::endl;
                // out_charge[iCh] = crystal[id].detectorSaturation[iDet].gain*(-crystal[id].detectorSaturation[iDet].saturation * TMath::Log(1.0 - ( (charge[iCh] - crystal[id].detectorSaturation[iDet].pedestal ) / (crystal[id].detectorSaturation[iDet].saturation) ) ));
                out_charge[iCh] = (-crystal[id].detectorSaturation[iDet].saturation * TMath::Log(1.0 - ( (charge[iCh] - crystal[id].detectorSaturation[iDet].pedestal ) / (crystal[id].detectorSaturation[iDet].saturation) ) ));

                // std::cout
                //           << crystal[id].detectorSaturation[iDet].saturation << " "
                //           << crystal[id].detectorSaturation[iDet].pedestal << " "
                //           << charge[iCh] << " " << out_charge[iCh]
                //           << std::endl;
              }
            }
          }
        }
        for(int iT = 0; iT < numOfT ; iT++)
        {
          out_timestamp[iT] = timeStamp[iT];
        }
        t1->Fill();
      }
    }



    //
    int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }

  std::cout << std::endl;
  // textfile.close();
  std::cout << counter << std::endl;
  std::cout << passCounter << std::endl;


  outputFile->cd();
  t1->Write();
  outputFile->Close();

  return 0;
}
