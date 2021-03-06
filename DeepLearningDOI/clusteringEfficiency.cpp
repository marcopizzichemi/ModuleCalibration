// compile with
// g++ -o ../../build/clusteringEfficiency clusteringEfficiency.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

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

// MPPC ID mapping
// D1      D2      D3      D4
// C1      C2      C3      C4
// B1      B2      B3      B4
// A1      A2      A3      A4
//
// Digitizer Channels mapping
// 18      17      22      21
// 19      16      23      20
// 24      27      28      31
// 25      26      29      30
//
// Crystal ID mapping
// 3       7       11      15
// 2       6       10      14
// 1       5       9       13
// 0       4       8       12




void usage()
{
  std::cout << "\t\t" << "[-i|--input] <file_prefix>  [-o|--output] <output.root> [-c|--calibration] <calibration.root> [-n|--crystal] <crystal number> [OPTIONS]" << std::endl
            << "\t\t" << "<file_prefix>                                     - prefix of TTree files to analyze"   << std::endl
            << "\t\t" << "<output.root>                                     - output file name"   << std::endl
            << "\t\t" << "<calibration.root>                                - output file name"   << std::endl
            << "\t\t" << "<crystal number>                                  - crystal to be analyzed"   << std::endl
            << "\t\t" << "--verbose                                         - "   << std::endl
            << "\t\t" << "--tagCut                                          - "   << std::endl
            << "\t\t" << "--triggerCut                                      - "   << std::endl
            << "\t\t" << "--broadCut                                        - "   << std::endl
            << "\t\t" << "--noiseCut                                        - "   << std::endl
            << "\t\t" << "--photopeakCut                                    - "   << std::endl
            << "\t\t" << "--cutgCut                                         - "   << std::endl
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
  std::string outputRootFileName = "";
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
      { "rootOutput", required_argument, 0, 0 },
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
    else if (c == 0 && optionIndex == 11){
      outputRootFileName = (char *)optarg;
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}


  // TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  std::ofstream outputFile;
  outputFile.open (outputFileName.c_str(),std::ofstream::out);
  TFile *outputRootFile = new TFile(outputRootFileName.c_str(),"RECREATE");

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
  // this time we add information on crystal selection
  ULong64_t out_DeltaTimeTag = 0;
  ULong64_t out_ExtendedTimeTag = 0;
  UShort_t  out_charge[64];
  Float_t   out_timestamp[64];
  UShort_t  out_tag = 0;
  UShort_t  out_clust = 0;
  UShort_t  out_ppeak = 0;
  UShort_t  out_max = 0;
  Short_t   out_crystalNumber = -1;
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
  // t1->Branch("zFromTag",&out_doiFromTag,"zFromTag/F");
  for (int i = 0 ; i < 64 ; i++)
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
  t1->Branch("tag",     &out_tag,   "tag/s");
  t1->Branch("clust",   &out_clust, "clust/s");
  t1->Branch("ppeak",    &out_ppeak,  "ppeak/s");
  t1->Branch("max",     &out_max,   "max/s");
  t1->Branch("crystal", &out_crystalNumber,   "crystal/S");


  //---------------------------------------------------------//


  //----------------------------------------------------------//
  //  Read and import from Calibration file                   //
  //----------------------------------------------------------//
  // open calibration file
  TFile *calibrationFile = new TFile(calibrationFileName.c_str());
  // prepare a list of TTreeFormula
  TList* formulasAnalysis = new TList();
  // prepare a vector of crystals (we will use only one)
  std::vector<Crystal_t> crystal;
  // read calibration file, fill the crystal structure with all relevant elements.
  // in particular, set the formulasAnalysis list, which contains the formulas used
  // in the loop to define the tag and crystal cuts. The tag cut is stored in
  // crystal[id].FormulaTagAnalysis
  // this includes only the tag photopeak condition. The crytal cut is in
  // crystal[id].FormulaAnalysis
  // and can be made of several cuts. The only relevant ones for real use here is Cutg
  // i.e. the 3D cut given by clustering. In general the user chooses what to include
  // in the formula (this routine is a general one used for other analysis), but here
  // the program should always be used with both --tagCut --cutgCut flags activated
  readCalibration(calibrationFile,          // calib file
                  tree,                     // input TChain (same for everyone)
                  formulasAnalysis,         // TList of all TTreeFormula
                  crystal,                  // structure of all crystals found in all calib lifes
                  UseTriggerChannelCut,     // include TriggerChannelCuthannel cut in crystalCut
                  UseBroadCut,              // include broadCut in crystalCut
                  UseCutNoise,              // include CutNoise in crystalCut
                  UsePhotopeakEnergyCut,    // include PhotopeakEnergyCut in crystalCut
                  UseCutgs                  // include CutGs in crystalCut
                 );

  // list the crystals with calibration data found
  std::cout << "Calibration data found for crystals: " << std::endl;
  for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  {
    if(crystal[i].accepted)
    {
      std::cout << crystal[i].number << std::endl;
    }
  }
  //----------------------------------------------------------//


  // switch off crystals if they were not chosen (only one chosen crystal allowed, sorry)
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
  //MAIN LOOP
  long long int nevent = tree->GetEntries();
  std::cout << "Total number of events = " << nevent << std::endl;
  // prepare counters
  long int counter = 0;
  long int tagCounter = 0;
  long int tagClustPPeakCounter = 0;
  long int tagPPeakCounter = 0;
  long int tagMaxCounter = 0;
  long int tagMaxPPeakCounter = 0;
  long int tagClustCounter = 0;
  // notify the formulas to the tree
  tree->SetNotify(formulasAnalysis);
  // loop on events
  for (long long int i=0;i<nevent;i++)
  {
    tree->GetEvent(i);              //read complete accepted event in memory
    counter++;                      // increase counter of events

    // copy ref time, charges and timestamps for all events. initialize flags
    out_DeltaTimeTag = ChainExtendedTimeTag;
    out_ExtendedTimeTag = ChainExtendedTimeTag;
    for (int iCh = 0 ; iCh < 64 ; iCh++)
    {
      out_charge[iCh]    = charge[iCh] ;
      out_timestamp[iCh] = timeStamp[iCh];
    }
    out_tag = 0;
    out_clust = 0;
    out_ppeak = 0;
    out_max = 0;
    out_crystalNumber = -1;
    // -----------

    // filter on TAG
    if(crystal[id].FormulaTagAnalysis->EvalInstance()) // filter on photopeak of tagging crystal
    {
      tagCounter++; // increase counter of events that pass the tagging photopeak cut (TAG)
      out_tag = 1; // write flag in output file

      // calculate the charge corrected by saturation, but just for the
      // "central" or trigger MPPC, since we are going to use single spectrum to cut photopeak
      UShort_t corrected_charge = 0;
      int iCh = crystal[id].detectorChannel;
      for(int iDet = 0 ; iDet < crystal[id].detectorSaturation.size(); iDet++)
      {
        // saturation correction is
        // - q_max * ln( 1 - (ch - ped)/(q_max))
        // where
        // q_max    = saturation
        // ch       = charge not corrected
        // ped      = pedestal (not corrected)
        // sat data is stored in the calibration file

        if( iCh == crystal[id].detectorSaturation[iDet].digitizerChannel)
        {
          corrected_charge = -crystal[id].detectorSaturation[iDet].saturation * TMath::Log(1.0 - ( (charge[iCh] - crystal[id].detectorSaturation[iDet].pedestal ) / (crystal[id].detectorSaturation[iDet].saturation) ) );
        }
      }

      // prapare flags. possibly not a great logical implementation, but it comes from stratification
      // and there's no time to change it now.
      bool isTriggerPhotopeak = false;
      bool isMax = true;

      // photopeak of trigger (aka central) mppc
      // luckily all single photopeaks are the same (no big surprise, since light sharing is ignored)
      // so the photopeak cut is beautifully hardcoded (shame!)
      if((corrected_charge > 30000) && (corrected_charge < 45000)  )
      {
        isTriggerPhotopeak = true; // set that trigger charge is in the photpeak of trigger spectrum
        out_ppeak = 1; // since this has already passed the TAG cut, strictly speaking the definition here is flawed, but it's not very important for now. they should be defined independently of TAG, but for sure the events not in TAG are ignored, so...
      }

      // charge in trigger is max
      // the loop is made on detectorSaturation structure just because it holds all crystals.
      // strictly speaking this is not really elegant, but anyway...
      for(int iDet = 0 ; iDet < crystal[id].detectorSaturation.size(); iDet++)
      {
        if( iCh == crystal[id].detectorSaturation[iDet].digitizerChannel)
        {
          // skip, it would compare charge with itself
        }
        else
        {
          // compare charges before saturation correction, it would make small difference
          // to correct before comparison (just detail)
          if(charge[iCh] < charge[crystal[id].detectorSaturation[iDet].digitizerChannel])
          {
            isMax = false; // set that trigger charge is not max charge in the event
          }
        }
      }

      if(isMax) // since this has already passed the TAG cut, strictly speaking the definition here is flawed, but it's not very important for now. they should be defined independently of TAG, but for sure the events not in TAG are ignored, so...
      {
        out_max = 1;
      }

      if(isTriggerPhotopeak)
      {
        tagPPeakCounter++; // increase counter of events in cut TAG + PPEAK
      }

      if(isMax)
      {
        tagMaxCounter++; // increase counter of events in cut TAG + MAX
      }

      if(isTriggerPhotopeak)
      {
        if(isMax)
        {
          tagMaxPPeakCounter++; // increase counter of events in cut TAG + MAX + PPEAK
        }
      }

      // apply full crystal cut (on the basis of what has been chosen by user with cmd line options)
      if(crystal[id].FormulaAnalysis->EvalInstance())  //cut of crystal
      {
        tagClustCounter++; // increase counter of events in cut TAG + CLUST
        out_clust = 1;
        out_crystalNumber = crystal[id].number;
        if(isTriggerPhotopeak)
        {
          tagClustPPeakCounter++; // increase counter of events in cut TAG + CLUST + PPEAK
        }
      }
    }

    // write the event in the output tree

    t1->Fill();

    int perc = ((100*counter)/nevent);
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }

  std::cout  << std::endl;
  // textfile.close();
  std::cout  << "# Legend " << std::endl;
  std::cout  << "# TAG   = events in photopeak of tagging crystal " << std::endl;
  std::cout  << "# CLUST = events in 3D cut, from clustering algorithm, for crystal aligned with tagging " << std::endl;
  std::cout  << "# PPEAK = events in photopeak SiPM0, i.e. the SiPM coupled to the crystal of the 4x4 array that is aligned with tagging" << std::endl;
  std::cout  << "# MAX   = events for which the charge of SiPM0 is greater than all the other charges recorded" << std::endl;
  std::cout  << std::endl;
  std::cout  << "All events                   = " << counter     << std::endl;
  std::cout  << "TAG                          = " << tagCounter  << std::endl;
  std::cout  << "TAG + CLUST                  = " << tagClustCounter << std::endl;
  std::cout  << "TAG + CLUST + PPEAK          = " << tagClustPPeakCounter << std::endl;
  std::cout  << "TAG + MAX                    = " << tagMaxCounter << std::endl;
  std::cout  << "TAG + MAX + PPEAK            = " << tagMaxPPeakCounter << std::endl;
  std::cout  << "TAG + PPEAK                  = " << tagPPeakCounter << std::endl;
  std::cout  << counter     << " "
             << tagCounter  << " "
             << tagClustCounter << " "
             << tagClustPPeakCounter << " "
             << tagMaxCounter << " "
             << tagMaxPPeakCounter << " "
             << tagPPeakCounter << std::endl;
  std::cout  << std::endl;

  // outputFile->cd();
  // t1->Write();
  // outputFile->Close();

  outputFile << "# Legend " << std::endl;
  outputFile << "# TAG   = events in photopeak of tagging crystal " << std::endl;
  outputFile << "# CLUST = events in 3D cut, from clustering algorithm, for crystal aligned with tagging " << std::endl;
  outputFile << "# PPEAK = events in photopeak SiPM0, i.e. the SiPM coupled to the crystal of the 4x4 array that is aligned with tagging" << std::endl;
  outputFile << "# MAX   = events for with the charge of SiPM0 is greater than all the other charges recorded" << std::endl;
  outputFile << std::endl;
  outputFile << "All events                   = " << counter     << std::endl;
  outputFile << "TAG                          = " << tagCounter  << std::endl;
  outputFile << "TAG + CLUST                  = " << tagClustCounter << std::endl;
  outputFile << "TAG + CLUST + PPEAK          = " << tagClustPPeakCounter << std::endl;
  outputFile << "TAG + MAX                    = " << tagMaxCounter << std::endl;
  outputFile << "TAG + MAX + PPEAK            = " << tagMaxPPeakCounter << std::endl;
  outputFile << "TAG + PPEAK                  = " << tagPPeakCounter << std::endl;
  outputFile << counter     << " "
             << tagCounter  << " "
             << tagClustCounter << " "
             << tagClustPPeakCounter << " "
             << tagMaxCounter << " "
             << tagMaxPPeakCounter << " "
             << tagPPeakCounter << std::endl;
  outputFile.close();


  outputRootFile->cd();
  t1->Write();
  outputRootFile->Close();

  return 0;
}
