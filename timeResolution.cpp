// template program to import n calibration data and loop on m acquisition ttrees
// use it as a base for analysis programs

// compile with
// g++ -o ../build/timeResolution timeResolution.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas


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

  // default args
  std::string calibrationFileNames = "";
  std::string inputFolderName = "./";
  std::string inputFilePrefix = "";
  std::string outputFileName = "outputBareboneFile.root";
  bool fitCorrection = false;
  Float_t doiFraction = 0.5;
  Float_t histoMin = -15e-9;//s
  Float_t histoMax = 15e-9;//s
  int histoBins = 500;
  int func = 0;
  Float_t fitPercMin = 5;
  Float_t fitPercMax = 6;
  int divs       = 10000;
  Float_t tagFwhm = 88.0e-12; //s //was 70.0e-12, then measured in 88.0e-12
  int excludeCh = -1;

  // parse command line arguments
  static struct option longOptions[] =
  {
      { "calibration", required_argument, 0, 0 },
      { "folder", required_argument, 0, 0 },
      { "prefix", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "fitCorrection", required_argument, 0, 0 },
      { "doiFraction", required_argument, 0, 0 },
      { "histoMin", required_argument, 0, 0 },
      { "histoMax", required_argument, 0, 0 },
      { "histoBins", required_argument, 0, 0 },
      { "func", required_argument, 0, 0 },
      { "fitPercMin", required_argument, 0, 0 },
      { "fitPercMax", required_argument, 0, 0 },
      { "divs", required_argument, 0, 0 },
      { "tagFwhm", required_argument, 0, 0 },
      { "exclude", required_argument, 0, 0 },

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
    else if (c == 0 && optionIndex == 4){
      fitCorrection = true;
    }
    else if (c == 0 && optionIndex == 5){
      doiFraction = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6){
      histoMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      histoMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      histoBins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      func = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      fitPercMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      fitPercMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 12){
      divs = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 13){
      tagFwhm = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 14){
      excludeCh = atoi((char *)optarg);
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

  //prepare output text file
  std::string textFileName = outputFileName.substr(0,outputFileName.size()-5);
  textFileName += ".txt";
  // std::cout << textFileName << std::endl;
  std::ofstream textfile;
  textfile.open (textFileName.c_str(),std::ofstream::out);



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
                    crystal,                  // structure of all crystals found in all calib files
                    false,                    // include TriggerChannelCuthannel cut in crystalCut
                    false,                    // include broadCut in crystalCut
                    false,                    // include CutNoise in crystalCut
                    true,                     // include PhotopeakEnergyCut in crystalCut
                    true                      // include CutGs in crystalCut
                    );



  }
  // optionally set w and z limits, and write values into crystal struct
  setWandZcuts(crystal);




  //prepare CTR histograms
  for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  {
    if(crystal[iCry].accepted)
    {
      std::stringstream sname;
      sname.str("");
      sname << "Central correction - Crystal " << crystal[iCry].number;
      crystal[iCry].centralCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
      sname.str("");
      sname << "Full correction - Crystal " << crystal[iCry].number;
      crystal[iCry].allCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
      sname.str("");
      sname << "No correction - Crystal " << crystal[iCry].number;
      crystal[iCry].simpleCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
      sname.str("");

      sname << "Polished correction - Crystal " << crystal[iCry].number;
      crystal[iCry].poliCorrCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
      sname.str("");

      sname << "CTR vs. W - Crystal " << crystal[iCry].number;
      crystal[iCry].ctrVSw = new TH2F(sname.str().c_str(),sname.str().c_str(),100,0,1,histoBins,histoMin,histoMax);
      sname.str("");

      // plots vs z

      sname << "Single ADC ch vs. Z - Crystal " << crystal[iCry].number;
      crystal[iCry].singleADCvsZ = new TH2F(sname.str().c_str(),sname.str().c_str(),100,0,crystal[iCry].length,200,0,200000);
      sname.str("");

      sname << "Total ADC ch vs. Z - Crystal " << crystal[iCry].number;
      crystal[iCry].totADCvsZ = new TH2F(sname.str().c_str(),sname.str().c_str(),100,0,crystal[iCry].length,200,0,200000);
      sname.str("");

      sname << "Basic CTR vs. Z - Crystal " << crystal[iCry].number;
      crystal[iCry].basicCTRvsZ = new TH2F(sname.str().c_str(),sname.str().c_str(),100,0,crystal[iCry].length,histoBins,histoMin,histoMax);
      sname.str("");

      sname << "Full CTR vs. Z - Crystal " << crystal[iCry].number;
      crystal[iCry].fullCTRvsZ = new TH2F(sname.str().c_str(),sname.str().c_str(),100,0,crystal[iCry].length,histoBins,histoMin,histoMax);
      sname.str("");
    }
  }


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

    bool keepEvent = true;
    if(keepEvent)
    {
      for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
      {
        if(crystal[iCry].accepted)
        {
          if(crystal[iCry].FormulaTagAnalysis->EvalInstance()) // if in photopeak of tagging crystal - or if in simulation
          {
            if(crystal[iCry].FormulaAnalysis->EvalInstance())  //if in global cut of crystal
            {
              goodEventsAnalysis++;

              //calculate FloodZ...
              float FloodZ = calculateFloodZ(charge,crystal[iCry]);
              // calculate reconstructed Z
              float z_reco = crystal[iCry].calibrationGraph->Eval(FloodZ);

              // calculate charge in trigger mppc (corrected by saturation)
              float centralCharge = calculate_trigger_charge(charge,crystal[iCry]);

              // calculate charge in the 9 relevant mppcs (corrected by saturation)
              float sumCharge = calculate_sum_charge(charge,crystal[iCry]);

              // Fill adc vs z scatter plots
              crystal[iCry].singleADCvsZ->Fill(z_reco,centralCharge);
              crystal[iCry].totADCvsZ->Fill(z_reco,sumCharge);


              //temp commented
              Float_t centralcorrection = 0.0;
              Float_t zeroCorrection    = 0.0;
              //no corr
              double simpleCTR = timeStamp[crystal[iCry].timingChannel] - timeStamp[crystal[iCry].taggingCrystalTimingChannel];

              if((timeStamp[crystal[iCry].timingChannel] != 0) && (timeStamp[crystal[iCry].taggingCrystalTimingChannel] != 0)) // no zeroes
              {
                crystal[iCry].simpleCTR->Fill(simpleCTR);
                crystal[iCry].basicCTRvsZ->Fill(z_reco,simpleCTR);
                // crystal[iCry].vSimple.push_back(simpleCTR);
              }

              if(crystal[iCry].tw_correction)
              {


                //skip event if is cut by min and maxAcceptedW
                // if(FloodZ > crystal[iCry].minAcceptedW && FloodZ < crystal[iCry].maxAcceptedW)
                if(true) // to avoid the w limit, it is actually useless
                {
                  double centralCTR;

                  if(fitCorrection)
                  {
                    centralcorrection = crystal[iCry].tw_correction_line->Eval(crystal[iCry].wz->Eval(crystal[iCry].length*doiFraction)) - crystal[iCry].tw_correction_line->Eval(FloodZ);
                  }
                  else
                  {
                    centralcorrection = crystal[iCry].tw_correction->Eval(crystal[iCry].wz->Eval(crystal[iCry].length*doiFraction)) - crystal[iCry].tw_correction->Eval(FloodZ);
                  }

                  centralCTR = (timeStamp[crystal[iCry].timingChannel] + (centralcorrection)) - timeStamp[crystal[iCry].taggingCrystalTimingChannel];

                  if((timeStamp[crystal[iCry].timingChannel] != 0) && (timeStamp[crystal[iCry].taggingCrystalTimingChannel] != 0)) // no zeroes
                  {
                    crystal[iCry].centralCTR->Fill(centralCTR);
                  }


                  // crystal[iCry].vCentral.push_back(centralCTR);

                  // modify to skip excluded channel
                  if(crystal[iCry].delay.size())
                  {
                    // begin of new way
                    float averageTimeStamp = 0.0;
                    float totalWeight = 0.0;

                    //first quickly check if there are zeroes
                    bool noZeroes = true;
                    if(timeStamp[crystal[iCry].taggingCrystalTimingChannel] == 0)
                    {
                      noZeroes = false;
                    }
                    for(unsigned int iDet = 0; iDet < crystal[iCry].correction_graphs.size(); iDet++)
                    {
                      int timingChannel = crystal[iCry].correction_graphs[iDet].timingChannel;
                      if(timingChannel == excludeCh)
                      {
                        //ignore
                      }
                      else
                      {
                        if(timeStamp[timingChannel] == 0)
                        {
                          noZeroes = false;
                        }
                      }

                    }
                    for(unsigned int iDet = 0; iDet < crystal[iCry].correction_graphs.size(); iDet++)
                    {
                      int timingChannel = crystal[iCry].correction_graphs[iDet].timingChannel;

                      if(timingChannel == excludeCh)
                      {
                        //ignore
                      }
                      else
                      {
                        float delay;

                        if(crystal[iCry].correction_graphs[iDet].isMainChannel)
                        {
                          delay = 0;
                        }
                        else
                        {
                          delay = crystal[iCry].correction_graphs[iDet].delay->Eval(FloodZ);
                        }
                        float correctedElement = timeStamp[timingChannel] - timeStamp[crystal[iCry].taggingCrystalTimingChannel] - delay;
                        if(correctedElement <=  histoMin || correctedElement >= histoMax )
                        {

                          noZeroes = false;
                        }
                      }

                    }

                    if(noZeroes)
                    {
                      for(unsigned int iDet = 0; iDet < crystal[iCry].correction_graphs.size(); iDet++)
                      {
                        //run on all the detectors, included the main one, but remember not to correct the main one for delay!
                        int timingChannel = crystal[iCry].correction_graphs[iDet].timingChannel;
                        if(timingChannel == excludeCh)
                        {
                           //ignore
                        }
                        else
                        {
                          float delay;
                          if(crystal[iCry].correction_graphs[iDet].isMainChannel)
                          {
                            delay = 0;
                          }
                          else
                          {
                            delay = crystal[iCry].correction_graphs[iDet].delay->Eval(FloodZ);
                          }
                          float delta = timeStamp[timingChannel] - timeStamp[crystal[iCry].taggingCrystalTimingChannel] - delay;
                          float weight = (1.0)/( TMath::Power(crystal[iCry].correction_graphs[iDet].rms->Eval(FloodZ),2) );

                          totalWeight += weight;
                          averageTimeStamp += delta*weight;
                        }


                      }
                      averageTimeStamp = averageTimeStamp/totalWeight;
                      crystal[iCry].ctrVSw->Fill(FloodZ,averageTimeStamp);
                      double allCTR = averageTimeStamp + centralcorrection;

                      crystal[iCry].allCTR->Fill(allCTR);
                      crystal[iCry].fullCTRvsZ->Fill(z_reco,allCTR);
                      // crystal[iCry].vAll.push_back(allCTR);

                    }
                  }
                }
              }

              if(crystal[iCry].polishedCorrection) //FIXME bugged now..
              {
                //central time stamp
                // std::cout << "----------------" << std::endl;

                float meanTimeStamp = 0.0;
                float sumWeight = 0.0;
                // std::cout << "crystal[iCry].fwhmForPolishedCorrection[0] = " << crystal[iCry].fwhmForPolishedCorrection[0]<< std::endl;

                bool noZeroes = true;
                if(timeStamp[crystal[iCry].taggingCrystalTimingChannel] == 0)
                {
                  noZeroes = false;
                }
                for(unsigned int iDet = 0; iDet < crystal[iCry].polished_correction.size(); iDet++)
                {
                  int timingChannel = crystal[iCry].polished_correction[iDet].timingChannel;
                  if(timingChannel == excludeCh)
                  {
                     //ignore
                  }
                  else
                  {
                    if(timeStamp[timingChannel] == 0)
                    {
                      noZeroes = false;
                    }
                  }

                }
                for(unsigned int iDet = 0; iDet < crystal[iCry].polished_correction.size(); iDet++)
                {
                  int timingChannel = crystal[iCry].polished_correction[iDet].timingChannel;
                  if(timingChannel == excludeCh)
                  {
                     //ignore
                  }
                  else
                  {
                    float delay = 0;
                    if(crystal[iCry].polished_correction[iDet].timingChannel == crystal[iCry].timingChannel)
                    {
                      delay = 0;
                    }
                    else
                    {
                      delay = crystal[iCry].polished_correction[iDet].mean;
                    }
                    float correctedElement = timeStamp[timingChannel] - timeStamp[crystal[iCry].taggingCrystalTimingChannel] - delay;
                    if(correctedElement <=  histoMin || correctedElement >= histoMax )
                    {
                      noZeroes = false;
                    }
                  }

                }

                if(noZeroes)
                {
                  // std::cout << "in " << crystal[iCry].polished_correction.size() <<  std::endl;
                  for(unsigned int iPoli = 0; iPoli < crystal[iCry].polished_correction.size(); iPoli++)
                  {
                    // std::cout << "ciao " << crystal[iCry].polished_correction.size() <<  std::endl;
                    // std::cout << iPoli <<" ";
                    if(crystal[iCry].polished_correction[iPoli].timingChannel == excludeCh)
                    {
                       //ignore
                    }
                    else
                    {

                      float delay = 0;
                      float weight = 0.0;
                      if(crystal[iCry].polished_correction[iPoli].timingChannel == crystal[iCry].timingChannel)
                      {
                        delay = 0;
                      }
                      else
                      {
                        delay = crystal[iCry].polished_correction[iPoli].mean;
                      }
                      float rms = crystal[iCry].polished_correction[iPoli].rms;
                      weight = pow(rms,-2);
                      sumWeight += weight;
                      float correctedTimepstamp = timeStamp[crystal[iCry].polished_correction[iPoli].timingChannel] - timeStamp[crystal[iCry].taggingCrystalTimingChannel] - delay;
                      meanTimeStamp += weight * correctedTimepstamp;
                    }

                  }
                  meanTimeStamp = meanTimeStamp/sumWeight;
                  double poliCorrCTR = meanTimeStamp;
                  crystal[iCry].poliCorrCTR->Fill(poliCorrCTR);
                  // crystal[iCry].vPoli.push_back(poliCorrCTR);
                }
              }
              // end of temp commented
            }
          }
        }
        // delete Formula;
      }
    }

    // LOOP ON ALL ACCEPTED CRYSTALS
    // I.E. ALL CRYSTALS IN THE CALIBRATION FILES
    // for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    // {
    //
    //   // do whatever you want
    //   // in particular, some tips:
    //   // 1) if(crystal[iCry].FormulaTagAnalysis->EvalInstance()){} --> this will select events in the photopeak of the external reference crystal, if there is one. Notice that each crystal can have its own FormulaTagAnalysis, since they can come from different calibration runs
    //   // 2) if(crystal[iCry].FormulaAnalysis->EvalInstance()) --> this will select events in a crystal (geometrical position in u-v-w and photopeak condition)
    //   // 3) to calculate w, if needed
    //   // float FloodZ = calculateFloodZ(charge,crystal[iCry]);
    //   // notice that the detectorSaturation values are written directly into the crystal struct, because they can be acquired in different conditions (but it has to be the same photodetector!)
    //
    //   //example here, selecting photopeak events in ext ref, photopeak events the crystal(s) found in calibration file(s), and calculating FloodZ for each
    //   if(crystal[iCry].accepted)
    //   {
    //     if(crystal[iCry].FormulaTagAnalysis->EvalInstance())
    //     {
    //       if(crystal[iCry].FormulaAnalysis->EvalInstance())
    //       {
    //         goodEventsAnalysis++;
    //         float FloodZ = calculateFloodZ(charge,crystal[iCry]);
    //         // std::cout << FloodZ << std::endl;
    //       }
    //     }
    //   }
    // }

    //LOOP COUNTER
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


  // do summary canvases for checking the fits
  int sqrtCrystals = ceil(sqrt( crystal.size() ) );
  TCanvas *cSumSimple  = new TCanvas("Summary Basic CTR","Summary Basic CTR",1200,1200);
  TCanvas *cSumCentral = new TCanvas("Summary Central CTR","Summary Central CTR",1200,1200);
  TCanvas *cSumAll     = new TCanvas("Summary Full CTR","Summary Full CTR",1200,1200);
  TCanvas *cPoliAll     = new TCanvas("Summary Polished CTR","Summary Polished CTR",1200,1200);
  cSumSimple ->Divide(sqrtCrystals,sqrtCrystals);
  cSumCentral->Divide(sqrtCrystals,sqrtCrystals);
  cSumAll->Divide(sqrtCrystals,sqrtCrystals);
  cPoliAll->Divide(sqrtCrystals,sqrtCrystals);


  // for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  // {
  //   TObjArray aSlices;
  //   crystal[iCry].ctrVSw->FitSlicesY(0,0,-1,0,"QNR",&aSlices);
  //
  //   crystal[iCry].aSlices = (TH1D*) aSlices[2];
  //
  //
  //   int nbinsx = crystal[iCry].aSlices->GetXaxis()->GetNbins();
  //   std::vector<double> x,y;
  //
  //   for (int i = 1; i < nbinsx ; i++)
  //   {
  //     if(crystal[iCry].aSlices->GetBinContent(i) > 0)
  //     {
  //       x.push_back(crystal[iCry].aSlices->GetBinCenter(i));
  //       y.push_back(1e12 * sqrt(2)*sqrt(pow(crystal[iCry].aSlices->GetBinContent(i),2)-pow(tagFwhm,2)))  ;
  //     }
  //   }
  //
  //    crystal[iCry].crtVSw_gr = new TGraph(x.size(),&x[0],&y[0]);
  //
  //
  //
  // }




  // crystal[iCry].ctrVSw->Write();


  std::cout << "Saving results..." << std::endl;
  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  // write whatever you want to save

  for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  {

    std::stringstream sname;



    Float_t realBasicCTRfwhm,realBasicCTRfwtm ;
    Float_t realCentralCTRfwhm,realCentralCTRfwtm;
    Float_t realAllCTRfwhm,realAllCTRfwtm;
    Float_t poliCorrCTRfwhm,poliCorrCTRfwtm;
    Float_t reallikeCTRfwhm,reallikeCTRfwtm;
    Float_t realhybridCTRfwhm,realhybridCTRfwtm;
    // double unbinnedSimpleCTR;
    // double unbinnedCentralCTR;
    // double unbinnedAllCTR;
    // double unbinnedPoliCTR;
    Double_t ret[4];
    Double_t fitRes[3];
    Float_t lightCentral;
    Float_t lightAll;

    // crystal[iCry].crtVSw_gr->Write();
    crystal[iCry].ctrVSw->Write();
    // crystal[iCry].aSlices->Write();
    // sigmas->Write();

    // get data on entries and light collected
    // std::cout << "quiiii" << std::endl;
    // crystal[iCry].basicCTRhisto->Write();

    // light central
    TF1 *gaussCentral = new TF1("gaussCentral","gaus");
    crystal[iCry].lightCentralHisto->Fit(gaussCentral,"Q");
    lightCentral = gaussCentral->GetParameter(1);
    crystal[iCry].lightCentralHisto->Write();
    // light central
    TF1 *gaussAll = new TF1("gaussAll","gaus");
    crystal[iCry].lightAllHisto->Fit(gaussAll,"Q");
    lightAll = gaussAll->GetParameter(1);
    crystal[iCry].lightAllHisto->Write();

    crystal[iCry].singleADCvsZ->Write();
    crystal[iCry].totADCvsZ->Write();
    crystal[iCry].basicCTRvsZ->Write();
    crystal[iCry].fullCTRvsZ->Write();




    if(crystal[iCry].simpleCTR)
    {
      Int_t CTRentries = crystal[iCry].simpleCTR->GetEntries();
      crystal[iCry].simpleCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].simpleCTR->SetFillStyle(3001);
      crystal[iCry].simpleCTR->SetFillColor(kGreen);
      crystal[iCry].simpleCTR->SetLineColor(kGreen);
      crystal[iCry].simpleCTR->SetStats(0);
      crystal[iCry].simpleCTR_norm = (TH1F*) crystal[iCry].simpleCTR->Clone();
      if(func == 0)
      {

        extractCTR(crystal[iCry].simpleCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractCTRWithEMG_withRef(crystal[iCry].simpleCTR,divs,tagFwhm,ret,fitRes);
          // extractWithGaussAndExp(crystal[iCry].simpleCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        // else
        // { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
        //   FindSmallestInterval(ret,crystal[iCry].simpleCTR,0.68,true,tagFwhm);
        // }
      }

      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      std::cout << "No corr    - cry " << crystal[iCry].number << "\t"
                << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << "\t"
                << CTRentries << "\t"
                << lightCentral << "\t"
                << lightAll << "\t"
                << fitRes[0] << "\t"
                << fitRes[1] << "\t"
                << fitRes[2] << "\t"
                << std::endl;


      //
      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      textfile  << "No corr    - cry " << crystal[iCry].number << "\t"
                << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << "\t"
                << CTRentries << "\t"
                << lightCentral << "\t"
                << lightAll << "\t"
                << fitRes[0] << "\t"
                << fitRes[1] << "\t"
                << fitRes[2] << "\t"
                << std::endl;

      realBasicCTRfwhm = ret[0]*1e12;
      realBasicCTRfwtm = ret[1]*1e12;
      // noCorr->Fill(ret[0]*1e12);

      crystal[iCry].simpleCTR->Write();
      cSumSimple->cd(iCry+1);
      crystal[iCry].simpleCTR->Draw();
      crystal[iCry].simpleCTR_norm->Scale(1.0/crystal[iCry].simpleCTR_norm->GetMaximum());

      // use unbinned method
      // if(unbinned)
      // {
      //   double mean,meanErr,min,max;
      //   double delta = FindSmallestInterval(mean,
      //                                       meanErr,
      //                                       min,
      //                                       max,
      //                                       crystal[iCry].vSimple,
      //                                       0.68,
      //                                       true);
      //   //now pass to fwhm
      //   double fwhm = 2.355 * (delta/2.0);
      //   unbinnedSimpleCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
      //   // unbinnednoCorr->Fill(unbinnedSimpleCTR);
      //
      //   std::cout << "Unbinned No corr    - cry " << crystal[iCry].number << "\t"
      //             << unbinnedSimpleCTR << "\t"
      //             << 0 << std::endl;
      //
      //   textfile  << "Unbinned No corr    - cry " << crystal[iCry].number << "\t"
      //             << unbinnedSimpleCTR << "\t"
      //             << 0 << std::endl;
      // }
    }

    if(crystal[iCry].centralCTR)
    {
      Int_t CTRentries = crystal[iCry].centralCTR->GetEntries();
      crystal[iCry].centralCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].centralCTR->SetFillStyle(3001);
      crystal[iCry].centralCTR->SetFillColor(kBlue);
      crystal[iCry].centralCTR->SetLineColor(kBlue);
      crystal[iCry].centralCTR->SetStats(0);
      crystal[iCry].centralCTR_norm = (TH1F*) crystal[iCry].centralCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractCTRWithEMG_withRef(crystal[iCry].centralCTR,divs,tagFwhm,ret,fitRes);
          // extractWithGaussAndExp(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        // else
        // { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
        //   FindSmallestInterval(ret,crystal[iCry].centralCTR,0.68,true,tagFwhm);
        // }
      }


      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      std::cout << "Central    - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Central    - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      realCentralCTRfwhm = ret[0]*1e12;
      realCentralCTRfwtm = ret[1]*1e12;
      // centralCorr->Fill(ret[0]*1e12);

      crystal[iCry].centralCTR->Write();
      cSumCentral->cd(iCry+1);
      crystal[iCry].centralCTR->Draw();
      // crystal[iCry].centralCTR->Scale(1.0/crystal[iCry].centralCTR->GetMaximum());
      crystal[iCry].centralCTR_norm->Scale(1.0/crystal[iCry].centralCTR_norm->GetMaximum());

      // use unbinned method
      // if(unbinned)
      // {
      //   double mean,meanErr,min,max;
      //   double delta = FindSmallestInterval(mean,
      //                                       meanErr,
      //                                       min,
      //                                       max,
      //                                       crystal[iCry].vCentral,
      //                                       0.68,
      //                                       true);
      //   //now pass to fwhm
      //   double fwhm = 2.355 * (delta/2.0);
      //   unbinnedCentralCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
      //   // unbinnedcentralCorr->Fill(unbinnedCentralCTR);
      //
      //   std::cout << "Unbinned Central    - cry " << crystal[iCry].number << "\t"
      //             << unbinnedCentralCTR << "\t"
      //             << 0 << std::endl;
      //
      //   textfile  << "Unbinned Central    - cry " << crystal[iCry].number << "\t"
      //             << unbinnedCentralCTR << "\t"
      //             << 0 << std::endl;
      // }
    }

    if(crystal[iCry].allCTR)
    {
      Int_t CTRentries = crystal[iCry].allCTR->GetEntries();
      crystal[iCry].allCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].allCTR->SetFillStyle(3001);
      crystal[iCry].allCTR->SetFillColor(kRed);
      crystal[iCry].allCTR->SetLineColor(kRed);
      crystal[iCry].allCTR->SetStats(0);

      crystal[iCry].allCTR_norm = (TH1F*) crystal[iCry].allCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractCTRWithEMG_withRef(crystal[iCry].allCTR,divs,tagFwhm,ret,fitRes);
          // extractWithGaussAndExp(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        // else
        // { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
        //   FindSmallestInterval(ret,crystal[iCry].allCTR,0.68,true,tagFwhm);
        // }
      }
      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      std::cout << "Full corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Full corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      realAllCTRfwhm = ret[0]*1e12;
      realAllCTRfwtm = ret[1]*1e12;
      // fullCorr->Fill(ret[0]*1e12);
      crystal[iCry].allCTR->Write();
      cSumAll->cd(iCry+1);
      crystal[iCry].allCTR->Draw();
      crystal[iCry].allCTR_norm->Scale(1.0/crystal[iCry].allCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());

      // use unbinned method
      // if(unbinned)
      // {
      //   double mean,meanErr,min,max;
      //   double delta = FindSmallestInterval(mean,
      //                                       meanErr,
      //                                       min,
      //                                       max,
      //                                       crystal[iCry].vAll,
      //                                       0.68,
      //                                       true);
      //   //now pass to fwhm
      //   double fwhm = 2.355 * (delta/2.0);
      //   unbinnedAllCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
      //   // unbinnedfullCorr->Fill(unbinnedAllCTR);
      //
      //   std::cout << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedAllCTR << "\t"
      //             << 0 << std::endl;
      //
      //   textfile  << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedAllCTR << "\t"
      //             << 0 << std::endl;
      // }
    }

    if(crystal[iCry].poliCorrCTR)
    {
      Int_t CTRentries = crystal[iCry].poliCorrCTR->GetEntries();
      crystal[iCry].poliCorrCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].poliCorrCTR->SetFillStyle(3001);
      crystal[iCry].poliCorrCTR->SetFillColor(kBlack);
      crystal[iCry].poliCorrCTR->SetLineColor(kBlack);
      crystal[iCry].poliCorrCTR->SetStats(0);

      crystal[iCry].poliCorrCTR_norm = (TH1F*) crystal[iCry].poliCorrCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractCTRWithEMG_withRef(crystal[iCry].poliCorrCTR,divs,tagFwhm,ret,fitRes);
          // extractWithGaussAndExp(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        // else
        // { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
        //   FindSmallestInterval(ret,crystal[iCry].poliCorrCTR,0.68,true,tagFwhm);
        // }
      }


      // std::cout << "Polished corr. - cry " << crystal[iCry].number << "\t"
      //           << ret[0]*1e12 << "\t"
      //           << ret[1]*1e12 << std::endl;
      // //
      // std::cout << "Polished FIT   - cry " << crystal[iCry].number << "\t"
      //           << fitRes[0] << "\t"
      //           << fitRes[1] << "\t"
      //           << fitRes[2] << "\t"
      //           << std::endl;
      //
      // textfile  << "Polished corr. - cry " << crystal[iCry].number << "\t"
      //           << ret[0]*1e12 << "\t"
      //           << ret[1]*1e12 << std::endl;
      // //
      // textfile << "Polished FIT   - cry " << crystal[iCry].number << "\t"
      //           << fitRes[0] << "\t"
      //           << fitRes[1] << "\t"
      //           << fitRes[2] << "\t"
      //           << std::endl;


      //
      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      std::cout << "Polished corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Polished corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      poliCorrCTRfwhm = ret[0]*1e12;
      poliCorrCTRfwtm = ret[1]*1e12;
      // poliCorr->Fill(ret[0]*1e12);
      crystal[iCry].poliCorrCTR->Write();
      cPoliAll->cd(iCry+1);
      crystal[iCry].poliCorrCTR->Draw();
      crystal[iCry].poliCorrCTR_norm->Scale(1.0/crystal[iCry].poliCorrCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());


      // use unbinned method
      // if(unbinned)
      // {
      //   double mean,meanErr,min,max;
      //   double delta = FindSmallestInterval(mean,
      //                                       meanErr,
      //                                       min,
      //                                       max,
      //                                       crystal[iCry].vPoli,
      //                                       0.68,
      //                                       true);
      //   //now pass to fwhm
      //   double fwhm = 2.355 * (delta/2.0);
      //   unbinnedPoliCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
      //   // unbinnedpoliCorr->Fill(unbinnedPoliCTR);
      //
      //   std::cout << "Unbinned Polished corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedPoliCTR << "\t"
      //             << 0 << std::endl;
      //
      //   textfile  << "Unbinned Polished corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedPoliCTR << "\t"
      //             << 0 << std::endl;
      // }
    }



    sname.str("");

    sname << "Summary - Crystal " << crystal[iCry].number;
    TCanvas* c_summary = new TCanvas(sname.str().c_str(),sname.str().c_str(),1200,800);
    c_summary->cd();
    THStack *hs = new THStack("hs","");
    hs->Add(crystal[iCry].simpleCTR_norm);
    hs->Add(crystal[iCry].centralCTR_norm);
    hs->Add(crystal[iCry].allCTR_norm);
    hs->Add(crystal[iCry].poliCorrCTR_norm);


    // std::cout << "Crystal " << crystal[iCry].number << std::endl;
    // std::cout << "CTR FWHM [ps] " << std::endl;
    hs->Draw("hist nostack");
    sname.str("");
    sname << "CTR - Crystal " << crystal[iCry].number << " - width in FWHM";
    hs->SetTitle(sname.str().c_str());
    hs->GetXaxis()->SetTitle("Time [s]");
    hs->GetXaxis()->SetTitleOffset(1);
    hs->GetXaxis()->SetTitleSize(0.045);
    hs->GetXaxis()->SetLabelSize(0.045);
    hs->GetYaxis()->SetLabelSize(0.045);
    TLegend *legend = new TLegend(0.54,0.62,0.89,0.89,"");
    legend->SetFillStyle(0);
    if(crystal[iCry].simpleCTR)
    {
      sname.str("");
      sname << "No correction        = " << realBasicCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].simpleCTR,sname.str().c_str(),"f");
      // std::cout << "No correction       = "<< realBasicCTRfwhm   << " ps" << std::endl;
    }
    if(crystal[iCry].centralCTR)
    {
      sname.str("");
      sname << "Central correction = " << realCentralCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].centralCTR,sname.str().c_str(),"f");
      // std::cout << "Central correction  = "<< realCentralCTRfwhm << " ps" << std::endl;
    }
    if(crystal[iCry].allCTR)
    {
      sname.str("");
      sname << "Full correction       = " << realAllCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].allCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }

    if(crystal[iCry].poliCorrCTR)
    {
      sname.str("");
      sname << "Polished correction       = " << poliCorrCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].poliCorrCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }



    sname.str("");
    legend->Draw();
    gStyle->SetOptTitle(0);
    TPaveLabel *title = new TPaveLabel(.11,.95,.35,.99,"new title","brndc");
    title->Draw();
    // std::cout << std::endl;

    c_summary->Write();


    TH1F* cloneBasic;
    TH1F* cloneCentral;
    TH1F* cloneAll;
    // TH1F* clonePoli;
    THStack *cloneHs = (THStack*) hs->Clone();
    TLegend *legend1 = new TLegend(0.15,0.69,0.49,0.89,"");
    legend1->SetFillStyle(0);
    sname.str("");
    sname << "Multi - Crystal " << crystal[iCry].number;
    TCanvas* c_multi = new TCanvas(sname.str().c_str(),sname.str().c_str(),1800,1400);
    c_multi->Divide(2,2);

    if(crystal[iCry].simpleCTR_norm)
    {
      cloneBasic   = (TH1F*) crystal[iCry].simpleCTR->Clone();
      c_multi->cd(1);
      sname.str("");
      sname << "No correction        = " << realBasicCTRfwhm << "ps";
      legend1->AddEntry(cloneBasic,sname.str().c_str(),"f");
      cloneBasic->Draw();
      legend1->Draw();
    }
    if(crystal[iCry].centralCTR_norm)
    {
      cloneCentral = (TH1F*) crystal[iCry].centralCTR->Clone();
      c_multi->cd(2);
      TLegend *legend2 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend2->SetFillStyle(0);
      sname.str("");
      sname << "Central correction   = " << realCentralCTRfwhm << "ps";
      legend2->AddEntry(cloneCentral,sname.str().c_str(),"f");
      cloneCentral->Draw();
      legend2->Draw();
    }
    if(crystal[iCry].allCTR_norm)
    {
      cloneAll     = (TH1F*) crystal[iCry].allCTR->Clone();
      c_multi->cd(3);
      TLegend *legend3 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend3->SetFillStyle(0);
      sname.str("");
      sname << "Full correction      = " << realAllCTRfwhm << "ps";
      legend3->AddEntry(cloneAll,sname.str().c_str(),"f");
      cloneAll->Draw();
      legend3->Draw();
    }

    c_multi->cd(4);
    c_multi->cd(4)->SetGrid();
    cloneHs->Draw("hist nostack");
    c_multi->Write();

  }

  // noCorr->Write();
  // centralCorr->Write();
  // fullCorr->Write();
  // poliCorr->Write();
  // unbinnednoCorr->Write();
  // unbinnedcentralCorr->Write();
  // unbinnedfullCorr->Write();
  // unbinnedpoliCorr->Write();
  cSumSimple ->Write();
  cSumCentral->Write();
  cSumAll->Write();
  cPoliAll->Write();

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
            << "[-o|--output] <output>\t [OPTIONS]"
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
            << "\t\t" << "--fitCorrection                                    - use line fit to perform correction   - default = not given (false)"  << std::endl
            << "\t\t" << "--doiFraction <value>                              - fraction of DOI length towards which the time stamps are corrected (from 0 to 1)"  << std::endl
            << "\t\t" << "--histoMin <value>                                 - lower limit of CTR spectra, in sec - default = -15e-9"  << std::endl
            << "\t\t" << "--histoMax <value>                                 - upper limit of CTR spectra, in sec - default = 15e-9"  << std::endl
            << "\t\t" << "--histoBins <value>                                - n of bins for CTR spectra - default = 500"  << std::endl
            << "\t\t" << "--func <value>                                     - function for fitting (default = 0)"  << std::endl
            << "\t\t" << "--fitPercMin <value>                               - time fit min is set to ((gauss fit mean) - fitPercMin*(gauss fit sigma))  - default = 5"  << std::endl
            << "\t\t" << "--fitPercMax <value>                               - time fit max is set to ((gauus fit mean) - fitPercMax*(gauss fit sigma))  - default = 6" << std::endl
            << "\t\t" << "--divs <value>                                     - n of divisions when looking for FWHM - default = 10000"  << std::endl
            << "\t\t" << "--tagFwhm <value>                                  - FWHM timing resolution of reference board, in sec - default = 88e-12"  << std::endl
            << "\t\t" << "--exclude <value>                                  - number of timing channel to exclude from corrections"  << std::endl


            << std::endl;
}
