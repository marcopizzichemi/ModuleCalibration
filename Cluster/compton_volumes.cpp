// template program to import n calibration data and loop on m acquisition ttrees
// use it as a base for analysis programs

// compile with
// g++ -o ../build/compton_volumes compton_volumes.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas


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

#include <iomanip>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "../libraries/CrystalStructs.h"     // Crystal_t , detector_t
#include "../libraries/Calibration.h"        // readTaggingData , readCalibration , setWandZcuts
#include "../libraries/Utilities.h"          // read_directory , invert_a_matrix
#include "../libraries/Extract.h"            // extractWithEMG , extractCTR , FindSmallestInterval

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
  std::string cutgFileName = "";
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
  // bool cluster = false;
  int Ncrystals = 0;

  // parse command line arguments
  static struct option longOptions[] =
  {
      { "calibration", required_argument, 0, 0 },
      { "folder", required_argument, 0, 0 },
      { "prefix", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "cutgfile", required_argument, 0, 0 },
      // { "crystals", required_argument, 0, 0 },
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
      cutgFileName = (char *)optarg;
    }
    // else if (c == 0 && optionIndex == 4){
    //   Ncrystals = atof((char *)optarg);
    // }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}
  
  // if(Ncrystals == 0)
  // {
  //   std::cout << std::endl;
  //   std::cout << "ERROR! You need to provide the number of crystals!" << std::endl;
  //   std::cout << "See program usage below..." << std::endl;
  //   std::cout << std::endl;
  //   std::cout << argv[0];
  //   usage();
  //   return 1;
  // }
  

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
  if(cutgFileName == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide the cutg file!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }

  // //prepare output text file
  // std::string textFileName = outputFileName.substr(0,outputFileName.size()-5);
  // textFileName += ".txt";
  // // std::cout << textFileName << std::endl;
  // std::ofstream textfile;
  // textfile.open (textFileName.c_str(),std::ofstream::out);



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
  // GET CUTG FILE                    //
  //----------------------------------//
  // it works like this.
  // 1. run standard modulecalibration without manualCutG option
  // 2. take a look at the 2d projections of Flood Histogram 3D for the mppc under analysis (the _zx and _zy) and the standard flood 2D
  // 3. create a graphical cut for one of them, like this:
  //    a) open toolbar in gui
  //    b) click on scissors symbol (Graphical Cut)
  //    c) draw a cut by left clicking on points, then close it with double click (cut has to be a closed line)
  //    d) retrieve from ROOT "world" and change name
  // 4. save the 3 cutgs to a root file
  // In details:
  //    - go to mppc folder
  //     - start with [...]_zx 
  //        // do graphical selection
  //        TCutG *mycutg1;
  //        mycutg1 = (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG");
  //        mycutg1->SetName("cutg_1");
  //        // open one of the crystal folder
  //        // adjust following cutg_xz_1_B2 changing 1 to other number 
  //        mycutg1->SetVarX(cutg_xz_1_B2->GetVarX())
  //        mycutg1->SetVarY(cutg_xz_1_B2->GetVarY())
  //     - then with [...]_zy 
  //        // do graphical selection
  //        TCutG *mycutg2;
  //        mycutg2 = (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG");
  //        mycutg2->SetName("cutg_2");
  //        // open again same crystal folder
  //        // adjust following cutg_xz_1_B2 changing 1 to other number 
  //        mycutg2->SetVarX(cutg_yz_1_B2->GetVarX())
  //        mycutg2->SetVarY(cutg_yz_1_B2->GetVarY())
  //     - finally with the flood histo 2D of this mppc
  //        // do graphical selection
  //        TCutG *mycutg3;
  //        mycutg3 = (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG");
  //        mycutg3->SetName("cutg_3");
  //        // open again same crystal folder
  //        // adjust following cutg_xz_1_B2 changing 1 to other number 
  //        mycutg3->SetVarX(cutg_xy_1_B2->GetVarX())
  //        mycutg3->SetVarY(cutg_xy_1_B2->GetVarY())
  //    g) save all in a tfile (cutgs.root)
  //        TFile *fOut = new TFile("cutgs.root","RECREATE")
  //        fOut->cd()
  //        mycutg1->Write()
  //        mycutg2->Write()
  //        mycutg3->Write()
  //        fOut->Close()
  TFile *fileCutg = TFile::Open(cutgFileName.c_str());
  TCutG **cutg_ext = new TCutG *[3];
  cutg_ext[0] = (TCutG*) gDirectory->Get("cutg_1");
  cutg_ext[1] = (TCutG*) gDirectory->Get("cutg_2");
  cutg_ext[2] = (TCutG*) gDirectory->Get("cutg_3");



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
  std::vector<int> crystal_numbers;
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
  int numOfcrystals = 0;
  std::string ch_prefix("ch");
  std::string t_prefix("t");
  std::string en_prefix("enCry");
  for(int i = 0 ; i < nLeaves ; i++)
  {
    if (!leavesName[i].compare(0, ch_prefix.size(), ch_prefix))
    {
      numOfCh++;
      detector_channels.push_back(atoi( (leavesName[i].substr(ch_prefix.size(),leavesName[i].size()-ch_prefix.size())).c_str() )) ;
    }
    if (!leavesName[i].compare(0, en_prefix.size(), en_prefix))
    {
      numOfcrystals++;
      crystal_numbers.push_back(atoi( (leavesName[i].substr(en_prefix.size(),leavesName[i].size()-en_prefix.size())).c_str() )) ;
    }
  }
  std::cout << "Detector Channels  \t= " << numOfCh << std::endl;
  std::cout << "Number of Crystals \t= " << numOfcrystals << std::endl;

  //set variables and branches
  ULong64_t     ChainExtendedTimeTag;                                // extended time tag
  ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous
  // Float_t      *charge;
  UShort_t      *charge;
  // Long_t      *timeStamp;
  Float_t      *timeStamp;
  Float_t      *energyPerCrystal;
  Short_t      CrystalsHit;
  Short_t      CrystalNumber;
  Float_t      TotalEnergyDeposited;
  Float_t      TotalEnergyDepositedInCrystals;
  TBranch      *bChainExtendedTimeTag;                               // branches for above data
  TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  TBranch      *bCrystalsHit;
  TBranch      *bCrystalNumber;
  TBranch      *bTotalEnergyDeposited;
  TBranch      *bTotalEnergyDepositedInCrystals;
  TBranch      **bCharge;
  TBranch      **btimeStamp;
  TBranch      **benergyPerCrystal;
  // charge = new Float_t[numOfCh];
  charge = new UShort_t[numOfCh];
  // timeStamp = new Long_t[numOfCh];
  timeStamp = new Float_t[numOfCh];
  energyPerCrystal = new Float_t[numOfcrystals];
  bCharge = new TBranch*[numOfCh];
  btimeStamp = new TBranch*[numOfCh];
  benergyPerCrystal = new TBranch*[numOfcrystals];
  // set branches for reading the input files
  tree->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag, &bChainExtendedTimeTag);
  tree->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag, &bChainDeltaTimeTag);
  for (int i = 0 ; i < detector_channels.size() ; i++)
  {
    std::stringstream sname;
    sname << "ch" << detector_channels[i];
    tree->SetBranchAddress(sname.str().c_str(),&charge[detector_channels[i]],&bCharge[detector_channels[i]]);
    sname.str("");
    sname << "t" << detector_channels[i];
    tree->SetBranchAddress(sname.str().c_str(),&timeStamp[detector_channels[i]],&btimeStamp[detector_channels[i]]);
    sname.str("");
  }
  for (int i = 0 ; i < crystal_numbers.size() ; i++)
  {
    std::stringstream sname;
    sname << "enCry" << crystal_numbers[i];
    tree->SetBranchAddress(sname.str().c_str(),&energyPerCrystal[detector_channels[i]],&benergyPerCrystal[detector_channels[i]]);
    sname.str("");
  }
  
  tree->SetBranchAddress("CrystalsHit", &CrystalsHit, &bCrystalsHit);
  tree->SetBranchAddress("CrystalNumber", &CrystalNumber, &bCrystalNumber);
  tree->SetBranchAddress("TotalEnergyDeposited", &TotalEnergyDeposited, &bTotalEnergyDeposited);
  tree->SetBranchAddress("TotalEnergyDepositedInCrystals", &TotalEnergyDepositedInCrystals, &bTotalEnergyDepositedInCrystals);
  
  
  

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
  // std::cout << "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" << std::endl;
  // optionally set w and z limits, and write values into crystal struct
  setWandZcuts(crystal);
  
  // list the crystals with calibration data found
  std::cout << "Calibration data found for crystals: " << std::endl;
  for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  {
    if(crystal[i].accepted)
    {
      crystal[i].polishedCorrection = false; // FIXME hardcode stop to polish correction
      std::cout << crystal[i].number << std::endl;
    }
  }
  
  //save the cutg anyway
  TCut ComptonCutG ;
  for(int gg = 0 ; gg < 3; gg++)
  {
    ComptonCutG += cutg_ext[gg]->GetName();
  }
  TTreeFormula* FormulaComptonCutG = new TTreeFormula("FormulaComptonCutG",ComptonCutG,tree);
  formulasAnalysis->Add(FormulaComptonCutG);

  // MAIN LOOP
  long long int counter = 0;
  tree->SetNotify(formulasAnalysis);
  long long int neventAnalysis = tree->GetEntries();
  std::cout << "Total number of events in analysis file = " << neventAnalysis << std::endl;
  long int goodEventsAnalysis = 0;
  long int counterAnalysis = 0;
  // cluster analysis, optional
  long int inCrystal = 0;
  long int inCrystalGood = 0;
  long int inCrystalBad = 0;
  long int outOfCrystalButGood = 0;
  

  goodEventsAnalysis = 0;
  counterAnalysis = 0;
  counter = 0;
  

  std::cout << "Saving results..." << std::endl;
  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");

  TTree *TreeCluster;

  TreeCluster = new TTree("cluster","cluster");
  int cryNum = -1;
  float totalEn = 0.;
  Short_t      out_CrystalsHit = 0;
  Short_t      out_CrystalNumber = 0;
  Float_t      out_TotalEnergyDeposited = 0;
  Float_t      out_TotalEnergyDepositedInCrystals = 0;
  // Float_t      out_u = 0.;
  // Float_t      out_v = 0.;
  // Float_t      out_w = 0.;
  Short_t cryNumAssignedByCluster = -1;
  Short_t inCry      = 0;
  Short_t inCompton  = 0;
  Short_t inCryGood  = 0;
  Short_t inCryBad   = 0;
  Short_t outCryGood = 0;
  // also save the ch and t data
  UShort_t *out_charge; //adc type
  out_charge = new UShort_t[numOfCh];
  Float_t *out_timestamp;
  out_timestamp = new Float_t[numOfCh];
  Float_t* out_energyPerCrystal;
  out_energyPerCrystal = new Float_t[numOfcrystals];

  TreeCluster->Branch("CrystalsHit",&out_CrystalsHit,"CrystalsHit/S");
  TreeCluster->Branch("CrystalNumber",&out_CrystalNumber,"CrystalNumber/S");
  TreeCluster->Branch("TotalEnergyDeposited",&out_TotalEnergyDeposited,"TotalEnergyDeposited/F");
  TreeCluster->Branch("TotalEnergyDepositedInCrystals",&out_TotalEnergyDepositedInCrystals,"TotalEnergyDepositedInCrystals/F");
  TreeCluster->Branch("cryNumAssignedByCluster",&cryNumAssignedByCluster,"cryNumAssignedByCluster/S");
  TreeCluster->Branch("inCrystal",&inCry,"inCrystal/S");
  TreeCluster->Branch("inCompton",&inCompton,"inCompton/S");
  TreeCluster->Branch("inCrystalGood",&inCryGood,"inCrystalGood/S");
  TreeCluster->Branch("inCrystalBad",&inCryBad,"inCrystalBad/S");
  TreeCluster->Branch("outOfCrystalButGood",&outCryGood,"outOfCrystalButGood/S");
  // TreeCluster->Branch("u",&out_u,"u/F");
  // TreeCluster->Branch("v",&out_v,"v/F");
  // TreeCluster->Branch("w",&out_w,"w/F");
  std::stringstream out_snames,out_stypes;
  for (int i = 0 ; i < numOfCh ; i++)
  {
    //empty the stringstreams
    out_snames.str("");
    out_stypes.str("");
    out_charge[i] = 0;

    out_snames << "ch" << i;
    out_stypes << "ch" << i << "/s";
    TreeCluster->Branch(out_snames.str().c_str(),&out_charge[i],out_stypes.str().c_str());
    out_snames.str("");
    out_stypes.str("");

  }
  for (int i = 0 ; i < numOfCh ; i++)
  {
    //empty the stringstreams
    out_snames.str("");
    out_stypes.str("");
    out_timestamp[i] = 0;
    out_snames << "t" << i;
    out_stypes << "t" << i << "/F";
    TreeCluster->Branch(out_snames.str().c_str(),&out_timestamp[i],out_stypes.str().c_str());
    out_snames.str("");
    out_stypes.str("");
  }
  for (int i = 0 ; i < numOfcrystals ; i++)
  {
    //empty the stringstreams
    out_snames.str("");
    out_stypes.str("");
    out_energyPerCrystal[i] = 0.;
    out_snames << "enCry" << i;
    out_stypes << "enCry" << i << "/F";
    TreeCluster->Branch(out_snames.str().c_str(),&out_energyPerCrystal[i],out_stypes.str().c_str());
    out_snames.str("");
    out_stypes.str("");
  }

  // MAIN LOOP
  
  for (long long int i=0;i<neventAnalysis;i++)
  {
    inCry      = 0;
    inCryGood  = 0;
    inCryBad   = 0;
    outCryGood = 0;
    tree->GetEvent(i);              //read complete accepted event in memory
    out_CrystalsHit = CrystalsHit;
    out_CrystalNumber = CrystalNumber;
    out_TotalEnergyDeposited = TotalEnergyDeposited;
    out_TotalEnergyDepositedInCrystals = TotalEnergyDepositedInCrystals;

    // assign out charges and timestamps
    for (int c = 0 ; c < numOfCh ; c++)
    {
      out_charge[c] = charge[c];
      out_timestamp[c] = timeStamp[c];
    }
    for (int c = 0 ; c < numOfcrystals ; c++)
    {
      out_energyPerCrystal[c] = energyPerCrystal[c];
    }
    // calculate (and save) always uvw
    // out_u      = 0.;
    // out_v      = 0.;
    // out_w      = calculateFloodZ_withoutCorrectingForSaturation(charge,crystal[iCry]);
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      // std::cout << "probing " << crystal[iCry].number << std::endl;
      if(crystal[iCry].accepted)
      {
        // std::cout << "accepted " << crystal[iCry].number << std::endl;
        if(crystal[iCry].FormulaJustCutG->EvalInstance())  //all crystal event
        {
          

          cryNumAssignedByCluster = crystal[iCry].number;
          // std::cout << "inCrystal " << crystal[iCry].number << std::endl;
          // cluster marked event in crystal
          inCrystal++;
          inCry = 1;
          if(cryNumAssignedByCluster == CrystalNumber)
          {
            inCrystalGood++;
            inCryGood = 1;
            // inCryBad = 0;
          }
          else
          {
            inCrystalBad++;
            // inCryGood = 0;
            inCryBad = 1;
          }
          // check num
        }
        else 
        {
          
          // inCry = 0;
          if(crystal[iCry].number == CrystalNumber)
          {
            outOfCrystalButGood++;
            outCryGood = 1;
          }
          else 
          {
            // outCryGood = 0;
          }

        }
      }
    }
    if(inCry == 0) // so if the event was not falling into the geometrical cut of any crystal
    {
      if(FormulaComptonCutG->EvalInstance())  //all compton cut
      {
        inCompton = 1;
        // inComptonCounter++;
      }
    }
    
    TreeCluster->Fill();
    inCompton  = 0;
    inCry      = 0;
    inCryGood  = 0;
    inCryBad   = 0;
    outCryGood = 0;
    // out_u      = 0.;
    // out_v      = 0.;
    // out_w      = 0.;
  }
  std::cout << "Good events         = " << goodEventsAnalysis << std::endl;
  std::cout << "inCrystal           = " << inCrystal           << std::endl;
  std::cout << "inCrystalGood       = " << inCrystalGood       << std::endl;
  std::cout << "inCrystalBad        = " << inCrystalBad        << std::endl;
  std::cout << "outOfCrystalButGood = " << outOfCrystalButGood << std::endl;
  
  outputFile->cd();
  TreeCluster->Write();
  outputFile->Close();
  
  std::cout << "Results saved in file " << outputFileName << std::endl;
  
  return 0;
}





// feedback to user
void usage()
{
  std::cout << "\t" << "[-c|--calibration] <list> " 
            << "\t" << "[-f|--folder] <folder> " 
            << "\t" << "[-p|--prefix] <prefix> "
            << "\t" << "[-o|--output] <output>" << std::endl
            // << "\t" << "[--crystals]  <crystals>" << std::endl
            << "\t\t" << "<list>           - csv list of calibration.root files (outputs of ModuleCalibration) " << std::endl
            << "\t\t" << "<folder>         - path to folder were input files are located - default = ./ " << std::endl
            << "\t\t" << "<prefix>         - prefix of input TTree files " << std::endl
            << "\t\t" << "<output>         - output file name - default = outputBareboneFile.root " << std::endl
            << "\t\t" << "<cutgfile>       - cutg file name" << std::endl
            // << "\t\t" << "<crystals>         - number of crystals - default = outputBareboneFile.root " << std::endl
            << std::endl;
}
