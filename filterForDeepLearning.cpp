// compile with
// g++ -o ../build/filterForDeepLearning filterForDeepLearning.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer

// small program to take doi bench data and filter events in one cyrstal, producing an output for deeplearning study

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
#include <dirent.h>

// #include "./libraries/Calibration.h"        // readTaggingData , readCalibration , setWandZcuts

// typedef std::vector<std::string> stringvec;
// list files in directory
// taken from
// http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
void read_directory(const std::string& name, std::vector<std::string> &v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}




struct detector_t
{
  int digitizerChannel;
  float saturation;
  float pedestal;
};



struct data_t
{
  double x;
  double y;
};


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
  Short_t      *charge;
  Float_t      *timeStamp;
  TBranch      *bChainExtendedTimeTag;                               // branches for above data
  TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  TBranch      **bCharge;
  TBranch      **btimeStamp;
  charge = new Short_t[numOfCh];
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
  UShort_t  out_charge[64];
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

  //----------------------------------------------------------//
  //  Read and import from Calibration file                   //
  //----------------------------------------------------------//
  TFile *calibrationFile = new TFile(calibrationFileName.c_str());

  //prepare individual relevat cuts
  // TCut* CrystalCut = NULL;
  TCut* taggingPhotopeakCut = NULL;
  TCut* TriggerChannelCut = NULL;
  TCut* broadCut = NULL;
  TCut* CutNoise = NULL;
  std::vector<TCutG*> cutg;
  TCut* PhotopeakEnergyCut = NULL;

  //enter main folder
  calibrationFile->cd("Module 0.0");

  // get data for saturaiton correction
  std::vector<int> *pChannels;
  std::vector<float> *pSaturation;
  std::vector<float> *pPedestal;
  gDirectory->GetObject("channels",pChannels);
  gDirectory->GetObject("saturation",pSaturation);
  gDirectory->GetObject("pedestal",pPedestal);
  std::vector<int> DetChannels = pChannels[0];
  std::vector<float> saturation = pSaturation[0];
  std::vector<float> pedestal = pPedestal[0];
  std::vector<detector_t> detectorSaturation;
  for(unsigned int iSat = 0; iSat < DetChannels.size(); iSat++)
  {
    detector_t tempDetector;
    tempDetector.digitizerChannel = DetChannels[iSat];
    tempDetector.saturation = saturation[iSat];
    tempDetector.pedestal = pedestal[iSat];
    detectorSaturation.push_back(tempDetector);
  }



  taggingPhotopeakCut = (TCut*) gDirectory->Get("taggingPhotopeakCut");

  std::stringstream snameCh;
  snameCh << ((TNamed*) gDirectory->Get("taggingPosition"))->GetTitle();
  float taggingPosition = atof(snameCh.str().c_str());

  snameCh.str("");
  snameCh << ((TNamed*) gDirectory->Get("taggingCrystalChannel"))->GetTitle();
  int taggingCrystalChannel = atoi(snameCh.str().c_str());


  // std::cout << taggingCrystalChannel << std::endl;


  // list keys
  TList *listModule = gDirectory->GetListOfKeys(); // list of keys
  int nKeysMod = listModule->GetEntries();  // number of keys
  std::vector<std::string> keysModName;
  // fill a vector with the leaves names
  std::string mppc_prefix("MPPC");
  for(int i = 0 ; i < nKeysMod ; i++)
  {
    keysModName.push_back(listModule->At(i)->GetName()); // put keys names in a vector of strings
  }
  // get MPPC folders
  std::vector<std::string> MPPCfolders; // check keys starting by "MPPC"
  for(unsigned int i = 0 ; i < keysModName.size() ; i++)
  {
    if (!keysModName[i].compare(0, mppc_prefix.size(), mppc_prefix))
    {
      MPPCfolders.push_back(keysModName[i]); // put MPPC keys name in list
    }
  }

  // run on all the mppc
  for(unsigned int iMppc = 0 ; iMppc < MPPCfolders.size() ; iMppc++)
  {
    // std::cout << MPPCfolders[iMppc] << std::endl;
    gDirectory->cd(MPPCfolders[iMppc].c_str());
    TList *listMppc = gDirectory->GetListOfKeys();
    int nKeysMppc = listMppc->GetEntries();
    std::vector<std::string> keysMppcName;
    // fill a vector with the leaves names
    std::string crystal_prefix("Crystal");

    for(int i = 0 ; i < nKeysMppc ; i++){
      keysMppcName.push_back(listMppc->At(i)->GetName());
    }

    std::vector<std::string> CrystalFolders;

    for(unsigned int i = 0 ; i < keysMppcName.size() ; i++)
    {
      if (!keysMppcName[i].compare(0, crystal_prefix.size(), crystal_prefix))
      {
        CrystalFolders.push_back(keysMppcName[i]);
      }
    }

    // enter crystal folders
    for(unsigned int iCry = 0 ; iCry < CrystalFolders.size() ; iCry++)
    {
      gDirectory->cd(CrystalFolders[iCry].c_str());
      int CryNumber = atoi((CrystalFolders[iCry].substr(crystal_prefix.size()+1,CrystalFolders[iCry].size()-crystal_prefix.size()-1)).c_str());


      if(CryNumber == selectedCrystal)
      {


        TList *listCry = gDirectory->GetListOfKeys();
        int nKeysCry = listCry->GetEntries();
        std::vector<std::string> keysCryName;
        if(nKeysCry) //if directory not empty
        {
          // CrystalCut = (TCut*) gDirectory->Get("CrystalCut");
          TriggerChannelCut  = (TCut*) gDirectory->Get("TriggerChannelCut");
          broadCut           = (TCut*) gDirectory->Get("broadCut");
          CutNoise           = (TCut*) gDirectory->Get("CutNoise");
          PhotopeakEnergyCut = (TCut*) gDirectory->Get("PhotopeakEnergyCut");

          std::string cutG_prefix("cutg");
          for(int i = 0 ; i < nKeysCry ; i++)
          {
            keysCryName.push_back(listCry->At(i)->GetName());
          }
          for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
          {
            if(!keysCryName[i].compare(0,cutG_prefix.size(),cutG_prefix)) // find tcutgs
            {
              TCutG* cut = (TCutG*) gDirectory->Get(keysCryName[i].c_str());
              cutg.push_back(cut);
            }
          }
        }
      }
      gDirectory->cd("..");
    }
    gDirectory->cd("..");
  }
  //----------------------------------------------------------//



  //----------------------------------------------------------//
  //  check                                                   //
  //----------------------------------------------------------//
  if(taggingPhotopeakCut)
  {
    if(verbose)
    {
      std::cout << "---------------------------------------" << std::endl;
      std::cout << "taggingPhotopeakCut: " << std::endl;
      std::cout << taggingPhotopeakCut->GetTitle() << std::endl;
    }
  }
  else
  {
    std::cout << "taggingPhotopeakCut not found! Aborting." << std::endl;
    return 1;
  }

  if(TriggerChannelCut)
  {
    if(verbose)
    {
      std::cout << "---------------------------------------" << std::endl;
      std::cout << "TriggerChannelCut: " << std::endl;
      std::cout << TriggerChannelCut->GetTitle() << std::endl;
    }
  }
  else
  {
    std::cout << "TriggerChannelCut not found! Aborting." << std::endl;
    return 1;
  }

  if(broadCut)
  {
    if(verbose)
    {
      std::cout << "---------------------------------------" << std::endl;
      std::cout << "broadCut: " << std::endl;
      std::cout << broadCut->GetTitle() << std::endl;
    }
  }
  else
  {
    std::cout << "broadCut not found! Aborting." << std::endl;
    return 1;
  }

  if(CutNoise)
  {
    if(verbose)
    {
      std::cout << "---------------------------------------" << std::endl;
      std::cout << "CutNoise: " << std::endl;
      std::cout << CutNoise->GetTitle() << std::endl;
    }
  }
  else
  {
    std::cout << "CutNoise not found! Aborting." << std::endl;
    return 1;
  }

  if(PhotopeakEnergyCut)
  {
    if(verbose)
    {
      std::cout << "---------------------------------------" << std::endl;
      std::cout << "PhotopeakEnergyCut: " << std::endl;
      std::cout << PhotopeakEnergyCut->GetTitle() << std::endl;
    }
  }
  else
  {
    std::cout << "PhotopeakEnergyCut not found! Aborting." << std::endl;
    return 1;
  }
  if(cutg.size() == 2)
  {
    if(verbose)
    {
      std::cout << "---------------------------------------" << std::endl;
      std::cout << "cutg[0]: " << std::endl;
      cutg[0]->Print();
    }
    if(verbose)
    {
      std::cout << "---------------------------------------" << std::endl;
      std::cout << "cutg[1]: " << std::endl;
      cutg[1]->Print();
    }
  }
  else
  {
    std::cout << "CutGs not found! Aborting." << std::endl;
    return 1;
  }
  //----------------------------------------------------------//

  //----------------------------------------------------------//
  //  set global cuts                                         //
  //----------------------------------------------------------//

  // create formula for cutting dataset
  TCut globalCut ;
  if(UseTaggingPhotopeakCut)     globalCut += taggingPhotopeakCut->GetTitle();
  if(UseTriggerChannelCut)       globalCut += TriggerChannelCut->GetTitle();
  if(UseBroadCut)                globalCut += broadCut->GetTitle();
  if(UseCutNoise)                globalCut += CutNoise->GetTitle();
  if(UsePhotopeakEnergyCut)      globalCut += PhotopeakEnergyCut->GetTitle();
  if(UseCutgs)
  {
    globalCut += cutg[0]->GetName();
    globalCut += cutg[1]->GetName();
  }

  TTreeFormula* FormulaAnalysis = new TTreeFormula("FormulaAnalysis",globalCut,tree);
  TList* formulasAnalysis = new TList();
  formulasAnalysis->Add(FormulaAnalysis);


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

    if(FormulaAnalysis->EvalInstance())  //if in global cut of crystal
    {
      passCounter++;
      out_doiFromTag      = taggingPosition;
      out_ExtendedTimeTag = ChainExtendedTimeTag;
      out_DeltaTimeTag    = ChainDeltaTimeTag;
      for(int iCh = 0; iCh < numOfCh ; iCh++)
      {
        if(charge[iCh] > 0) // old readout could give negative integrals (why?) anyway filter them out
        {
          if(iCh == taggingCrystalChannel)
          {
            out_charge[iCh] = charge[iCh];
          }
          for(int iDet = 0 ; iDet < detectorSaturation.size(); iDet++)
          {
            // sat correction is
            // - q_max * ln( 1 - (ch - ped)/(q_max))
            // where
            // q_max    = saturation
            // ch       = charge not corrected
            // ped      = pedestal (not corrected)
            if(iCh == detectorSaturation[iDet].digitizerChannel)
            {
              out_charge[iCh] = -detectorSaturation[iDet].saturation * TMath::Log(1.0 - ( (charge[iCh] - detectorSaturation[iDet].pedestal ) / (detectorSaturation[iDet].saturation) ) );
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
