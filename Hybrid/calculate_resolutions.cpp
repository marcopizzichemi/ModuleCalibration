// compile with
// g++ -o ../../build/calculate_resolutions calculate_resolutions.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas



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
  std::cout << "\t\t" << "[-i|--input] <file_prefix>  [-o|--output] <output.txt> [-c|--calibration] <calibration.root> [-n|--crystal] number [--avg] avgFile.txt" << std::endl
            << "\t\t" << "<file_prefix>                                     - prefix of TTree files to analyze"   << std::endl
            << "\t\t" << "<output.root>                                     - output file name"   << std::endl
            << "\t\t" << "<calibration.root>                                - calibration file name"   << std::endl
            << "\t\t" << "<avgFile.txt>                                     - average file name"   << std::endl
            << "\t\t" << "<number>                                          - crystal number"   << std::endl
            << "\t\t" << "--linear                                          - use linear fits (default false if flag not given)"   << std::endl
            << "\t\t" << "--light                                           - use light sharing for weights (default false if flag not given)"   << std::endl

            << "\t\t" << std::endl;
}

struct correction_gr_t
{
  int timingChannel;
  int digitizerChannel;
  TGraphErrors *delay;
  TGraphErrors *rms;
  TGraphErrors *light;
  TGraphErrors *lightErr;
  TF1* delay_line;
  TF1* rms_line;
  TF1 *light_line;
  bool isMainChannel;
};

// class of input points from doi tag bench
class inputFile_t
{
public:
  double z;
  double w;
  double delta;
  double deltaErr;
  double deltaRMS;
  double deltaRMSerr;
  double lightCentral;
  double lightCentralErr;
  double lightSum;

  std::vector<int> ch;
  std::vector<double> delay;
  std::vector<double> delayErr;
  std::vector<double> delayRMS;
  std::vector<double> delayRMSerr;
  std::vector<double> lightSharing;
  std::vector<double> lightSharingErr;
  int nOfNeigh;
  inputFile_t(int a){ nOfNeigh = a;};
  inputFile_t(){};
  void clear()
  {
    ch.clear();
    delay.clear();
    delayErr.clear();
    delayRMS.clear();
    delayRMSerr.clear();
    lightSharing.clear();
    lightSharingErr.clear();
  };
  void setnOfNeigh(int a) {nOfNeigh = a;};
  friend std::istream& operator>>(std::istream& input, inputFile_t& s)
  {
    input >> s.z;
    input >> s.w;
    input >> s.delta;
    input >> s.deltaErr;
    input >> s.deltaRMS;
    input >> s.deltaRMSerr;
    input >> s.lightCentral;
    input >> s.lightCentralErr;
    input >> s.lightSum;
    for(int p = 0; p < s.nOfNeigh; p++)
    {
      double ch_v,delay_v,delayErr_v,delayRMS_v,delayRMSerror_v,ls_v,lse_v;
      input >> ch_v >> delay_v >> delayErr_v >> delayRMS_v >> delayRMSerror_v >> ls_v >> lse_v;
      s.ch.push_back(ch_v);
      s.delay.push_back(delay_v);
      s.delayErr.push_back(delayErr_v);
      s.delayRMS.push_back(delayRMS_v);
      s.delayRMSerr.push_back(delayRMSerror_v);
      s.lightSharing.push_back(ls_v);
      s.lightSharingErr.push_back(lse_v);
    }
    return input;
  }
};



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
  std::string textFileName = "";
  std::string calibrationFileName = "";
  std::string avgFile = "";
  int selectedCrystal = -1;
  bool simulation = false;
  bool linear = false;
  bool light = false;
  bool oldLight = false;
  int tChannel = 28;
  int dChannel = 28;


  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "calibration", required_argument, 0, 0 },
      { "crystal", required_argument, 0, 0 },
      { "avg", required_argument, 0, 0 },
      { "linear", no_argument, 0, 0 },
      { "light", no_argument, 0, 0 },
      { "oldLight", no_argument, 0, 0 },
      { "tChannel", no_argument, 0, 0 },
      { "dChannel", no_argument, 0, 0 },
      // { "sim",no_argument,0,0},
      // { "verbose", no_argument, 0, 0 },
      // { "tagCut", no_argument, 0, 0 },
      // { "triggerCut", no_argument, 0, 0 },
      // { "broadCut", no_argument, 0, 0 },
      // { "noiseCut", no_argument, 0, 0 },
      // { "photopeakCut", no_argument, 0, 0 },
      // { "cutgCut", no_argument, 0, 0 },
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
      avgFile = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 5){
      linear = true;
    }
    else if (c == 0 && optionIndex == 6){
      light = true;
    }
    else if (c == 0 && optionIndex == 7){
      oldLight = true;
    }
    else if (c == 0 && optionIndex == 8){
      tChannel = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      dChannel = atoi((char *)optarg);
    }
    // else if (c == 0 && optionIndex == 4){
    //   simulation = true;
    // }
    // else if (c == 0 && optionIndex == 4){
    //   verbose = true;
    // }
    // else if (c == 0 && optionIndex == 5){
    //   UseAllCuts = false;
    //   UseTaggingPhotopeakCut = true;
    // }
    // else if (c == 0 && optionIndex == 6){
    //   UseAllCuts = false;
    //   UseTriggerChannelCut = true;
    // }
    // else if (c == 0 && optionIndex == 7){
    //   UseAllCuts = false;
    //   UseBroadCut = true;
    // }
    // else if (c == 0 && optionIndex == 8){
    //   UseAllCuts = false;
    //   UseCutNoise = true;
    // }
    // else if (c == 0 && optionIndex == 9){
    //   UseAllCuts = false;
    //   UsePhotopeakEnergyCut = true;
    // }
    // else if (c == 0 && optionIndex == 10){
    //   UseAllCuts = false;
    //   UseCutgs = true;
    // }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}

  if(oldLight)
  {
    light = true;
  }

  textFileName = outputFileName + ".txt";
  outputFileName += ".root";
  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");

  // if(UseAllCuts)
  // {
  //   UseTaggingPhotopeakCut = true;
  //   UseTriggerChannelCut   = true;
  //   UseBroadCut            = true;
  //   UseCutNoise            = true;
  //   UsePhotopeakEnergyCut  = true;
  //   UseCutgs               = true;
  // }

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
  UShort_t      *charge;
  Short_t      *s_charge;
  Float_t      *timeStamp;
  TBranch      *bChainExtendedTimeTag;                               // branches for above data
  TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  TBranch      **bCharge;
  TBranch      **btimeStamp;
  TBranch      *b_zFromTag;
  charge = new UShort_t[numOfCh];
  s_charge = new Short_t[numOfCh];
  timeStamp = new Float_t[numOfT];
  bCharge = new TBranch*[numOfCh];
  btimeStamp = new TBranch*[numOfT];
  Float_t      zFromTag;

  // sim part
  Float_t RealX;
  Float_t RealY;
  Float_t RealZ;
  Short_t CrystalsHit;
  Short_t NumbOfInteractions;
  TBranch *bRealX;
  TBranch *bRealY;
  TBranch *bRealZ;
  TBranch *bCrystalsHit;
  TBranch *bNumbOfInteractions;



  // set branches for reading the input files
  tree->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag, &bChainExtendedTimeTag);
  tree->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag, &bChainDeltaTimeTag);
  // if(simulation)
  // {
  //   for (int i = 0 ; i < detector_channels.size() ; i++)
  //   {
  //     std::stringstream sname;
  //     sname << "ch" << detector_channels[i];
  //     tree->SetBranchAddress(sname.str().c_str(),&s_charge[detector_channels[i]],&bCharge[detector_channels[i]]);
  //     sname.str("");
  //   }
  // }
  // else
  // {
    for (int i = 0 ; i < detector_channels.size() ; i++)
    {
      std::stringstream sname;
      sname << "ch" << detector_channels[i];
      tree->SetBranchAddress(sname.str().c_str(),&charge[detector_channels[i]],&bCharge[detector_channels[i]]);
      sname.str("");
    }
  // }


  // if(simulation)
  // {
  //   tree->SetBranchAddress("RealX", &RealX, &bRealX);
  //   tree->SetBranchAddress("RealY", &RealY, &bRealY);
  //   tree->SetBranchAddress("RealZ", &RealZ, &bRealZ);
  //   tree->SetBranchAddress("CrystalsHit",&CrystalsHit, &bCrystalsHit);
  //   tree->SetBranchAddress("NumbOfInteractions",&NumbOfInteractions, &bNumbOfInteractions);
  // }
  // else
  // {
    tree->SetBranchAddress("zFromTag", &zFromTag, &b_zFromTag);
  // }


  for (int i = 0 ; i < t_channels.size() ; i++)
  {
    std::stringstream sname;
    sname << "t" << t_channels[i];
    tree->SetBranchAddress(sname.str().c_str(),&timeStamp[t_channels[i]],&btimeStamp[t_channels[i]]);
    sname.str("");
  }
  //----------------------------------------------------------//


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
                  false,                        // include TriggerChannelCuthannel cut in crystalCut
                  false,                        // include broadCut in crystalCut
                  false,                               // include CutNoise in crystalCut
                  true,                          // include PhotopeakEnergyCut in crystalCut
                  true                                   // include CutGs in crystalCut
                 );

  // optionally set w and z limits, and write values into crystal struct
  // setWandZcuts(crystal);
  tree->SetNotify(formulasAnalysis);

  // list the crystals with calibration data found
  std::cout << "Calibration data found for crystals: " << std::endl;
  for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  {
    if(crystal[i].accepted)
    {
      std::cout << crystal[i].number << std::endl;
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

  long long int nevent = tree->GetEntries();
  std::cout << "Total number of events = " << nevent << std::endl;
  long int counter = 0;


  //----------------------------------------------------------//
  //  LINEAR REG LOOP                                         //
  //----------------------------------------------------------//
  // fill scatter plot
  // TH2F *h2 = new TH2F("h2","z vs w",1000,0,1,10000,0,15);
  // TH2F *h2_delta = new TH2F("h2_delta","(t_cry-t_tag) vs w",1000,0,1,10000,0,15);
  // TH2F *h2_delay = new TH2F("h2_delay","(t_neigh - t_cry) vs w",1000,0,1,10000,0,15);
  // for (long long int i=0;i<nevent;i++)
  // {
  //   // std::cout << "Event " << i << std::endl;
  //   tree->GetEvent(i);              //read complete accepted event in memory
  //
  //     if(crystal[id].accepted)
  //     {
  //
  //         Float_t w = calculateFloodZ_withoutCorrectingForSaturation(charge,crystal[id]);
  //         h2->Fill(w,zFromTag);
  //     }
  //
  //
  //   //Progress
  //   counter++;
  //   int perc = ((100*counter)/nevent);
  //   if( (perc % 10) == 0 )
  //   {
  //     std::cout << "\r";
  //     std::cout << perc << "% done... ";
  //   }
  // }
  // std::cout << std::endl;

  // open avg file and import data
  std::vector<inputFile_t> inputData;
  inputFile_t tempInput(8);
  std::ifstream inputFile;
  inputFile.open(avgFile.c_str(), std::ios::in);
  while(inputFile >> tempInput)
  {
    inputData.push_back(tempInput);
    tempInput.clear();
  }



  // int sigmaCounter=0;
  // float avgSigma = 0.0;
  // float calculatedDoiRes = 0.0;

  // for(unsigned int i = 0 ; i < inputData.size(); i++ )
  // {
  //   std::cout << inputData[i].z << "\t";
  //   std::cout << inputData[i].w << "\t";
  //   std::cout << inputData[i].delta << "\t";
  //   std::cout << inputData[i].deltaErr << "\t";
  //   std::cout << inputData[i].deltaRMS << "\t";
  //   std::cout << inputData[i].deltaRMSerr << "\t";
  //   std::cout << inputData[i].lightCentral << "\t";
  //   std::cout << inputData[i].lightCentralErr << "\t";
  //   std::cout << inputData[i].lightSum << "\t";
  //   for(int k = 0; k < inputData[i].nOfNeigh ; k++)
  //   {
  //     std::cout << inputData[i].ch[k] << "\t";
  //     std::cout << inputData[i].delay[k] << "\t";
  //     std::cout << inputData[i].delayErr[k] << "\t";
  //     std::cout << inputData[i].delayRMS[k] << "\t";
  //     std::cout << inputData[i].delayRMSerr[k] << "\t";
  //     std::cout << inputData[i].lightSharing[k] << "\t";
  //     std::cout << inputData[i].lightSharingErr[k] << "\t";
  //   }
  //   std::cout << std::endl;
  // }
  //
  // std::cout << std::endl;

  // for(int k = 0; k < inputData[0].nOfNeigh ; k++)
  // {
  //   for(unsigned int i = 0 ; i < inputData.size(); i++ )
  //   {
  //     std::cout << inputData[i].ch[k] << "\t";
  //     std::cout << inputData[i].w << "\t";
  //     std::cout << inputData[i].delay[k] << std::endl;
  //   }
  //   std::cout << std::endl;
  //
  // }
  // std::cout << std::endl;


  // we need graphs for 2 corrections
  // 1. FULL -> the correction to use all the light collected from 9 mppcs, so that we get the best estimator for t0
  // 2. CENTRAL -> the correction of doi bias
  // we start from 2.
  // need to get the so called "delta" graphs, i.e. the (t_cry - t_tag) as a function of w.
  // this is store in the txt file and not in inputData vector, in delta, with it's error deltaErr
  //
  // then, for 1. we need the delay graphs, i.e. for each neighbour the graph of (t_neigh - t_cry)
  // as a function of w. These are in delay data. Then we also need the estimation of sigma to calc
  // the weights, both for the central timestamp and for the neighbours.
  // For central, the error as a function of w is in deltaRMS
  // For the neighbours, it's in delayRMS
  // both have their errors, for doing error bars for fitting if needed

  // we also fit all graphs with lines, to be able to fully correct data
  // line forced to be horizontal for RMS graphs -> average

  // get the delta graph and the deltaRMS graph
  std::vector<float> w,z,ew,ez;
  std::vector<float> delta,deltaErr,deltaRMS,deltaRMSerr,lightCentral,lightCentralErr,lightSum,lightSumErr;

  for(unsigned int i = 0 ; i < inputData.size(); i++ )
  {
    // z vs w
    w.push_back(inputData[i].w);
    z.push_back(inputData[i].z);
    ew.push_back(0);
    ez.push_back(0);
    // delta vs w
    delta.push_back(inputData[i].delta);
    deltaErr.push_back(inputData[i].deltaErr);
    // delta RMS vs w
    deltaRMS.push_back(inputData[i].deltaRMS);
    deltaRMSerr.push_back(inputData[i].deltaRMSerr);
    //light data
    lightCentral.push_back(inputData[i].lightCentral);
    lightCentralErr.push_back(inputData[i].lightCentralErr);
    lightSum.push_back(inputData[i].lightSum);
    lightSumErr.push_back(0);
  }

  // build a graph of z(w)
  TGraphErrors *g_z_w ;
  g_z_w = new TGraphErrors(w.size(),&w[0],&z[0],&ew[0],&ew[0]);
  g_z_w->SetName("z(w)");
  g_z_w->GetXaxis()->SetTitle("w");
  g_z_w->GetYaxis()->SetTitle("z [mm]");
  g_z_w->SetMarkerStyle(21);
  g_z_w->SetMarkerColor(kBlue);
  // fit it with a line
  TF1 *fLine_zw = new TF1("fLine_gzw","pol1",0,1);
  g_z_w->Fit(fLine_zw,"Q");

  // build a graph of w(z)
  TGraphErrors *g_w_z = new TGraphErrors(w.size(),&z[0],&w[0],&ez[0],&ew[0]);
  g_w_z->SetName("w(z)");
  g_w_z->GetYaxis()->SetTitle("w");
  g_w_z->GetXaxis()->SetTitle("z [mm]");
  g_w_z->SetMarkerStyle(21);
  g_w_z->SetMarkerColor(kBlue);

  TF1 *fLine_wz = new TF1("fLine_wz","pol1",0,15);
  g_w_z->Fit(fLine_wz,"Q");

  // now built the delta graph
  TGraphErrors *gDelta = new TGraphErrors(w.size(),&w[0],&delta[0],&ew[0],&deltaErr[0]);
  gDelta->SetName("gDelta");
  gDelta->GetXaxis()->SetTitle("w");
  gDelta->GetYaxis()->SetTitle("t_cry - t_tag [s]");
  gDelta->SetMarkerStyle(21);
  gDelta->SetMarkerColor(kBlue);

  TF1 *fLine_delta = new TF1("fLine_delta","pol1",0,1);
  gDelta->Fit(fLine_delta,"Q");

  // now built the delta RMS graph
  TGraphErrors *gDeltaRMS = new TGraphErrors(w.size(),&w[0],&deltaRMS[0],&ew[0],&deltaRMSerr[0]);
  std::stringstream sname;
  sname << "gDeltaRMS t"  << tChannel;
  gDeltaRMS->SetName(sname.str().c_str());
  gDeltaRMS->GetXaxis()->SetTitle("w");
  gDeltaRMS->GetYaxis()->SetTitle("RMS [s]");
  gDeltaRMS->SetMarkerStyle(21);
  gDeltaRMS->SetMarkerColor(kBlue);
  sname.str("");

  TF1 *fLine_deltaRMS = new TF1("fLine_deltaRMS","pol0",0,1);
  gDeltaRMS->Fit(fLine_deltaRMS,"Q");

  // build a light graph
  TGraphErrors *gLight = new TGraphErrors(w.size(),&w[0],&lightCentral[0],&ew[0],&lightCentralErr[0]);
  sname << "gLight c"  << dChannel;
  gLight->SetName(sname.str().c_str());
  gLight->GetXaxis()->SetTitle("w");
  gLight->GetYaxis()->SetTitle("Light [ADC ch]");
  gLight->SetMarkerStyle(21);
  gLight->SetMarkerColor(kBlue);
  sname.str("");

  TF1 *fLine_Light = new TF1("fLine_Light","pol1",0,1);
  gLight->Fit(fLine_Light,"Q");

  // build a light err graph
  TGraphErrors *gLightErr = new TGraphErrors(w.size(),&w[0],&lightCentralErr[0],&ew[0],&ew[0]);
  sname << "gLight err ch"  << dChannel;
  gLightErr->SetName(sname.str().c_str());
  gLightErr->GetXaxis()->SetTitle("w");
  gLightErr->GetYaxis()->SetTitle("Light [ADC ch]");
  gLightErr->SetMarkerStyle(21);
  gLightErr->SetMarkerColor(kBlue);
  sname.str("");

  // build a light sum graph
  TGraphErrors *gLightSum = new TGraphErrors(w.size(),&w[0],&lightSum[0],&ew[0],&lightSumErr[0]);
  gLightSum->SetName("gLightSum");
  gLightSum->GetXaxis()->SetTitle("w");
  gLightSum->GetYaxis()->SetTitle("Light [ADC ch]");
  gLightSum->SetMarkerStyle(21);
  gLightSum->SetMarkerColor(kBlue);

  // TF1 *fLine_Light = new TF1("fLine_Light","pol1",0,1);
  // gLightSum->Fit(fLine_Light,"Q");


  // now the delay graphs
  // we store them in a struct, to avoid mess
  // to be more general, the struct includes also the central mppc.
  // obvioulsy the central mppc has no delay wrt to itself, but
  // we can store its RMS graph here

  // create the struct
  std::vector<correction_gr_t> corr_gr;


  // fill the first one, i.e. the central mppc
  correction_gr_t tempCorr;

  tempCorr.timingChannel = tChannel;
  tempCorr.digitizerChannel = dChannel;
  tempCorr.isMainChannel = true;
  tempCorr.delay = NULL;
  tempCorr.rms = gDeltaRMS;
  tempCorr.light = gLight;
  tempCorr.lightErr = gLightErr;
  tempCorr.delay_line = NULL;
  tempCorr.rms_line = fLine_deltaRMS;
  tempCorr.light_line = fLine_Light;
  corr_gr.push_back(tempCorr);


  // for(unsigned int i = 0 ; i < inputData.size(); i++ )
  // {
  //   for(int k = 0; k < inputData[i].nOfNeigh ; k++)
  //   {
  //     corr_gr[k+1]...
  //     // std::cout << inputData[i].ch[k] << "\t";
  // //     std::cout << inputData[i].delay[k] << "\t";
  // //     std::cout << inputData[i].delayErr[k] << "\t";
  // //     std::cout << inputData[i].delayRMS[k] << "\t";
  // //     std::cout << inputData[i].delayRMSerror[k] << "\t";
  //   }
  // }

  // now fill the others
  for(int k = 0; k < inputData[0].nOfNeigh ; k++)
  {
    // prepare the TGraphErrors and the vectors
    TGraphErrors *gDelay;
    TGraphErrors *gDelayRMS;
    TGraphErrors *gLightSharing;
    TGraphErrors *gLightSharingErr;
    std::vector<float> delay,delayErr,delayRMS,delayRMSerr,lSharing,lSharingErr;
    int channel = inputData[0].ch[k]; // take channel from first row, they are all the same in a column

    for(unsigned int i = 0 ; i < inputData.size(); i++ ) // run on rows (i.e. on doi scan points)
    {
      // std::cout << channel << "\t";
      // std::cout << inputData[i].w << "\t";
      // std::cout << inputData[i].delay[k] << std::endl;
      delay.push_back(inputData[i].delay[k]);
      delayErr.push_back(inputData[i].delayErr[k]);
      delayRMS.push_back(inputData[i].delayRMS[k]);
      delayRMSerr.push_back(inputData[i].delayRMSerr[k]);
      lSharing.push_back(inputData[i].lightSharing[k]);
      lSharingErr.push_back(inputData[i].lightSharingErr[k]);
    }

    gDelay = new TGraphErrors(w.size(),&w[0],&delay[0],&ew[0],&delayErr[0]);
    sname << "gDelay t" << channel;
    gDelay->SetName(sname.str().c_str());
    sname.str("");
    gDelay->GetXaxis()->SetTitle("w");
    sname << "t" << channel << " - t_cry [s]";
    gDelay->GetYaxis()->SetTitle(sname.str().c_str());
    sname.str("");
    gDelay->SetMarkerStyle(21);
    gDelay->SetMarkerColor(kBlue);
    TF1 *fLine_delay = new TF1("fLine_delay","pol1",0,1);
    gDelay->Fit(fLine_delay,"Q");

    gDelayRMS = new TGraphErrors(w.size(),&w[0],&delayRMS[0],&ew[0],&delayRMSerr[0]);
    sname << "gDelayRMS t" << channel;
    gDelayRMS->SetName(sname.str().c_str());
    sname.str("");
    gDelayRMS->GetXaxis()->SetTitle("w");
    sname << "RMS (t" << channel << " - t_cry) [s]";
    gDelayRMS->GetYaxis()->SetTitle(sname.str().c_str());
    sname.str("");
    gDelayRMS->SetMarkerStyle(21);
    gDelayRMS->SetMarkerColor(kBlue);
    TF1 *fLine_delayRMS = new TF1("fLine_delayRMS","pol0",0,1);
    gDelayRMS->Fit(fLine_delayRMS,"Q");

    gLightSharing = new TGraphErrors(w.size(),&w[0],&lSharing[0],&ew[0],&lSharingErr[0]);
    sname << "gLight ch" << channel;
    gLightSharing->SetName(sname.str().c_str());
    sname.str("");
    gLightSharing->GetXaxis()->SetTitle("w");
    sname << "Light [ADC ch]";
    gLightSharing->GetYaxis()->SetTitle(sname.str().c_str());
    sname.str("");
    gLightSharing->SetMarkerStyle(21);
    gLightSharing->SetMarkerColor(kBlue);
    TF1 *fLine_LightSharing = new TF1("fLine_LightSharing","pol1",0,1);
    gLightSharing->Fit(fLine_LightSharing,"Q");

    gLightSharingErr = new TGraphErrors(w.size(),&w[0],&lSharingErr[0],&ew[0],&ew[0]);
    sname << "gLightErr ch" << channel;
    gLightSharingErr->SetName(sname.str().c_str());
    sname.str("");
    gLightSharingErr->GetXaxis()->SetTitle("w");
    sname << "Light [ADC ch]";
    gLightSharingErr->GetYaxis()->SetTitle(sname.str().c_str());
    sname.str("");
    gLightSharingErr->SetMarkerStyle(21);
    gLightSharingErr->SetMarkerColor(kBlue);

    correction_gr_t tempCorrNeigh;
    tempCorrNeigh.timingChannel    = channel;            // set timing channel
    tempCorrNeigh.digitizerChannel = channel;            // set digitizerChannel channel
    tempCorrNeigh.isMainChannel    = false;              // set flag to not main
    tempCorrNeigh.delay            = gDelay;             // set the delay gr
    tempCorrNeigh.rms              = gDelayRMS;          // set the rms gr
    tempCorrNeigh.light            = gLightSharing;      // set the light gr
    tempCorrNeigh.lightErr         = gLightSharingErr;   // set the light gr
    tempCorrNeigh.delay_line       = fLine_delay;        // set the delay line
    tempCorrNeigh.rms_line         = fLine_delayRMS;     // set the rms line
    tempCorrNeigh.light_line       = fLine_LightSharing; // set the light line
    corr_gr.push_back(tempCorrNeigh);
  }









  //----------------------------------------------------------//
  //  RES LOOP                                                //
  //----------------------------------------------------------//
  // DOI
  // build the histograms to calculate the DOI RES with second and third method, i.e.
  // -> for each event, calculate the distance between real DOI given by tag and the one derived using the linear interpolation
  // use both RMS and gauss fit to compute final values
  TH1F *hDeltaToLine = new TH1F("hDeltaToLine","hDeltaToLine",100,-10,10);
  TH1F* rawCTR_H = new TH1F("rawCTR_H","rawCTR_H",100,-1e-9,1e-9);
  TH1F* baseCTR_H = new TH1F("baseCTR_H","baseCTR_H",100,-1e-9,1e-9);
  TH1F* centralCorrH = new TH1F("centralCorrH","centralCorrH",100,-1e-9,1e-9);
  TH1F* bestT0corrH = new TH1F("bestT0corrH","bestT0corrH",100,-1e-9,1e-9);
  TH1F* fullCorrH = new TH1F("fullCorrH","fullCorrH",100,-1e-9,1e-9);

  //MAIN LOOP
  // long int passCounter = 0;
  tree->SetNotify(formulasAnalysis);
  counter = 0;
  for (long long int i=0;i<nevent;i++)
  {
    // std::cout << "Event " << i << std::endl;
    tree->GetEvent(i);              //read complete accepted event in

    // calc w
    Float_t FloodZ = calculateFloodZ_withoutCorrectingForSaturation(charge,crystal[id]);

    // fill one histo with no constrains
    rawCTR_H->Fill(timeStamp[crystal[id].timingChannel] - timeStamp[crystal[id].taggingCrystalTimingChannel]);

    if(crystal[id].accepted)
    {
      bool noZeroes = true;
      // first check if there are no zeroes
      // this means that events where the DAQ failed to compute the timestamp of
      // one of the 9 detectors are ignored
      // can be improved by improving DAQ!
      for(unsigned int iCorr = 0; iCorr < corr_gr.size();iCorr ++) // run on all relevant channels
      {
        if(timeStamp[corr_gr[iCorr].timingChannel] == 0)
        {
          noZeroes = false;
        }
      }
      // at the same time, sometimes the corrected elements are beyond crazy,
      // probably because the DAQ was not aligned. ignored also these events
      // i.e. if one of the channels is triggering extremely far from the others (more than 15 ns)
      for(unsigned int iCorr = 0; iCorr < corr_gr.size();iCorr ++) // run on all relevant channels
      {
        // calc the delay
        float delay = 0.0;
        if(corr_gr[iCorr].isMainChannel)
        {
          delay = 0.0;
        }
        else
        {
          if(linear)
          {
            delay = corr_gr[iCorr].delay_line->Eval(FloodZ);
          }
          else
          {
            delay = corr_gr[iCorr].delay->Eval(FloodZ);
          }
        }
        float delta = timeStamp[corr_gr[iCorr].timingChannel] - timeStamp[crystal[id].taggingCrystalTimingChannel] - delay;
        if(delta <=  -15e-9 || delta >= 15e-9 )
        {
          noZeroes = false;
        }
      }
      if(noZeroes)
      {

        // fill base ctr histo
        baseCTR_H->Fill(timeStamp[crystal[id].timingChannel] - timeStamp[crystal[id].taggingCrystalTimingChannel]);

        // calc reconstructed doi
        float doiFromLine;
        if(linear)
        {
          doiFromLine =  fLine_zw->Eval(FloodZ);
        }
        else
        {
          doiFromLine = g_z_w->Eval(FloodZ);
        }

        // calc distance to "real" z (i.e. to tagging setup position)
        float deltaToLine        = zFromTag - doiFromLine;
        hDeltaToLine       ->Fill(deltaToLine); // fill histogram to calc resolution


        // calc doi t correction (CENTRAL)
        float correction;
        if(linear)
        {
          correction = fLine_delta->Eval(fLine_wz->Eval(7.5)) - fLine_delta->Eval(FloodZ);
        }
        else
        {
          correction = gDelta->Eval(g_w_z->Eval(7.5)) - gDelta->Eval(FloodZ);
        }
         // = fLine_delta->Eval(fLine_wz->Eval(7.5)) - fLine_delta->Eval(FloodZ);

        // calc corrected t stamp with only central correction
        float correctedStamp = (timeStamp[crystal[id].timingChannel] + (correction)) - timeStamp[crystal[id].taggingCrystalTimingChannel];
        centralCorrH->Fill(correctedStamp); // fill histo

        // calc best estimator of t0 (first part of FULL)
        float averageTimeStamp = 0.0;
        float totalWeight = 0.0;


        for(unsigned int iCorr = 0; iCorr < corr_gr.size();iCorr ++) // run on all relevant channels
        {

          // calc the delay
          float delay = 0.0;
          if(corr_gr[iCorr].isMainChannel)
          {
            delay = 0.0;
          }
          else
          {
            if(linear)
            {
              delay = corr_gr[iCorr].delay_line->Eval(FloodZ);
            }
            else
            {
              delay = corr_gr[iCorr].delay->Eval(FloodZ);
            }
          }


          float delta = timeStamp[corr_gr[iCorr].timingChannel] - timeStamp[crystal[id].taggingCrystalTimingChannel] - delay;

          float weight;
          if(light)
          {
            if(oldLight)
            {
              float lValue;
              if(linear)
              {
                lValue = corr_gr[iCorr].light_line->Eval(FloodZ);
              }
              else
              {
                lValue = corr_gr[iCorr].light->Eval(FloodZ);
              }
              float lErr =1;
              // float lErr = corr_gr[iCorr].lightErr->Eval(FloodZ);

              float lSigma = -0.5 * pow(lValue,-1.5)*lErr;

              weight = (1.0)/( TMath::Power(lSigma,2) );
            }
            else
            {
              // new correction for light
              // the weight is calculated on the assumption that sigma(t) goes with 1/sqrt(L), where
              // L is the light on the detector. Light is taken as the value corrected by saturation
              // For main channel, the value whose sigma needs to be calculated is
              // t_0 - t_ref
              // so the sigma is
              // sqrt( (1/L_0) + (1/L_ref) )
              // for the other channels is
              // t_N - t_ref - delay
              // so the sigma(delay) needs to be evaluated in terms of light
              // we do it like this (remember delay = t_N -t_0 )
              // sigma(delay) = sqrt( (1/avg(L_N)) + (1/avg(L_0)) )

              // so the formula is

              //sigma = sqrt( A + B + C + D ), with
              // A = 1/L_N
              // B = 1/L_ref
              // C = 1/avg(L_N)
              // D = 1/avg(L_0)

              // but C = D = 0 for central channels (delay = 0)

              // so first calc 1/L values
              float a = 1.0/(charge[corr_gr[iCorr].digitizerChannel]);
              float b = 1.0/(charge[crystal[id].taggingCrystalChannel]);
              float c = 0;
              float d = 0;

              if(corr_gr[iCorr].isMainChannel)
              {
                // nothing, c and d stay = 0
              }
              else
              {
                if(linear)
                {
                  c = 1.0/(corr_gr[iCorr].light_line->Eval(FloodZ));
                  d = 1.0/(corr_gr[0].light_line->Eval(FloodZ)); // main channel always first. This is seriously dangerous coding... :D
                }
                else
                {
                  c = 1.0/(corr_gr[iCorr].light->Eval(FloodZ));
                  d = 1.0/(corr_gr[0].light->Eval(FloodZ)); // main channel always first. This is seriously dangerous coding... :D
                }
              }
              // std::cout << a << " "
              //           << b << " "
              //           << c << " "
              //           << d << "\n";






              float lSigma = sqrt( a + b + c + d);



              // float lErr =1;
              // float lErr = corr_gr[iCorr].lightErr->Eval(FloodZ);

              // float lSigma = -0.5 * pow(lValue,-1.5)*lErr;

              weight = (1.0)/( TMath::Power(lSigma,2) );
            }
          }
          else
          {
            if(linear)
            {
              weight = (1.0)/( TMath::Power(corr_gr[iCorr].rms_line->Eval(FloodZ),2) );
            }
            else
            {
              weight = (1.0)/( TMath::Power(corr_gr[iCorr].rms->Eval(FloodZ),2) );
            }
          }

          totalWeight += weight;
          averageTimeStamp += delta*weight;
        }
        averageTimeStamp = averageTimeStamp/totalWeight;
        // fill histogram with just best estimator of t0
        bestT0corrH->Fill(averageTimeStamp); // fill histo
        float fullCorrectionCTR = averageTimeStamp + correction;

        // fill histogram with full correction
        fullCorrH->Fill(fullCorrectionCTR); // fill histo

      }

    }

    //Progress
    counter++;
    int perc = ((100*counter)/nevent);
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }
  std::cout << std::endl;

  // TF1 *fGaussCalib = new TF1("fGaussCalib","gaus",0,1);
  // hDeltaToCalibration->Fit(fGaussCalib,"Q");

  TF1 *fGaussLinear = new TF1("fGaussLinear","gaus",0,1);
  hDeltaToLine->Fit(fGaussLinear,"Q");

  // fit the ctr plots
  TF1 *fGaussRaw = new TF1("fGaussRaw","gaus",-1e-9,1e-9);
  TF1 *fGaussBase = new TF1("fGaussBase","gaus",-1e-9,1e-9);
  TF1 *fGaussCentral = new TF1("fGaussCentral","gaus",-1e-9,1e-9);
  TF1 *fGaussBestT0 = new TF1("fGaussBestT0","gaus",-1e-9,1e-9);
  TF1 *fGaussFull = new TF1("fGaussFull","gaus",-1e-9,1e-9);

  rawCTR_H->Fit(fGaussRaw,"Q");
  baseCTR_H->Fit(fGaussBase,"Q");
  centralCorrH->Fit(fGaussCentral,"Q");
  bestT0corrH->Fit(fGaussBestT0,"Q");
  fullCorrH->Fit(fGaussFull,"Q");

  float doiResFWHM = fGaussLinear->GetParameter(2) *2.355;
  float rawCTR    = sqrt(2)*sqrt(pow(fGaussRaw    ->GetParameter(2)*2.355,2)-pow(90e-12,2));
  float baseCTR   = sqrt(2)*sqrt(pow(fGaussBase   ->GetParameter(2)*2.355,2)-pow(90e-12,2));
  float centraCTR = sqrt(2)*sqrt(pow(fGaussCentral->GetParameter(2)*2.355,2)-pow(90e-12,2));
  float bestT0CTR = sqrt(2)*sqrt(pow(fGaussBestT0 ->GetParameter(2)*2.355,2)-pow(90e-12,2));
  float fullCTR   = sqrt(2)*sqrt(pow(fGaussFull   ->GetParameter(2)*2.355,2)-pow(90e-12,2));

  // write in text file
  // crystalNumber DOIres1 DOIres2(rms) DOIres2(fit) DOIres3(rms) DOIres3(fit)

  std::ofstream textfile (textFileName.c_str(), std::ofstream::out);
  textfile << "# crystal doiRES rawCTR baseCTR centraCTR bestT0CTR fullCTR (all FWHM)" << std::endl;
  textfile  << crystal[id].number << " "
            << doiResFWHM << " "
            << rawCTR     << " "
            << baseCTR    << " "
            << centraCTR  << " "
            << bestT0CTR  << " "
            << fullCTR    << " "
            << std::endl;
  textfile.close();

  std::cout << "# crystal doiRES rawCTR baseCTR centraCTR bestT0CTR fullCTR (all FWHM)" << std::endl;
  std::cout << crystal[id].number << " "
            << doiResFWHM << " "
            << rawCTR     << " "
            << baseCTR   << " "
            << centraCTR << " "
            << bestT0CTR << " "
            << fullCTR   << " "
            << std::endl;

  outputFile->cd();
  // hDeltaToCalibration->Write();
  hDeltaToLine->Write();
  baseCTR_H->Write();
  centralCorrH->Write();
  bestT0corrH->Write();
  fullCorrH->Write();
  // h2->Write();
  // if(!simulation)
  // {
  g_z_w->Write();
  g_w_z->Write();
  gDelta->Write();
  gLightSum->Write();
  // gDeltaRMS->Write();

  for (unsigned int i = 0; i < corr_gr.size(); i++)
  {
    if (i!=0) corr_gr[i].delay->Write();
    corr_gr[i].rms->Write();
    corr_gr[i].light->Write();
    corr_gr[i].lightErr->Write();
  }
  // }

  outputFile->Close();

  return 0;
}
