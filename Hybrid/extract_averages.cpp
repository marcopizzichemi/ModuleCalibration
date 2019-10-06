// compile with
// g++ -o ../../build/extract_averages extract_averages.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

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

#include "../libraries/CrystalStructs.h"        // Crystal_t , detector_t
#include "../libraries/Calibration.h"           // readTaggingData , readCalibration , setWandZcuts
#include "../libraries/Utilities.h"             // read_directory , invert_a_matrix



void usage()
{
  std::cout << "\t\t" << "[-i|--input] <input.root>  [-o|--output] <output.txt> [-c|--crystal] number" << std::endl
            << "\t\t" << "<input.root>                - input file name"  << std::endl
            << "\t\t" << "<output.txt>                - output file name" << std::endl
            // << "\t\t" << "<number>                    - crystal number"   << std::endl
            << "\t\t" << std::endl;
}

main (int argc, char** argv)
{

  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }
  std::string inputFileName = "";
  std::string outputFileName = "";
  // std::string textFileName = "";
  // std::string calibrationFileName = "";
  int selectedCrystal = -1;
  // bool simulation = false;


  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      // { "calibration", required_argument, 0, 0 },
      // { "crystal", required_argument, 0, 0 },
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
		int c = getopt_long(argc, argv, "i:o:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			inputFileName = (char *)optarg;
    }
		else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
    // else if (c == 'c'){
    //   calibrationFileName = (char *)optarg;
    // }
    // else if (c == 'c'){
    //   selectedCrystal = atoi((char *)optarg);
    // }
    else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      outputFileName = (char *)optarg;
    }
    // else if (c == 0 && optionIndex == 2){
    //   calibrationFileName = (char *)optarg;
    // }
    // else if (c == 0 && optionIndex == 2){
    //   selectedCrystal = atoi((char *)optarg);
    // }
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


  TFile *_file0 = TFile::Open(inputFileName.c_str());

  _file0->cd("/Module 0.0/");

  std::stringstream snameCh;
  snameCh << ((TNamed*) gDirectory->Get("taggingPosition"))->GetTitle();
  float taggingPosition = atof(snameCh.str().c_str());

  gDirectory->cd("MPPC B2 - 0.0-1.1/Crystal 5");


  //get w histo
  TCanvas *c_w = (TCanvas*) gDirectory->Get("W histogram - Crystal 5");
  TH1F* wHisto = (TH1F*) c_w->GetPrimitive("W histogram - Crystal 5");
  // TF1* w_gauss = new TF1("w_gauss","gaus",0,1);
  // wHisto->Fit(w_gauss);

  // get basic ctr
  double meanW = wHisto->GetMean();


  TH1F* basicCtr = (TH1F*) gDirectory->Get("Basic CTR - Crystal 5");
  double meanDelta    = basicCtr->GetMean();
  double meanDeltaErr = basicCtr->GetMeanError();
  double meanDeltaRMS    = basicCtr->GetRMS();
  double meanDeltaRMSErr = basicCtr->GetRMSError();
  // TF1* bc_gauss = new TF1("bc_gauss","gaus",0,1);
  // basicCtr->Fit(bc_gauss);

  // get the delay values
  gDirectory->cd("Plots");
  float neigh_ch[8] = {16,19,23,24,25,26,28,29};

  std::vector<double> delay;
  std::vector<double> delayError;
  std::vector<double> lightSharing;
  std::vector<double> lightSharingErr;
  // std::vector<double> RMSdelay;
  // std::vector<double> RMSdelay;
  // std::vector<TH1F*> delayHistos;

  // std::cout << "w" << "\t"
  //           << "t0-tag " << "\t"
  //           << "s_t0-tag " << "\t";

  for(int i = 0 ; i < 8 ; i++)
  {
    std::stringstream sname;
    sname << "T_Channel_" << neigh_ch[i] << " - T_Crystal_5";
    TH1F* delayHisto = (TH1F*) gDirectory->Get(sname.str().c_str());

    delay.push_back(delayHisto->GetMean());
    delayError.push_back(delayHisto->GetMeanError());

  }

  gDirectory->cd("../RawCTRs");
  std::vector<double> delayRMS;
  std::vector<double> delayRMSerror;

  for(int i = 0 ; i < 8 ; i++)
  {
    std::stringstream sname;
    sname << "T_Channel_" << neigh_ch[i] << " - T_tag";
    TH1F* delayHisto = (TH1F*) gDirectory->Get(sname.str().c_str());

    delayRMS.push_back(delayHisto->GetRMS());
    delayRMSerror.push_back(delayHisto->GetRMSError());

  }

  gDirectory->cd("../LightSharingPlots");
  // TH1F* LScentral = (TH1F*) gDirectory->Get("Light collected in trigger crystal");
  lightSharing.push_back( ((TH1F*) gDirectory->Get("Light collected in trigger crystal"))->GetMean() );
  lightSharing.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC C2 for photoelectric event in crystal 5"))->GetMean() );
  lightSharing.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC C1 for photoelectric event in crystal 5"))->GetMean() );
  lightSharing.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC C3 for photoelectric event in crystal 5"))->GetMean() );
  lightSharing.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC B1 for photoelectric event in crystal 5"))->GetMean() );
  lightSharing.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC A1 for photoelectric event in crystal 5"))->GetMean() );
  lightSharing.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC A2 for photoelectric event in crystal 5"))->GetMean() );
  lightSharing.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC B3 for photoelectric event in crystal 5"))->GetMean() );
  lightSharing.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC A3 for photoelectric event in crystal 5"))->GetMean() );

  lightSharingErr.push_back( ((TH1F*) gDirectory->Get("Light collected in trigger crystal"))->GetMeanError() );
  lightSharingErr.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC C2 for photoelectric event in crystal 5"))->GetMeanError() );
  lightSharingErr.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC C1 for photoelectric event in crystal 5"))->GetMeanError() );
  lightSharingErr.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC C3 for photoelectric event in crystal 5"))->GetMeanError() );
  lightSharingErr.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC B1 for photoelectric event in crystal 5"))->GetMeanError() );
  lightSharingErr.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC A1 for photoelectric event in crystal 5"))->GetMeanError() );
  lightSharingErr.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC A2 for photoelectric event in crystal 5"))->GetMeanError() );
  lightSharingErr.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC B3 for photoelectric event in crystal 5"))->GetMeanError() );
  lightSharingErr.push_back( ((TH1F*) gDirectory->Get("Light collected on MPPC A3 for photoelectric event in crystal 5"))->GetMeanError() );

  int sumLight = std::accumulate(lightSharing.begin(), lightSharing.end(), 0);

  std::ofstream outFile;
  outFile.open(outputFileName.c_str(), std::ofstream::out);

  outFile   << taggingPosition << "\t"
            << meanW     << "\t"
            << meanDelta << "\t"
            << meanDeltaErr << "\t"
            << meanDeltaRMS << "\t"
            << meanDeltaRMSErr << "\t"
            << lightSharing[0] << "\t"
            << lightSharingErr[0] << "\t"
            << sumLight << "\t";
  for(int i = 0 ; i < 8 ; i++)
  {
    outFile   << neigh_ch[i]     << "\t"
              << delay[i]    << "\t"
              << delayError[i]    << "\t"
              << delayRMS[i]    << "\t"
              << delayRMSerror[i]    << "\t"
              << lightSharing[i+1] << "\t"
              << lightSharingErr[i+1] << "\t";
  }
  outFile << std::endl;
  outFile.close();
  // TFile *fHisto = new TFile("histos.root","RECREATE");
  // fHisto->cd();
  // wHisto->Write();
  // basicCtr->Write();
  //
  // fHisto->Close();



  return 0;

}
