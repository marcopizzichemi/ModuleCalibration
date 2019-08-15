// compile with
// g++ -o ../../build/testDOIres testDOIres.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

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
#include <sys/stat.h>
#include <dirent.h>

#include"./testDOIres.h"



void usage()
{
  std::cout << "\t\t" << "[-i|--input] <inputFile.root>  [-o|--output] <output> [-c|--calibration] <calibration.root> [-n|--crystal] <number> " << std::endl
            << "\t\t" << "<inputFile.root>                           - input file name"   << std::endl
            << "\t\t" << "<output>                                   - prefix of output files"   << std::endl
            << "\t\t" << "<calibration.root>                         - output file name"   << std::endl
            << "\t\t" << "<number>                                   - crystal nuber"   << std::endl
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
  std::string textFileName = "";
  std::string calibrationFileName = "";
  int selectedCrystal = -1;

  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "calibration", required_argument, 0, 0 },
      { "crystal", required_argument, 0, 0 },
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
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}
  //----------------------------------------------------------//


  //output file names
  textFileName = outputFileName + ".txt";
  outputFileName += ".root";
  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  //----------------------------------------------------------//

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
  Float_t      *timeStamp;
  TBranch      *bChainExtendedTimeTag;                               // branches for above data
  TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  TBranch      **bCharge;
  TBranch      **btimeStamp;
  TBranch      *b_zFromTag;
  charge = new UShort_t[numOfCh];
  timeStamp = new Float_t[numOfT];
  bCharge = new TBranch*[numOfCh];
  btimeStamp = new TBranch*[numOfT];
  Float_t      zFromTag;
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
  tree->SetBranchAddress("zFromTag", &zFromTag, &b_zFromTag);
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
                  0,                        // can be ignored here
                  0,                        // can be ignored here
                  0,                        // can be ignored here
                  1,                        // can be ignored here
                  1                         // can be ignored here
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
  // switch off crystals if they were not chosen (only one chosen allowed)
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
  //----------------------------------------------------------//

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
  TH2F *h2 = new TH2F("h2","z vs w",1000,0,1,100,0,15);
  for (long long int i=0;i<nevent;i++)
  {
    // std::cout << "Event " << i << std::endl;
    tree->GetEvent(i);              //read complete accepted event in memory
    Float_t w = calculateW(charge,crystal[id]);
    h2->Fill(w,zFromTag);

    //Progress...
    counter++;
    int perc = ((100*counter)/nevent);
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }
  std::cout << std::endl;

  // fit w histograms
  std::vector<float> x,y,ex,ey;
  int sigmaCounter=0;
  float avgSigma = 0.0;
  for (int i = 1 ; i < h2->GetNbinsY(); i++)
  {
    TH1D* pjx = (TH1D*) h2->ProjectionX("pjx",i,i,"");
    int entries = pjx->GetEntries();
    if(entries)
    {
      sigmaCounter++;
      TF1 *fGauss = new TF1("fGauss","gaus",0,1);
      pjx->Fit(fGauss,"Q");
      x.push_back(fGauss->GetParameter(1));
      ex.push_back(fGauss->GetParError(1));
      y.push_back(h2->GetYaxis()->GetBinCenter(i));
      float error = (h2->GetYaxis()->GetBinWidth(i)) / (3.4641016151) ;
      ey.push_back(error);
      avgSigma += fGauss->GetParameter(2);
    }
  }

  TGraphErrors *g = new TGraphErrors(x.size(),&x[0],&y[0],&ex[0],&ey[0]);
  g->SetMarkerStyle(21);
  g->SetMarkerColor(kBlue);
  TF1 *fLine = new TF1("fLine","pol1",0,1);
  g->Fit(fLine,"Q");
  // METHOD_1 to calc DOI RES:
  // Calculate uncertainty in w estimation
  avgSigma = avgSigma / sigmaCounter;
  float m = fLine->GetParameter(1);
  float calculatedDoiRes = fabs(avgSigma * m) ;

  //----------------------------------------------------------//
  //  METHOD_2 LOOP                                           //
  //----------------------------------------------------------//
  // build the histograms to calculate the DOI RES with method_2 i.e.
  // -> for each event, calculate the distance between real DOI given by tag and computed doi given linear interpolation
  // -> gauss fit to compute final values
  TH1F *hDeltaToLine = new TH1F("hDeltaToLine","hDeltaToLine",100,-10,10);
  counter = 0;
  for (long long int i=0;i<nevent;i++)
  {
    tree->GetEvent(i);              //read complete accepted event in memory

    Float_t w = calculateW(charge,crystal[id]);
    float doiFromLine        =  fLine->Eval(w);
    float deltaToLine        = zFromTag - doiFromLine;
    hDeltaToLine       ->Fill(deltaToLine);

    //Progress...
    counter++;
    int perc = ((100*counter)/nevent);
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }
  std::cout << std::endl;

  TF1 *fGaussLinear = new TF1("fGaussLinear","gaus",0,1);
  hDeltaToLine->Fit(fGaussLinear,"Q");

  // write in text file
  // i j crystalNumber DOIres1 DOIres2

  std::ofstream textfile (textFileName.c_str(), std::ofstream::out);
  textfile << (int) (crystal[id].number / 8) << " "
           << (int) (crystal[id].number % 8) << " "
           << crystal[id].number             << " "
           << calculatedDoiRes              *2.355 << " "
           << fGaussLinear->GetParameter(2) *2.355 << " "
           << std::endl;
  textfile.close();

  std::cout<< (int) (crystal[id].number / 8) << " "
           << (int) (crystal[id].number % 8) << " "
           << crystal[id].number             << " "
           << calculatedDoiRes              *2.355 << " "
           << fGaussLinear->GetParameter(2) *2.355 << " "
           << std::endl;

  outputFile->cd();
  hDeltaToLine->Write();
  outputFile->Close();

  return 0;
}
