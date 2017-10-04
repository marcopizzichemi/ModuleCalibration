// compile with
// g++ -o ../build/doiResolutionWithCalibration doiResolutionWithCalibration.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer

// small program to extrac config file and pc info from the results of ModuleCalibration.
// run extractConfiguration without arguments for usage info

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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <getopt.h>
#include <algorithm>    // std::sort


struct Crystal_t
{
  int number;
  TCut *CrystalCut;
  TCut *CrystalCutWithoutCutG;
  TCut *PhotopeakEnergyCut;
  std::vector<TCutG*> cutg;
  TGraph *calibrationGraph;
  TH1F *doiResolution;
  TTreeFormula *Formula;
  std::vector<double> z;
};

// class of real z positions of tag points (after check of alignment)
class inputZ_t
{
public:
  int pointsFromDoi;
  int i;
  int j;
  std::vector<double> z;
  inputZ_t(int a){ pointsFromDoi = a;};
  void clear(){z.clear();};

  friend std::istream& operator>>(std::istream& input, inputZ_t& s)
  {
    input >> s.i; //read i
    input >> s.j;           //read
    for(int p = 0; p < s.pointsFromDoi; p++)
    {
      double zValue;
      input >> zValue;
      s.z.push_back(zValue);
    }
    return input;
  }

};



bool compareByNumber(const Crystal_t &a,const Crystal_t  &b)
{
  return a.number < b.number;
}

void usage()
{
  std::cout << "\t\t" << "[-i|--input] <temp.root>  [-o|--output] <output.root> [-c|--calibration] calibration.root   [-z|--zpos] zPositions.txt --pointsFromDoi N" << std::endl
            << "\t\t" << "temp.root              - complete dataset (analysis ttree) of tagging run (obtained with hadd)"   << std::endl
            << "\t\t" << "output.root            - output file name"   << std::endl
            << "\t\t" << "calibration.root       - calibration file, obtained by running ModuleCalibration on temp.root with calcDoiResWithCalibration option and providing a calibration file obtained on a top (or lat or bg) irradiation" << std::endl
            << "\t\t" << "zPositions.txt         - output file name"   << std::endl
            << "\t\t" << "N                      - points taken from doi bench - MANDATORY if --zpos is given!"   << std::endl
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
  std::string zPositionsFileName = "";
  int pointsFromDoi = 0;
  bool zPosGiven = false;

  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "calibration", required_argument, 0, 0 },
      { "zpos", required_argument, 0, 0 },
      { "pointsFromDoi", required_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:c:", longOptions, &optionIndex);
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
      zPositionsFileName = (char *)optarg;
      zPosGiven = true;
    }
    else if (c == 0 && optionIndex == 4){
      pointsFromDoi = atoi((char *)optarg);
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}

  if(zPosGiven)
  {
    if(pointsFromDoi == 0)
    {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
    }
  }
  // std::string letter[16] = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","0","P"};  //standard ordering
  // std::string number[16] = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"};  //standard ordering
  // int nmodulex = 1;
  // int nmoduley = 1;
  // int nmppcx = 4;
  // int nmppcy = 4;
  // int ncrystalsx = 2;
  // int ncrystalsy = 2;


  //read zPositions from file if given
  std::ifstream fZPos;
  fZPos.open(zPositionsFileName.c_str(),std::ios::in);
  std::vector<inputZ_t> inputZ;
  inputZ_t tempInput(pointsFromDoi);
  while(fZPos >> tempInput)
  {
    inputZ.push_back(tempInput);
    tempInput.clear();
  }

  //DEBUG
  std::cout << "InputZ -------------" << std::endl;
  for(int i = 0; i < inputZ.size(); i++)
  {
    std::cout << inputZ[i].i << " " << inputZ[i].j << " ";
    for(int j = 0 ; j < inputZ[i].z.size(); j++)
      std::cout << inputZ[i].z[j] << " ";
    std::cout << std::endl;
  }
  std::cout << "-------------" << std::endl;
  std::cout << std::endl;


  TFile *treeFile = new TFile(inputFileName.c_str());
  TTree* tree = (TTree*) treeFile->Get("adc");

  std::vector<Crystal_t> crystal;

  TFile *calibrationFile = new TFile(calibrationFileName.c_str());
  calibrationFile->cd("Module 0.0");
  TList *listModule = gDirectory->GetListOfKeys();
  int nKeysMod = listModule->GetEntries();
  std::vector<std::string> keysModName;
  // fill a vector with the leaves names
  std::string mppc_prefix("MPPC");
  for(int i = 0 ; i < nKeysMod ; i++){
    keysModName.push_back(listModule->At(i)->GetName());
  }

  std::vector<std::string> MPPCfolders;
  for(unsigned int i = 0 ; i < keysModName.size() ; i++)
  {
    if (!keysModName[i].compare(0, mppc_prefix.size(), mppc_prefix))
    {
      MPPCfolders.push_back(keysModName[i]);
    }
  }
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
    for(unsigned int iCry = 0 ; iCry < CrystalFolders.size() ; iCry++)
    {
      //  std::cout << CrystalFolders[iCry] << std::endl;
       gDirectory->cd(CrystalFolders[iCry].c_str());

       Crystal_t temp_crystal;
       temp_crystal.number = -1;
       temp_crystal.CrystalCut = NULL;
       temp_crystal.CrystalCutWithoutCutG = NULL;
       temp_crystal.PhotopeakEnergyCut = NULL;
       temp_crystal.calibrationGraph = NULL;
       temp_crystal.doiResolution = NULL;

       //get crystal number
      //  std::cout << atoi((CrystalFolders[iCry].substr(crystal_prefix.size()+1,CrystalFolders[iCry].size()-crystal_prefix.size()-1)).c_str()) << std::endl;
       temp_crystal.number = atoi((CrystalFolders[iCry].substr(crystal_prefix.size()+1,CrystalFolders[iCry].size()-crystal_prefix.size()-1)).c_str());

       if(zPosGiven) //take the z vector from input file
       {
         // run on the inputZs
         for(int iZ = 0; iZ < inputZ.size(); iZ++)
         {
           //calc crystal number. this is terrible since it's not general at all..
           int cryNum = inputZ[iZ].i * 8 + inputZ[iZ].j;
           if(cryNum == temp_crystal.number)
           {
             temp_crystal.z = inputZ[iZ].z;
           }
         }
       }
       TList *listCry = gDirectory->GetListOfKeys();
       int nKeysCry = listCry->GetEntries();
       std::vector<std::string> keysCryName;
       if(nKeysCry) //if directory not empty
       {

         for(int i = 0 ; i < nKeysCry ; i++){
           keysCryName.push_back(listCry->At(i)->GetName());
         }
         std::string CalibName;
         std::string CutName;
         std::vector<std::string> cutgNames;
         std::string cutG_prefix("cutg");
         std::string calibration_prefix("Calibration");
         std::string crystalCut_prefix("CrystalCut");
         std::string crystalCutWithoutCutG_prefix("CrystalCutWithoutCutG");
         std::string photopeakEnergyCut_prefix("PhotopeakEnergyCut");
         for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
         {
           if(!keysCryName[i].compare(0,calibration_prefix.size(),calibration_prefix)) //find calibration graph
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCanvas* C_graph = NULL;
             TGraph *calibGraph = NULL;
             C_graph = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
             if(C_graph)
               calibGraph = (TGraph*) C_graph->GetPrimitive(keysCryName[i].c_str());
             if(calibGraph)
               temp_crystal.calibrationGraph = calibGraph;
           }

           if(!keysCryName[i].compare(0,crystalCut_prefix.size(),crystalCut_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.CrystalCut = cut;
           }
           if(!keysCryName[i].compare(0,crystalCutWithoutCutG_prefix.size(),crystalCutWithoutCutG_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.CrystalCutWithoutCutG = cut;
           }
           if(!keysCryName[i].compare(0,photopeakEnergyCut_prefix.size(),photopeakEnergyCut_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.PhotopeakEnergyCut = cut;
           }

           if(!keysCryName[i].compare(0,cutG_prefix.size(),cutG_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCutG* cutg = (TCutG*) gDirectory->Get(keysCryName[i].c_str());

             temp_crystal.cutg.push_back(cutg);
           }

         }

         std::stringstream sname;
         sname << "Doi Res Cry " << temp_crystal.number;
         temp_crystal.doiResolution = new TH1F(sname.str().c_str(),sname.str().c_str(),100,-10,10);
         sname.str("");

         TCut globalCut; // the cut for the formula
         globalCut += temp_crystal.CrystalCutWithoutCutG;     // this is BasicCut (XYZ and taggingPhotopeak) + CutTrigger (TriggerChannel and broadcut)
         globalCut += temp_crystal.PhotopeakEnergyCut;        // this is the cut on photopeak energy of the corrected spectrum for this crystal
         for(unsigned int iCutg = 0; iCutg < temp_crystal.cutg.size(); iCutg++)
         {
           globalCut += temp_crystal.cutg[iCutg]->GetName();              // these are the two cutg for this crystal
         }

         sname << "Formula" << temp_crystal.number;
         TTreeFormula* Formula = new TTreeFormula(sname.str().c_str(),globalCut,tree);
         temp_crystal.Formula = Formula;

         if(temp_crystal.calibrationGraph && temp_crystal.CrystalCutWithoutCutG && temp_crystal.PhotopeakEnergyCut && (temp_crystal.cutg.size() == 2))
           crystal.push_back(temp_crystal);
       }



       gDirectory->cd("..");
    }
    calibrationFile->cd("Module 0.0");
  }

  // //BEGIN of DEBUG
  // for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  // {
  //   std::cout << crystal[i].number << "\t"
  //             << crystal[i].cut->GetTitle() << "\t"
  //             << crystal[i].calibrationGraph->GetName() << "\t";
  //   for(unsigned int j = 0 ; j < crystal[i].cutg.size(); j++)
  //   {
  //     std::cout << crystal[i].cutg[j]->GetName() << "\t";
  //   }
  //   std::cout << std::endl;
  // }
  // //END of DEBUG



  Float_t FloodX,FloodY,FloodZ;
  Float_t ZPosition;
  TBranch      *bFloodX;
  TBranch      *bFloodY;
  TBranch      *bFloodZ;
  TBranch      *bZPosition;



  tree->SetBranchAddress("FloodX", &FloodX, &bFloodX);
  tree->SetBranchAddress("FloodY", &FloodY, &bFloodY);
  tree->SetBranchAddress("FloodZ", &FloodZ, &bFloodZ);
  tree->SetBranchAddress("ZPosition", &ZPosition, &bZPosition);

  //MAIN LOOP
  long long int nevent = tree->GetEntries();
  std::cout << "Total number of events = " << nevent << std::endl;

  // for (long long int i=0;i<1000000;i++)
  for (long long int i=0;i<nevent;i++)
  {
    tree->GetEvent(i);              //read complete accepted event in memory
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {

      if(crystal[iCry].Formula->EvalInstance())  //evaluate TFormula of TCut
      {
        // std::cout << crystal[iCry].number << "\t";
        // std::cout << crystal[iCry].calibrationGraph->Eval(FloodZ) << "\t" << (crystal[iCry].calibrationGraph->Eval(FloodZ) - ZPosition) << std::endl;
        Float_t realZ = ZPosition;
        if(zPosGiven)//check ZPosition, then take the z vector from the crystal and use the closest z
        {
          float distance = INFINITY; // in mm, i doubt it will ever be bigger...
          int posMarker = -1;
          for(unsigned int iDistance = 0 ; iDistance < crystal[iCry].z.size(); iDistance++)
          {
            if( fabs(crystal[iCry].z[iDistance] - ZPosition) < distance )
            {
              posMarker = iDistance;
              distance = fabs(crystal[iCry].z[iDistance] - ZPosition);
            }
          }
          realZ = crystal[iCry].z[posMarker];
        }

        crystal[iCry].doiResolution->Fill(( crystal[iCry].calibrationGraph->Eval(FloodZ) - realZ ));
      }
      // delete Formula;
    }
  }

  std::sort(crystal.begin(), crystal.end(), compareByNumber);

  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  {
    std::stringstream sname;
    sname << "Fit Cry " << crystal[iCry].number;
    TF1 *gauss = new TF1(sname.str().c_str(),  "[0]*exp(-0.5*((x-[1])/[2])**2)",-10,10);
    gauss->SetParameter(1,crystal[iCry].doiResolution->GetMean());
    gauss->SetParameter(2,crystal[iCry].doiResolution->GetRMS());
    crystal[iCry].doiResolution->Fit(sname.str().c_str(),"Q","",crystal[iCry].doiResolution->GetMean()-3.0,crystal[iCry].doiResolution->GetMean()+2.0);
    std::cout << crystal[iCry].number << "\t" << fabs(gauss->GetParameter(2) * 2.355) << std::endl;
    crystal[iCry].doiResolution->Write();
  }

  treeFile->Close();
  calibrationFile->Close();
  outputFile->Close();
  return 0;
}
