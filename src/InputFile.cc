// InputFile.cc
// this class reads the input tchain and generates the TTree that will be used for the analysis /plot production
// Reading of config file and creation of TChain is performed directly in the class constructor, then CreateTree
// makes the TTree and FillElements is used to fill the modules, mppcs, and crystals with the information stored in
// the config file. Mainly it's the parent sons structure.

#include <iostream>
#include <sstream>

#include "TChain.h"

#include "InputFile.h"
#include <cmath>
#include "TMath.h"
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include "TError.h"
#include "TEllipse.h"
#include "TROOT.h"

// struct Point // definition fo a point (it is used in the binary output)
// {
//   float x;
//   float y;
//   float z;
// } __attribute__((__packed__));


// default constructor.
// - Reads the necessary info from the config file
// - Opens the input TChain
// - Prepares the analysis TTree
InputFile::InputFile (ConfigFile& config)
{

  gErrorIgnoreLevel = kError;
  GoodCounter = 0;
  badEvents = 0;
  counter = 0;
  //read configuration file
  fname                       = config.read<std::string>("chainName","adc");
  ncrystalsx                  = config.read<int>("ncrystalsx",2);
  ncrystalsy                  = config.read<int>("ncrystalsy");
  nmppcx                      = config.read<int>("nmppcx");
  nmppcy                      = config.read<int>("nmppcy");
  nmodulex                    = config.read<int>("nmodulex");
  nmoduley                    = config.read<int>("nmoduley");
  taggingPosition             = config.read<float>("taggingPosition",0);
  usingTaggingBench           = config.read<bool>("usingTaggingBench",0);
  taggingForTiming            = config.read<bool>("taggingForTiming",0);
  taggingCrystalChannel       = config.read<int>("taggingCrystalChannel",16);
  taggingCrystalTimingChannel = config.read<int>("taggingCrystalTimingChannel",16);
  usingRealSimData            = config.read<bool>("usingRealSimData");
  correctingSaturation        = config.read<bool>("correctingSaturation");
  saturationRun               = config.read<bool>("saturationRun",0);
  digitizer_s                 = config.read<std::string>("digitizer");
  mppc_s                      = config.read<std::string>("mppc");
  plotPositions_s             = config.read<std::string>("plotPositions");
  xPositions_s                = config.read<std::string>("xPositions");
  yPositions_s                = config.read<std::string>("yPositions");
  pedestal_s                  = config.read<std::string>("pedestal","0");
  noise_s                     = config.read<std::string>("noise","0");
  if(correctingSaturation)
  {
    saturation_s                = config.read<std::string>("saturation");
  }

  adcChannels                 = config.read<int>("digitizerTotalCh");
  nclock                      = config.read<double>("nclock",0);
  crystalx                    = config.read<double>("crystalx",1.53);
  crystaly                    = config.read<double>("crystaly",1.53);
  crystalz                    = config.read<double>("crystalz",15);
  chargeBinningADC            = config.read<double>("chargeBinningADC",156e-15);  // adc charge binning
  saturationFormat            = config.read<int>("saturationFormat",0);   // format of saturation data
  esrThickness                = config.read<double>("esrThickness",0.07);
  usingAllChannels            = config.read<bool>("usingAllChannels",0);
  wAllChannels                = config.read<bool>("wAllChannels",0);                  // whether we use the sum of all channels to compute w (true = 1) of just the neighbours (false = 0). Deafult to false.
  digitizerType               = config.read<int>("digitizerType",0);       // type of digitizer. 0 = desktop, 1 = vme, 2 = petiroc
  approximateTDCbinning       = config.read<float>("approximateTDCbinning",35.0);
  TDCcalculationEntries       = config.read<int>("TDCcalculationEntries",0);
  calculateTDCbinning         = config.read<bool>("calculateTDCbinning",0);
  minDeltaForFT               = config.read<float>("minDeltaForFT",10.0);
  pedestalTag                  = config.read<float>("pedestalTag",0.0);
  smearTaggingTime = config.read<bool>("smearTaggingTime",0);// whether to smear the time stamp of external tagging. Needed for simulations, where the tagging time stamp is always 0 (i.e. the gamma emission time) - default = 0
  sigmaTimeTag  = config.read<float>("SigmaTimeTag",30.0); // sigma for the smearing of tagging time [ps]. it's the time resolution of an hypothetical external short crystal + fast sipm - default = 30.0, which corresponds to Hamamatsu MPPC + 2x2x3 m3 LSO-Ca codoped crystal (100ps FWHM CTR, see Mythra poster)
  randGen = new TRandom3(0);

  if(digitizerType == 1 || digitizerType == 2)
  {
    timingCh_s                 = config.read<std::string>("timing","0");
    if(timingCh_s.compare("0") == 0)
    {
      std::cout << "Using same digitizer channels ordering for timing channels" << std::endl;
      timingCh_s = digitizer_s;

    }
    config.split( timingCh_f, timingCh_s, "," );
    for(unsigned int i = 0 ; i < timingCh_f.size() ; i++)
    {
      config.trim(timingCh_f[i]);
      timingCh.push_back(atoi(timingCh_f[i].c_str()));
    }
  }



  //split them using the config file class
  config.split( digitizer_f, digitizer_s, "," );
  config.split( mppc_f, mppc_s, "," );
  config.split( plotPositions_f, plotPositions_s, "," );
  config.split( xPositions_f, xPositions_s, "," );
  config.split( yPositions_f, yPositions_s, "," );
  if(correctingSaturation)
  {
    config.split( saturation_f, saturation_s, "," );
  }

  //trim them using the config file class (i.e. remove spaces)
  //and at the same time put in vectors with numbers for the ones that are numbers
  for(unsigned int i = 0 ; i < digitizer_f.size() ; i++)
  {
    config.trim(digitizer_f[i]);
    digitizer.push_back(atoi(digitizer_f[i].c_str()));
  }
  for(unsigned int i = 0 ; i < mppc_f.size() ; i++)
  {
    config.trim(mppc_f[i]);
    mppc_label.push_back(mppc_f[i]);
  }
  for(unsigned int i = 0 ; i < plotPositions_f.size() ; i++)
  {
    config.trim(plotPositions_f[i]);
    plotPositions.push_back(atoi(plotPositions_f[i].c_str()));
  }
  for(unsigned int i = 0 ; i < xPositions_f.size() ; i++)
  {
    config.trim(xPositions_f[i]);
    xPositions.push_back(atof(xPositions_f[i].c_str()));
  }
  for(unsigned int i = 0 ; i < yPositions_f.size() ; i++)
  {
    config.trim(yPositions_f[i]);
    yPositions.push_back(atof(yPositions_f[i].c_str()));
  }
  if(correctingSaturation)
  {
    for(unsigned int i = 0 ; i < saturation_f.size() ; i++)
    {
      config.trim(saturation_f[i]);
      saturation.push_back(atof(saturation_f[i].c_str()));
    }
  }

  if (pedestal_s.compare("0") == 0) // default pedestal key, means that all pedestals are set to 0 automatically
  {
    for(unsigned int i = 0; i <  digitizer.size(); i++)
    {
      pedestal.push_back(0.0);
    }
  }
  else
  {
    config.split( pedestal_f, pedestal_s, "," );
    for(unsigned int i = 0 ; i < digitizer_f.size() ; i++)
    {
      config.trim(pedestal_f[i]);
      pedestal.push_back(atof(pedestal_f[i].c_str()));
    }
  }

  if (noise_s.compare("0") == 0) // default pedestal key, means that all pedestals are set to 0 automatically
  {
    for(unsigned int i = 0; i <  digitizer.size(); i++)
    {
      noise.push_back(0.0);
    }
  }
  else
  {
    config.split( noise_f, noise_s, "," );
    for(unsigned int i = 0 ; i < digitizer_f.size() ; i++)
    {
      config.trim(noise_f[i]);
      noise.push_back(atof(noise_f[i].c_str()));
    }
  }


  //check if the vectors just built have the same size
  assert( (digitizer.size() == mppc_label.size() ) &&
          (digitizer.size() == plotPositions.size()) &&
          (digitizer.size() == xPositions.size()) &&
          (digitizer.size() == yPositions.size()) &&
          (digitizer.size() == noise.size()) &&
          (digitizer.size() == pedestal.size()));
  if(correctingSaturation)
  {
    assert(digitizer.size() == saturation.size());
  }
  if(digitizerType == 1 || digitizerType == 2)
  {
    assert(digitizer.size() == timingCh.size());
  }

  //on or off for modular analysis
  mppcOFF_s                      = config.read<std::string>("mppcOFF","");
  config.split( mppcOFF_f, mppcOFF_s, "," );
  for(unsigned int i = 0 ; i < mppcOFF_f.size() ; i++)
  {
    config.trim(mppcOFF_f[i]);
    mppcOFF.push_back(mppcOFF_f[i]);
  }


  std::vector<float> xPositionsSorted;
  std::vector<float> yPositionsSorted;
  //push the first values
  xPositionsSorted.push_back(xPositions[0]);
  yPositionsSorted.push_back(yPositions[0]);
  for(unsigned int iFill = 1; iFill < xPositions.size(); iFill++)
  {
    bool XalreadyThere = false;
    for(unsigned int iSorted = 0 ; iSorted < xPositionsSorted.size() ; iSorted++)
    {
      if(xPositions[iFill] == xPositionsSorted[iSorted])
      {
        XalreadyThere = true;
      }
    }
    if(!XalreadyThere)
    {
      xPositionsSorted.push_back(xPositions[iFill]);
    }
    bool YalreadyThere = false;
    for(unsigned int iSorted = 0 ; iSorted < yPositionsSorted.size() ; iSorted++)
    {
      if(yPositions[iFill] == yPositionsSorted[iSorted])
      {
        YalreadyThere = true;
      }
    }
    if(!YalreadyThere)
    {
      yPositionsSorted.push_back(yPositions[iFill]);
    }
  }
  //sort the position vectors
  std::sort (xPositionsSorted.begin(),xPositionsSorted.end());
  std::sort (yPositionsSorted.begin(),yPositionsSorted.end());



  for(unsigned int i = 0 ; i < digitizer.size() ; i++)
  {
    detector_t det;
    det.digitizerChannel = digitizer[i];
    if(digitizerType == 1 || digitizerType == 2)
    {
      det.timingChannel = timingCh[i];
    }
    det.label            = mppc_label[i];
    //convert saturation data to ADC channels, if necessary
    if(correctingSaturation)
    {
      if(saturationFormat == 0)
      {
        det.saturation       = saturation[i];
      }
      else
      {
        det.saturation       = saturation[i] / chargeBinningADC;
      }
    }
    // std::cout << det.saturation << std::endl;
    det.plotPosition     = plotPositions[i];
    det.xPosition        = xPositions[i];
    det.yPosition        = yPositions[i];

    //also pedestals and noise need to be corrected for saturation!
    if(correctingSaturation)
    {
      // TreeAdcChannel[iDet] = (Float_t) (-detector[iDet].saturation * TMath::Log(1.0 - ( (ADCminusPedestal)/((Float_t) detector[iDet].saturation) )));
      det.pedestal = (float) (-det.saturation * TMath::Log(1.0 - ( (pedestal[i])/((Float_t) det.saturation) )));
      det.noise    = (float) (-det.saturation * TMath::Log(1.0 - ( (noise[i])/((Float_t) det.saturation) )));
    }
    else
    {
      det.pedestal         = pedestal[i];
      det.noise            = noise[i];
    }


    det.OnForDOI         = 0;

    det.OnForModular     = true;
    for(unsigned int modCounter = 0; modCounter < mppcOFF.size(); modCounter++)
    {
      if(mppcOFF[modCounter].compare(det.label) == 0) det.OnForModular = false;
    }
    //find I and J
    int posI = -1;
    int posJ = -1;
    for(unsigned int pos = 0 ; pos < xPositionsSorted.size(); pos++)
    {
      if(det.xPosition == xPositionsSorted[pos])
      {
        posI = pos;
      }
      if(det.yPosition == yPositionsSorted[pos])
      {
        posJ = pos;
      }
    }
    det.i = posI;
    det.j = posJ;
    detector.push_back(det);
  }

  //find which channels are the neighbours
  for(unsigned int i = 0 ; i < detector.size() ; i++)
  {
     std::vector<int> neighbourChannels;
     std::vector<int> neighbourI;  //collection of i of neighbour channels
     std::vector<int> neighbourJ;  //collection of j of neighbour channels
     neighbourI.push_back(detector[i].i);
     neighbourJ.push_back(detector[i].j);

     if(detector[i].i != 0)
     {
       neighbourI.push_back(detector[i].i - 1);
     }
     if(detector[i].i != ( (int) (xPositionsSorted.size() - 1)))
     {
       neighbourI.push_back(detector[i].i + 1);
     }

     if(detector[i].j != 0)
     {
       neighbourJ.push_back(detector[i].j - 1);
     }
     if(detector[i].j != ( (int)(yPositionsSorted.size()-1)))
     {
       neighbourJ.push_back(detector[i].j + 1);
     }

     for(unsigned int iDet = 0; iDet < detector.size(); iDet++)
     {
       for(unsigned int iCheck = 0; iCheck < neighbourI.size() ; iCheck++)
       {
         if(detector[iDet].i == neighbourI[iCheck])
         {
           for(unsigned int jCheck = 0; jCheck < neighbourJ.size() ; jCheck++)
           {
             if(detector[iDet].j == neighbourJ[jCheck])
             {
               //both i and j are neighbour, then just check is not the same and add
               if(detector[iDet].digitizerChannel != detector[i].digitizerChannel)
               {
                 //both i and j are neighbour, then just check is not the same and add
                 detector[i].neighbourChannels.push_back(detector[iDet].digitizerChannel);
               }
             }
           }
         }
       }
     }
  }

  //read string for doi analysis channels
  digitizerDoi_s = config.read<std::string>("digiChannelsForDoi","8,9,10,11");
  config.split( digitizerDoi_f, digitizerDoi_s, "," );
  for(unsigned int i = 0 ; i < digitizerDoi_f.size() ; i++)
  {
    config.trim(digitizerDoi_f[i]);
    digitizerDoi.push_back(atoi(digitizerDoi_f[i].c_str()));
  }
  for(unsigned int i = 0 ; i < digitizerDoi.size() ; i++)
  {
    detector[digitizerDoi[i]].OnForDOI = 1;
  }


  crystalOFF_s                      = config.read<std::string>("crystalOFF","-1");
  config.split( crystalOFF_f, crystalOFF_s, "," );
  for(unsigned int i = 0 ; i < crystalOFF_f.size() ; i++)
  {
    config.trim(crystalOFF_f[i]);
    crystalOFF.push_back(atoi(crystalOFF_f[i].c_str()));
  }

  // global 3d plots variables for single mppcs
  global_histo3DchannelBin = config.read<int>("histo3DchannelBin",100);
  global_histo3Dmin        = config.read<float>("histo3Dmin",0);
  global_histo3Dmax        = config.read<float>("histo3Dmax",1);
  global_div               = config.read<int>("clusterLevelPrecision",10);
  global_clusterVolumeCut  = config.read<double>("clusterVolumeCut",0.001);

  //specific variables for mppcs
  specificMPPC_s          = config.read<std::string>("specificMPPCname","");
  config.split( specificMPPC_f, specificMPPC_s, "," );
  for(unsigned int i = 0 ; i < specificMPPC_f.size() ; i++)
  {
    config.trim(specificMPPC_f[i]);
    specificMPPC.push_back(specificMPPC_f[i]);
  }
  specificBin_s          = config.read<std::string>("specificBin","");
  config.split( specificBin_f, specificBin_s, "," );
  for(unsigned int i = 0 ; i < specificBin_f.size() ; i++)
  {
    config.trim(specificBin_f[i]);
    specificBin.push_back(atoi(specificBin_f[i].c_str()));
  }
  specificPrecision_s          = config.read<std::string>("specificPrecision","");
  config.split( specificPrecision_f, specificPrecision_s, "," );
  for(unsigned int i = 0 ; i < specificPrecision_f.size() ; i++)
  {
    config.trim(specificPrecision_f[i]);
    specificPrecision.push_back(atoi(specificPrecision_f[i].c_str()));
  }
  specificCut_s          = config.read<std::string>("specificCut","");
  config.split( specificCut_f, specificCut_s, "," );
  for(unsigned int i = 0 ; i < specificCut_f.size() ; i++)
  {
    config.trim(specificCut_f[i]);
    specificCut.push_back(atof(specificCut_f[i].c_str()));
  }



  crystalIsOn = new bool*[ncrystalsx*nmppcx*nmodulex];
  for(int i = 0 ; i < ncrystalsx*nmppcx*nmodulex ; i++)
  {
    crystalIsOn[i] = new bool[ncrystalsy*nmppcy*nmoduley];
    for(int j = 0 ; j < ncrystalsy*nmppcy*nmoduley ; j++)
    {
      crystalIsOn[i][j] = false;
    }
  }

  // specific MPPC strings






  //   if(digitizer.size() > 16) //FIXME is this necessary?
  //   {
  //     std::cout << "ERROR: Only one module can be analyzed at a time! Set 16 or less input channels in the config file!" << std::endl;
  //   }
  //feedback to the user
  std::cout << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << " Channels configuration " << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << "ADC input\tTime Ch\tMPPC ch\tCanvas\tx[mm]\ty[mm]\tPedestal[ADC ch]\tNoise[ADC ch]\tNeighbour channels" << std::endl;
  std::cout << "------------------------" << std::endl;
  for(unsigned int i = 0 ; i < digitizer.size() ; i++)
  {
    std::cout << "Channel[" << detector[i].digitizerChannel << "] = \t";
    if(digitizerType == 1 || digitizerType == 2)
    {
      std::cout <<  detector[i].timingChannel << "\t";
    }
    else
    {
      std::cout <<  "\t" << "\t";
    }
    std::cout <<  detector[i].label << "\t"
              << detector[i].plotPosition << "\t"
              << detector[i].xPosition << "\t"
              << detector[i].yPosition << "\t"
              << detector[i].pedestal << "\t\t\t"
              << detector[i].noise << "\t\t\t";
    for(unsigned int iNeighbour = 0; iNeighbour < detector[i].neighbourChannels.size(); iNeighbour++)
        std::cout << detector[i].neighbourChannels[iNeighbour] << " ";
    std::cout << std::endl;
  }
  std::cout << "------------------------" << std::endl;
  std::cout << std::endl;



  if(saturationRun) // if this is a saturationRun and the user has left correctingSaturation flag ON by mistake, force it to OFF (in saturationRun correcting for saturation makes no sense...)
  {
    correctingSaturation = false;
  }

  //------------------------------------------------------------------------------------------//
  //  opens the Tchain, set its branches, create the TTree that will be used for the analysis //
  //------------------------------------------------------------------------------------------//
  inputChannels          = detector.size();
  fchain                 = new TChain(fname.c_str());  // create the input tchain and the analysis ttree
  ftree                  = new TTree(fname.c_str(),fname.c_str());
  // first, create the adc channels variables and branches
  ChainAdcChannel        = new Int_t [adcChannels];
  ChainDesktopAdcChannel = new Short_t [adcChannels]; // input from ADC desktop
  ChainVMEadcChannel     = new UShort_t [adcChannels]; // input from ADC desktop
  ChainTimeStamp         = new Float_t[adcChannels];
  TDCBinning             = new Float_t[adcChannels];
  DigitizerChannelOn     = new bool[adcChannels];
  bChainAdcChannel       = new TBranch* [adcChannels];
  bChainTimeStamp        = new TBranch* [adcChannels];
  TreeAdcChannel         = new Float_t [inputChannels]; // channels analyzed
  // if(digitizerType  == 1 || digitizerType == 2)
  // {
  TreeTimeStamp          = new Float_t[inputChannels];
  // }
  for(int i = 0; i < adcChannels; i++)
  {
    DigitizerChannelOn[i] = false;
  }
  for(int i = 0; i < inputChannels; i++)
  {
    DigitizerChannelOn[detector[i].digitizerChannel] = true;
  }

}



void InputFile::ImportTChain(int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>");

  // fill the tchain with input files
  if(std::string(argv[1]) == std::string("-c")) // first argument is -c, then the config file name is passed by command line
  {
    for (int i = 3; i < argc ; i++) // run on the remaining arguments to add all the input files
    {
      std::cout << "Adding file " << argv[i] << std::endl;
      fchain->Add(argv[i]);
    }
  }
  else // the config file was indeed the default one
  {
    for (int i = 1; i < argc ; i++) // run on the remaining arguments to add all the input files
    {
      std::cout << "Adding file " << argv[i] << std::endl;
      fchain->Add(argv[i]);
    }
  }
  // set branches for reading the input files
  fchain->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag, &bChainExtendedTimeTag);
  fchain->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag, &bChainDeltaTimeTag);
  if(usingRealSimData)
  {
    fchain->SetBranchAddress("RealX", &RealX, &bRealX);
    fchain->SetBranchAddress("RealY", &RealY, &bRealY);
    fchain->SetBranchAddress("RealZ", &RealZ, &bRealZ);
    // fchain->SetBranchAddress("Tagging", &simTaggingCharge, &bsimTaggingCharge);
    // fchain->SetBranchAddress("TaggingTimeStamp", &simTaggingTime, &bsimTaggingTime);
    fchain->SetBranchAddress("CrystalsHit",&CrystalsHit, &bCrystalsHit);
    fchain->SetBranchAddress("NumbOfInteractions",&NumbOfInteractions, &bNumbOfInteractions);
    // fchain->SetBranchAddress("TotalCryEnergy",&TotalCryEnergy, &bTotalCryEnergy);
  }
  for(int i=0; i<adcChannels; i++)
  {

    if(digitizerType == 0)
    {
      std::stringstream sname;
      sname << "ch" << i;
      fchain->SetBranchAddress(sname.str().c_str(), &ChainDesktopAdcChannel[i], &bChainAdcChannel[i]);
    }
    if(digitizerType == 1)
    {
      std::stringstream sname;
      sname << "ch" << i;
      fchain->SetBranchAddress(sname.str().c_str(), &ChainVMEadcChannel[i], &bChainAdcChannel[i]);
      sname.str("");
      sname << "t" << i;
      fchain->SetBranchAddress(sname.str().c_str(), &ChainTimeStamp[i],&bChainTimeStamp[i]);
    }
  }
}

void InputFile::PrepareTTree()
{
  //set branches also for the analysis ttree
  ftree->Branch("ExtendedTimeTag",&TreeExtendedTimeTag,"ExtendedTimeTag/l");
  ftree->Branch("DeltaTimeTag",&TreeDeltaTimeTag,"DeltaTimeTag/l");
  //branches of the channels data
  for (int i = 0 ; i < inputChannels ; i++)
  {
    //empty the stringstreams
    std::stringstream sname,stype;
    sname << "ch" << detector[i].digitizerChannel;
    stype << "ch" << detector[i].digitizerChannel << "/F";
    ftree->Branch(sname.str().c_str(),&TreeAdcChannel[i],stype.str().c_str());
    sname.str("");
    stype.str("");
    if(digitizerType == 1)
    {
      sname << "t" << detector[i].digitizerChannel;
      stype << "t" << detector[i].digitizerChannel << "/F";
      ftree->Branch(sname.str().c_str(),&TreeTimeStamp[i],stype.str().c_str());
      sname.str("");
      stype.str("");
    }
    if(digitizerType == 2)
    {
      sname << "t" << detector[i].digitizerChannel;
      stype << "t" << detector[i].digitizerChannel << "/F";
      ftree->Branch(sname.str().c_str(),&TreeTimeStamp[i],stype.str().c_str());
      sname.str("");
      stype.str("");
    }
  }
  ftree->Branch("TriggerChannel",&TreeTriggerChannel,"TriggerChannel/I");
  ftree->Branch("FloodX",&TreeFloodX,"FloodX/F");
  ftree->Branch("FloodY",&TreeFloodY,"FloodY/F");
  ftree->Branch("FloodZ",&TreeFloodZ,"FloodZ/F");
  ftree->Branch("BadEvent",&TreeBadevent,"BadEvent/O");
  if(usingTaggingBench || taggingForTiming)
  {
    ftree->Branch("Tagging",&TreeTagging,"Tagging/F");
    ftree->Branch("TaggingTimeStamp",&TaggingTimeStamp,"TaggingTimeStamp/F");
    ftree->Branch("ZPosition",&TreeZPosition,"ZPosition/F");
  }
  // ftree->Branch("Theta",&TreeTheta,"Theta/F");
  // ftree->Branch("Phi",&TreePhi,"Phi/F");
  if(usingRealSimData)
  {
    // pTreeTotalCryEnergy = &TreeTotalCryEnergy;
    ftree->Branch("RealX",&TreeRealX,"RealX/F");
    ftree->Branch("RealY",&TreeRealY,"RealY/F");
    ftree->Branch("RealZ",&TreeRealZ,"RealZ/F");
    ftree->Branch("CrystalsHit",&TreeCrystalsHit,"CrystalsHit/S");
    ftree->Branch("NumbOfInteractions",&TreeNumbOfInteractions,"NumbOfInteractions/S");
    // ftree->Branch("TotalCryEnergy","std::vector<float>",&pTreeTotalCryEnergy);
  }
}


void InputFile::FillEvent()
{
  double maxCharge = 0;
  // double secondCharge = 0;
  float columnsum= 0;
  float rowsum= 0;
  float total=  0;
  // float totalForFloodZ = 0;
  TreeBadevent = false;
  TreeExtendedTimeTag = ChainExtendedTimeTag;
  TreeDeltaTimeTag = ChainDeltaTimeTag;



  // int TreeEntryCounter = 0;
  // int detectorCounter = 0;
  int TriggerID = 0;
  //loop to fill the channel
  for (int j = 0 ; j < adcChannels ; j++)
  {
    for(unsigned int iDet = 0 ; iDet < detector.size(); iDet++)
    {
      if(j == detector[iDet].digitizerChannel)
      {
        //charge part
        Float_t ADCminusPedestal = ChainAdcChannel[j] - detector[iDet].pedestal;
        // std::cout << ADCminusPedestal << std::endl;
        if(correctingSaturation)
        {
          if( ADCminusPedestal > detector[iDet].saturation)
          {
            TreeBadevent = true;
          }
          TreeAdcChannel[iDet] = (Float_t) (-detector[iDet].saturation * TMath::Log(1.0 - ( (ADCminusPedestal)/((Float_t) detector[iDet].saturation) )));
        }
        else
          TreeAdcChannel[iDet] = (Float_t) ADCminusPedestal;
            //find the max charge and therefore the TriggerChannel
        if (TreeAdcChannel[iDet] > maxCharge)
        {
          maxCharge = TreeAdcChannel[iDet];
          TreeTriggerChannel = detector[iDet].digitizerChannel;
          TriggerID = iDet;
        }
      }
      if(j == detector[iDet].timingChannel)
      {
        //timing part
        TreeTimeStamp[iDet] = (Float_t) ChainTimeStamp[j];
      }
    }
    if(usingTaggingBench || taggingForTiming)
    {
      if( j == taggingCrystalChannel)
      {
        TreeTagging = (Float_t) (ChainAdcChannel[j] - pedestalTag); // no saturation correction for the tagging crystal..
      }
      if( j == taggingCrystalTimingChannel)
      {
        TaggingTimeStamp = (Float_t) ChainTimeStamp[j];//timing part
      }
    }
  }

  //localize trigger channel detector
  //loop to calculate u,v
  // int counterFill = 0;

  for(unsigned int iDet = 0 ; iDet < detector.size(); iDet++)
  {
    bool acceptedChannel = false;
    if(usingAllChannels) //all channels for u and v, so accept all the channels
    {
      acceptedChannel = true;
    }
    else // only neighbour channels (and trigger channel itself) are accepted
    {
      if(detector[iDet].digitizerChannel == detector[TriggerID].digitizerChannel) // accept the detector if this is the trigger channel
      {
        acceptedChannel = true;
      }
      else
      {
        for(unsigned int iNeighbour = 0; iNeighbour < detector[TriggerID].neighbourChannels.size(); iNeighbour++)
        {
          if(detector[iDet].digitizerChannel == detector[TriggerID].neighbourChannels[iNeighbour]) // check if this channel is in the list of neighbours of the trigger channel
          acceptedChannel = true;
        }
      }
    }

    bool acceptedChannelW = false;
    if(wAllChannels) //all channels for w
    {
      acceptedChannelW = true;
    }
    else
    {
      if(detector[iDet].digitizerChannel == detector[TriggerID].digitizerChannel) // accept the detector if this is the trigger channel
      {
        acceptedChannelW = true;
      }
      else
      {
        for(unsigned int iNeighbour = 0; iNeighbour < detector[TriggerID].neighbourChannels.size(); iNeighbour++)
        {
          if(detector[iDet].digitizerChannel == detector[TriggerID].neighbourChannels[iNeighbour]) // check if this channel is in the list of neighbours of the trigger channel
          acceptedChannelW = true;
        }
      }
    }

    if(acceptedChannel)
    {
      rowsum    += TreeAdcChannel[iDet]*detector[iDet].xPosition;
      columnsum += TreeAdcChannel[iDet]*detector[iDet].yPosition;
    }
    if(acceptedChannelW)
    {
      total += TreeAdcChannel[iDet];
    }
  }


  //compute u,v,w
  // near channels vs. total channels depending on what decided before
  TreeFloodX = rowsum/total;
  TreeFloodY = columnsum/total;
  TreeFloodZ =  maxCharge/total;
  //     TreeFloodZ =  maxCharge/totalForFloodZ;

  if(usingTaggingBench || taggingForTiming) TreeZPosition = taggingPosition;

  if(usingRealSimData)
  {
    TreeRealX = RealX;
    TreeRealY = RealY;
    TreeRealZ = RealZ;
    TreeTagging = 1;      //FIXME for now like this...

    TaggingTimeStamp = 0; //FIXME for now like this...
    if(smearTaggingTime) // smear time tag of "tagging bench"
    {
      TaggingTimeStamp = (Float_t) ((gRandom->Gaus(0,sigmaTimeTag))*1e-12);
    }
    TreeCrystalsHit = CrystalsHit;
    TreeNumbOfInteractions = NumbOfInteractions;
  }

  if(TreeExtendedTimeTag >= nclock)
  {
    if(!TreeBadevent)
    {
      ftree->Fill();
      GoodCounter++;
    }
    else
    {
      badEvents++;
    }
  }
  //counter to give a feedback to the user
  // counter++;
}





void InputFile::FillTreeCAEN()
{
  long long int nevent = fchain->GetEntries();

  for (long long int i=0;i<nevent;i++)
  {
    fchain->GetEvent(i);
    //copy the input charges to the larger type array
    for (int j = 0 ; j < adcChannels ; j++)
    {
      if(digitizerType == 0)
      {
        ChainAdcChannel[j] = (Int_t) ChainDesktopAdcChannel[j];
      }
      else
      {
        ChainAdcChannel[j] = (Int_t) ChainVMEadcChannel[j];
      }
    }
    FillEvent();
    int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
    counter++;
  }

}





void InputFile::FillTreePetiroc(int argc, char** argv)
{
  //read the input file from Petiroc
  //open stream input
  std::ifstream inputFilePetiroc;
  std::string PetirocFileName;
  if(std::string(argv[1]) == std::string("-c")) // first argument is -c, then the config file name is passed by command line
  {
    PetirocFileName = argv[3];
  }
  else // the config file was indeed the default one
  {
    PetirocFileName = argv[1];
  }

  if(calculateTDCbinning)
  {
    std::cout << "Scanning file " << PetirocFileName << " to calculate TDC binnings..." << std::endl;
    inputFilePetiroc.open (PetirocFileName.c_str(), std::ifstream::in);
    string dummyLine;
    getline(inputFilePetiroc, dummyLine);
    inputPetirocFile_t inputPetiroc(adcChannels);
    int *maxFT;
    int *minFT;
    maxFT = new int[adcChannels];
    minFT = new int[adcChannels];
    for(int j = 0; j < adcChannels; j++)
    {
      maxFT[j] = -30000;
      minFT[j] = 30000;
    }


    // read file and update min e max
    while(inputFilePetiroc >> inputPetiroc)
    {
      if( (TDCcalculationEntries != 0) && (counter > TDCcalculationEntries)) //if TDCcalculationEntries == 0, whole input is used, otherwise up to TDCcalculationEntries events
      {
        break;
      }
      for(int j = 0 ; j < adcChannels ; j++)
      {
        // std::cout << inputPetiroc.FineTime[j]<< " minFT[" << j << "] = " <<  minFT[j] << "\t" << "maxFT[" << j << "] = " <<  maxFT[j] << std::endl;
        if(inputPetiroc.FineTime[j] > maxFT[j])
        {
          maxFT[j] = inputPetiroc.FineTime[j];
        }
        if(inputPetiroc.FineTime[j] < minFT[j])
        {
          minFT[j] = inputPetiroc.FineTime[j];
        }
      }
      inputPetiroc.clear();
      counter++;
    }
    inputFilePetiroc.close();

    //calculate the TDC binnings
    for(int j = 0 ; j < adcChannels ; j++)
    {
      // std::cout << "minFT[" << j << "] = " <<  minFT[j] << "\t" << "maxFT[" << j << "] = " <<  maxFT[j] << std::endl;
      if( (maxFT[j] - minFT[j]) > minDeltaForFT )
      {
        TDCBinning[j] = (25e-9) / (maxFT[j] - minFT[j]);
      }
      else
      {
        TDCBinning[j] = approximateTDCbinning*1e-12;
      }
      std::cout << "TDCBinning[" << j << "] = " <<  TDCBinning[j] << std::endl;
    }




    counter = 0;

  }
  else
  {
    std::cout << "Using approximate TDC binning = " << approximateTDCbinning << " ps for all channels"<< std::endl;
    for (int j = 0 ; j < adcChannels ; j++)
    {
      TDCBinning[j] = approximateTDCbinning*1e-12;
    }
  }



  std::cout << "Reading file " << PetirocFileName << std::endl;
  inputFilePetiroc.open (PetirocFileName.c_str(), std::ifstream::in);
  //skip the first line
  string dummyLine;
  getline(inputFilePetiroc, dummyLine);
  inputPetirocFile_t inputPetiroc(adcChannels);
  while(inputFilePetiroc >> inputPetiroc)
  {
    ChainExtendedTimeTag = counter; //FIXME is there an absolute time tag in this ADC?
    ChainDeltaTimeTag = 1;
    for (int j = 0 ; j < adcChannels ; j++)
    {
      ChainAdcChannel[j] = (Int_t) inputPetiroc.Charge[j];
      ChainTimeStamp[j] =  (Float_t) ( ( (inputPetiroc.CoarseTime[j] + 1)*25e-9) - (inputPetiroc.FineTime[j]*TDCBinning[j]) );
    }
    FillEvent();
    // int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
    if( (counter % 1000) == 0 )
    {
      std::cout << "\r";
      std::cout << counter << " events done... ";
      //std::cout << counter << std::endl;
    }
    inputPetiroc.clear();
    counter++;
  }
  FillEvent();
}

// Runs on the input TChain elements
// and fills the analysis TTree
void InputFile::FillTree(int argc, char** argv)
{
  //creates the TTree from the input Tchain
  std::cout << "Filling the TTree for the analysis... " << std::endl;

  if(digitizerType == 0 || digitizerType == 1)
  {
    FillTreeCAEN();
  }
  if(digitizerType == 2)
  {
    FillTreePetiroc(argc,argv);
  }
  //some feedback...
  std::cout << std::endl;
  std::cout << "Tot events = \t" << counter << std::endl;
  std::cout << "Accepted events = \t" << GoodCounter << std::endl;
  std::cout << "Bad events = \t" << badEvents << std::endl;
}


// takes the Elements created in the main file and fill them with information (hierarchy, positions, names, etc..)
void InputFile::FillElements(Module*** module,Mppc*** mppc,Crystal*** crystal)
{

  // changing i and j, no more from top left and y x, but x y from bottom left
  std::stringstream sname;
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      sname << "Module " << iModule << "." << jModule;
      module[iModule][jModule] = new Module(); // creates a default module
      module[iModule][jModule]->SetName(sname.str().c_str());          // assign a name
      sname.str("");
      sname << iModule << "." << jModule;
      module[iModule][jModule]->SetExtendedID(sname.str().c_str());
      sname.str("");
      module[iModule][jModule]->SetID(iModule*nmoduley+jModule);
      module[iModule][jModule]->SetI(iModule);
      module[iModule][jModule]->SetJ(jModule);
      module[iModule][jModule]->SetChildrenI(nmppcx);
      module[iModule][jModule]->SetChildrenJ(nmppcy);
      module[iModule][jModule]->SetSeed(randGen->GetSeed());

      std::vector<int> moduleChannels;
      for(unsigned int iDet = 0; iDet < detector.size(); iDet++)
      {
        moduleChannels.push_back(detector[iDet].digitizerChannel);
      }
      module[iModule][jModule]->SetChannels(moduleChannels);
      module[iModule][jModule]->SetDetector(detector);

      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
        for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
        {
          //calc I and J
          int mppcI = (iModule*nmppcx)+iMppc;
          int mppcJ = (jModule*nmppcy)+jMppc;

          //find corresponding detector
          int detID = 0;
          for(unsigned int iDet = 0; iDet < detector.size(); iDet++)
          {
            if(detector[iDet].i == mppcI && detector[iDet].j == mppcJ)
              detID = iDet;
          }
          // compute the position of this mppc element from the plotPositions taken from config file
          //
          //  	  int pos = (mppcI*nmppcx + mppcJ) +1 ; // +1 is because of root multicanvas starting from 1...
          // int pos = 1 + (nmppcx * nmppcy) - (mppcJ * nmppcx) - ((nmppcy-mppcI));
          // 	  std::cout << mppcI << "." << mppcJ << " " << pos << std::endl;
          // int posID;
          // for(int arrayCounter = 0 ; arrayCounter < digitizer.size() ; arrayCounter++) // get which input mppc goes here...
          // {
            // if(plotPositions[arrayCounter] == pos) posID = arrayCounter;
          // }


          sname << "MPPC " << mppcI << "." << mppcJ;

          mppc[mppcI][mppcJ] = new Mppc();
          mppc[mppcI][mppcJ]->SetName(sname.str().c_str());
          mppc[mppcI][mppcJ]->SetLabel( detector[detID].label ); // ----
          sname.str("");
          sname << iModule << "." << jModule << "-" << iMppc << "." << jMppc;

          mppc[mppcI][mppcJ]->SetIsOnForModular(detector[detID].OnForModular);
          // mppc[mppcI][mppcJ]->SetIsOnForModular(false);
          // for(int modCounter = 0; modCounter < mppcOFF.size(); modCounter++)
          // {
          //   if(mppcOFF[modCounter].compare(detector[detID].label) == 0)
          // }

          //set the global mppc variables for 3d plots, then override if specified
          mppc[mppcI][mppcJ]->SetHisto3DchannelBin(global_histo3DchannelBin);
          mppc[mppcI][mppcJ]->SetHisto3Dmin(global_histo3Dmin);
          mppc[mppcI][mppcJ]->SetHisto3Dmax(global_histo3Dmax);
          mppc[mppcI][mppcJ]->SetClusterLevelPrecision(global_div);
          mppc[mppcI][mppcJ]->SetClusterVolumeCut(global_clusterVolumeCut);
          //now override
          for(unsigned int modCounter = 0; modCounter < specificMPPC.size(); modCounter++)
          {
            if(specificMPPC[modCounter].compare(detector[detID].label) == 0)
            {
              mppc[mppcI][mppcJ]->SetHisto3DchannelBin(specificBin[modCounter]);
              mppc[mppcI][mppcJ]->SetClusterLevelPrecision(specificPrecision[modCounter]);
              mppc[mppcI][mppcJ]->SetClusterVolumeCut(specificCut[modCounter]);
            }
          }


          mppc[mppcI][mppcJ]->SetExtendedID(sname.str().c_str());
          mppc[mppcI][mppcJ]->SetID(detID);                  // ----
          mppc[mppcI][mppcJ]->SetI(mppcI);
          mppc[mppcI][mppcJ]->SetJ(mppcJ);
          mppc[mppcI][mppcJ]->SetChildrenI(ncrystalsx);
          mppc[mppcI][mppcJ]->SetChildrenJ(ncrystalsy);
          mppc[mppcI][mppcJ]->SetPosition(detector[detID].xPosition,detector[detID].yPosition,0);
          mppc[mppcI][mppcJ]->SetDigitizerChannel(detector[detID].digitizerChannel);
          mppc[mppcI][mppcJ]->SetCanvasPosition(detector[detID].plotPosition);
          mppc[mppcI][mppcJ]->SetParentName(module[iModule][jModule]->GetName());

          mppc[mppcI][mppcJ]->SetIsOnForDoi(detector[detID].OnForDOI);
          mppc[mppcI][mppcJ]->SetNeighbours(detector[detID].neighbourChannels);
          // for(int iDoi = 0 ; iDoi < digitizerDoi.size(); iDoi++)
          // {
          //   if(mppc[mppcI][mppcJ]->GetDigitizerChannel() == digitizerDoi[iDoi])
          //     mppc[mppcI][mppcJ]->SetIsOnForDoi(true);
          // }

          // 	  mppc[mppcI][mppcJ]->Print();
          sname.str("");

          module[iModule][jModule]->AddChild( mppc[mppcI][mppcJ]->GetName() );

          for(int iCry = 0; iCry < ncrystalsx ; iCry++)
          {
            for(int jCry = 0; jCry < ncrystalsy ; jCry++)
            {
              int cryI = (iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry);
              int cryJ = (jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry);
              sname << "Crystal " << cryI << "." << cryJ;
              crystal[cryI][cryJ] = new Crystal();
              crystal[cryI][cryJ]->SetName(sname.str().c_str());   // assign a name
              sname.str("");
              sname << iModule << "." << jModule << "-" << iMppc << "." << jMppc << "-" << iCry << "." << jCry;
              crystal[cryI][cryJ]->SetExtendedID(sname.str().c_str());
              sname.str("");
              crystal[cryI][cryJ]->SetID(cryI*ncrystalsx*nmppcx + cryJ);          // assign an ID  //----

              crystal[cryI][cryJ]->SetIsOnForModular(true);
              for(unsigned int modCounter = 0; modCounter < crystalOFF.size(); modCounter++)
              {
                if(crystalOFF[modCounter] == (cryI*ncrystalsx*nmppcx + cryJ) ) crystal[cryI][cryJ]->SetIsOnForModular(false);
              }

              crystal[cryI][cryJ]->SetI(cryI);
              crystal[cryI][cryJ]->SetJ(cryJ);
              crystal[cryI][cryJ]->SetParentName(mppc[mppcI][mppcJ]->GetName());
              crystal[cryI][cryJ]->SetCrystalOn(crystalIsOn[cryI][cryJ]);
              crystal[cryI][cryJ]->SetPosition(
                detector[detID].xPosition + iCry*(crystalx+esrThickness) - (ncrystalsx-1)*((crystalx+esrThickness)/2.0)
                , detector[detID].yPosition + jCry*(crystaly+esrThickness) - (ncrystalsy-1)*((crystaly+esrThickness)/2.0)
                ,  0  ); //FIXME z not useful so set to 0, but maybe we should put the real one..
              crystal[cryI][cryJ]->SetDimension(crystalx,crystaly,crystalz);
              mppc[mppcI][mppcJ]->AddChild( crystal[cryI][cryJ]->GetName() );
            }
          }
        }
      }
    }
  }

  // summarize the situation to the user
  std::cout << std::endl;
  std::cout << "MPPC ID mapping ";
  std::cout << std::endl;
  for(int jMppc = nmppcy-1; jMppc >=0 ; jMppc--)
  {
    for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
    {
      std::cout <<  mppc[iMppc][jMppc]->GetLabel() << "\t";

      //       std::cout <<  mppc[iMppc][jMppc]->GetLabel() << " " <<  mppc[iMppc][jMppc]->GetI() << "." << mppc[iMppc][jMppc]->GetJ() << " " << mppc[iMppc][jMppc]->GetDigitizerChannel();  /*<< "\t"*/;
      //       std::cout <<  " " <<mppc[iMppc][jMppc]->GetX() << " " <<  mppc[iMppc][jMppc]->GetY();
      //       std::cout << std::endl;
    }  std::cout << std::endl;

  }
  std::cout << std::endl;
  std::cout << "Digitizer Channels mapping ";
  std::cout << std::endl;
  for(int jMppc = nmppcy-1; jMppc >=0 ; jMppc--)
  {
    for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
    {
      std::cout << mppc[iMppc][jMppc]->GetDigitizerChannel() << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Crystal ID mapping ";
  std::cout << std::endl;
  for(int jCrystal = ncrystalsy*nmppcy*nmoduley -1 ; jCrystal >=0  ; jCrystal--)
  {
    for(int iCrystal = 0 ; iCrystal < ncrystalsx*nmppcx*nmodulex   ; iCrystal++)
    {
      std::cout << crystal[iCrystal][jCrystal]->GetID() << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Crystal Positions mapping ";
  std::cout << std::endl;
  for(int jCrystal = ncrystalsy*nmppcy*nmoduley -1 ; jCrystal >=0  ; jCrystal--)
  {
    for(int iCrystal = 0 ; iCrystal < ncrystalsx*nmppcx*nmodulex   ; iCrystal++)
    {
      std::cout << "(" << crystal[iCrystal][jCrystal]->GetX() << "," << crystal[iCrystal][jCrystal]->GetY() << ")"  << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;



}
