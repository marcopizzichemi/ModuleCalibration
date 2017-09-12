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

struct Point // definition fo a point (it is used in the binary output)
{
  float x;
  float y;
  float z;
} __attribute__((__packed__));


// default constructor.
// - Reads the necessary info from the config file
// - Opens the input TChain
// - Prepares the analysis TTree
InputFile::InputFile (int argc, char** argv, ConfigFile& config)
{

  gErrorIgnoreLevel = kError;
  //read configuration file
  fname                       = config.read<std::string>("chainName");
  ncrystalsx                  = config.read<int>("ncrystalsx",2);
  ncrystalsy                  = config.read<int>("ncrystalsy");
  nmppcx                      = config.read<int>("nmppcx");
  nmppcy                      = config.read<int>("nmppcy");
  nmodulex                    = config.read<int>("nmodulex");
  nmoduley                    = config.read<int>("nmoduley");
  taggingPosition             = config.read<float>("taggingPosition",0);
  usingTaggingBench           = config.read<bool>("usingTaggingBench",0);
  taggingCrystalChannel       = config.read<int>("taggingCrystalChannel",16);
  usingRealSimData            = config.read<bool>("usingRealSimData");
  binary                      = config.read<bool>("binary");
  correctingSaturation        = config.read<bool>("correctingSaturation");
  saturationRun               = config.read<bool>("saturationRun",0);
  BinaryOutputFileName        = config.read<std::string>("output");
  BinaryOutputFileName       += ".bin";
  //read the strings that describe the input channels
  digitizer_s                 = config.read<std::string>("digitizer");
  mppc_s                      = config.read<std::string>("mppc");
  plotPositions_s             = config.read<std::string>("plotPositions");
  xPositions_s                = config.read<std::string>("xPositions");
  yPositions_s                = config.read<std::string>("yPositions");
  saturation_s                = config.read<std::string>("saturation");
  adcChannels                 = config.read<int>("digitizerTotalCh");
  nclock                      = config.read<double>("nclock",0);
  crystalx                    = config.read<double>("crystalx",1.53);
  crystaly                    = config.read<double>("crystaly",1.53);
  crystalz                    = config.read<double>("crystalz",15);
  chargeBinningADC            = config.read<double>("chargeBinningADC",156e-15);  // adc charge binning
  saturationFormat            = config.read<int>("saturationFormat",1);   // format of saturation data
  esrThickness                = config.read<double>("esrThickness",0.07);
  usingAllChannels            = config.read<bool>("usingAllChannels",1);
  wAllChannels                = config.read<bool>("wAllChannels",0);                  // whether we use the sum of all channels to compute w (true = 1) of just the neighbours (false = 0). Deafult to false.

  //split them using the config file class
  config.split( digitizer_f, digitizer_s, "," );
  config.split( mppc_f, mppc_s, "," );
  config.split( plotPositions_f, plotPositions_s, "," );
  config.split( xPositions_f, xPositions_s, "," );
  config.split( yPositions_f, yPositions_s, "," );
  config.split( saturation_f, saturation_s, "," );
  //trim them using the config file class (i.e. remove spaces)
  //and at the same time put in vectors with numbers for the ones that are numbers
  for(int i = 0 ; i < digitizer_f.size() ; i++)
  {
    config.trim(digitizer_f[i]);
    digitizer.push_back(atoi(digitizer_f[i].c_str()));
  }
  for(int i = 0 ; i < mppc_f.size() ; i++)
  {
    config.trim(mppc_f[i]);
    mppc_label.push_back(mppc_f[i]);
  }
  for(int i = 0 ; i < plotPositions_f.size() ; i++)
  {
    config.trim(plotPositions_f[i]);
    plotPositions.push_back(atoi(plotPositions_f[i].c_str()));
  }
  for(int i = 0 ; i < xPositions_f.size() ; i++)
  {
    config.trim(xPositions_f[i]);
    xPositions.push_back(atof(xPositions_f[i].c_str()));
  }
  for(int i = 0 ; i < yPositions_f.size() ; i++)
  {
    config.trim(yPositions_f[i]);
    yPositions.push_back(atof(yPositions_f[i].c_str()));
  }
  for(int i = 0 ; i < saturation_f.size() ; i++)
  {
    config.trim(saturation_f[i]);
    saturation.push_back(atof(saturation_f[i].c_str()));
  }
  //check if the vectors just built have the same size
  assert( (digitizer.size() == mppc_label.size() ) && (digitizer.size() == plotPositions.size()) && (digitizer.size() == xPositions.size()) && (digitizer.size() == yPositions.size()) && (digitizer.size() == saturation.size()) );




  //on or off for modular analysis
  mppcOFF_s                      = config.read<std::string>("mppcOFF","");
  config.split( mppcOFF_f, mppcOFF_s, "," );
  for(int i = 0 ; i < mppcOFF_f.size() ; i++)
  {
    config.trim(mppcOFF_f[i]);
    mppcOFF.push_back(mppcOFF_f[i]);
  }


  std::vector<float> xPositionsSorted;
  std::vector<float> yPositionsSorted;
  //push the first values
  xPositionsSorted.push_back(xPositions[0]);
  yPositionsSorted.push_back(yPositions[0]);
  for(int iFill = 1; iFill < xPositions.size(); iFill++)
  {
    bool XalreadyThere = false;
    for(int iSorted = 0 ; iSorted < xPositionsSorted.size() ; iSorted++)
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
    for(int iSorted = 0 ; iSorted < yPositionsSorted.size() ; iSorted++)
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



  for(int i = 0 ; i < digitizer.size() ; i++)
  {
    detector_t det;
    det.digitizerChannel = digitizer[i];
    det.label            = mppc_label[i];
    //convert saturation data to ADC channels, if necessary
    if(saturationFormat == 0)
    {
      det.saturation       = saturation[i];
    }
    else
    {
      det.saturation       = saturation[i] / chargeBinningADC;
    }
    // std::cout << det.saturation << std::endl;
    det.plotPosition     = plotPositions[i];
    det.xPosition        = xPositions[i];
    det.yPosition        = yPositions[i];
    det.OnForDOI         = 0;

    det.OnForModular     = true;
    for(int modCounter = 0; modCounter < mppcOFF.size(); modCounter++)
    {
      if(mppcOFF[modCounter].compare(det.label) == 0) det.OnForModular = false;
    }
    //find I and J
    int posI = -1;
    int posJ = -1;
    for(int pos = 0 ; pos < xPositionsSorted.size(); pos++)
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
  for(int i = 0 ; i < detector.size() ; i++)
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
     if(detector[i].i != (xPositionsSorted.size() - 1))
     {
       neighbourI.push_back(detector[i].i + 1);
     }

     if(detector[i].j != 0)
     {
       neighbourJ.push_back(detector[i].j - 1);
     }
     if(detector[i].j != (yPositionsSorted.size()-1))
     {
       neighbourJ.push_back(detector[i].j + 1);
     }

     for(int iDet = 0; iDet < detector.size(); iDet++)
     {
       for(int iCheck = 0; iCheck < neighbourI.size() ; iCheck++)
       {
         if(detector[iDet].i == neighbourI[iCheck])
         {
           for(int jCheck = 0; jCheck < neighbourJ.size() ; jCheck++)
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
  for(int i = 0 ; i < digitizerDoi_f.size() ; i++)
  {
    config.trim(digitizerDoi_f[i]);
    digitizerDoi.push_back(atoi(digitizerDoi_f[i].c_str()));
  }
  for(int i = 0 ; i < digitizerDoi.size() ; i++)
  {
    detector[digitizerDoi[i]].OnForDOI = 1;
  }


  crystalOFF_s                      = config.read<std::string>("crystalOFF","-1");
  config.split( crystalOFF_f, crystalOFF_s, "," );
  for(int i = 0 ; i < crystalOFF_f.size() ; i++)
  {
    config.trim(crystalOFF_f[i]);
    crystalOFF.push_back(atoi(crystalOFF_f[i].c_str()));
  }

  // global 3d plots variables for single mppcs
  global_histo3DchannelBin = config.read<int>("histo3DchannelBin",100);
  global_div               = config.read<int>("clusterLevelPrecision",10);
  global_clusterVolumeCut  = config.read<double>("clusterVolumeCut",0.001);

  //specific variables for mppcs
  specificMPPC_s          = config.read<std::string>("specificMPPCname","");
  config.split( specificMPPC_f, specificMPPC_s, "," );
  for(int i = 0 ; i < specificMPPC_f.size() ; i++)
  {
    config.trim(specificMPPC_f[i]);
    specificMPPC.push_back(specificMPPC_f[i]);
  }
  specificBin_s          = config.read<std::string>("specificBin","");
  config.split( specificBin_f, specificBin_s, "," );
  for(int i = 0 ; i < specificBin_f.size() ; i++)
  {
    config.trim(specificBin_f[i]);
    specificBin.push_back(atoi(specificBin_f[i].c_str()));
  }
  specificPrecision_s          = config.read<std::string>("specificPrecision","");
  config.split( specificPrecision_f, specificPrecision_s, "," );
  for(int i = 0 ; i < specificPrecision_f.size() ; i++)
  {
    config.trim(specificPrecision_f[i]);
    specificPrecision.push_back(atoi(specificPrecision_f[i].c_str()));
  }
  specificCut_s          = config.read<std::string>("specificCut","");
  config.split( specificCut_f, specificCut_s, "," );
  for(int i = 0 ; i < specificCut_f.size() ; i++)
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
  std::cout << "ADC input\tMPPC ch\tCanvas\tx[mm]\ty[mm]\tNeighbour channels" << std::endl;
  std::cout << "------------------------" << std::endl;
  for(int i = 0 ; i < digitizer.size() ; i++)
  {
    std::cout << "Channel[" << detector[i].digitizerChannel << "] = \t" <<  detector[i].label << "\t" << detector[i].plotPosition << "\t" << detector[i].xPosition << "\t" << detector[i].yPosition << "\t";
    for(int iNeighbour = 0; iNeighbour < detector[i].neighbourChannels.size(); iNeighbour++)
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
  inputChannels       = detector.size();
  fchain              = new TChain(fname.c_str());  // create the input tchain and the analysis ttree
  ftree               = new TTree(fname.c_str(),fname.c_str());
  // first, create the adc channels variables and branches
  ChainAdcChannel     = new Short_t [adcChannels]; // input from ADC
  DigitizerChannelOn  = new bool[adcChannels];
  bChainAdcChannel    = new TBranch* [adcChannels];
  TreeAdcChannel      = new Short_t [inputChannels]; // channels analyzed
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
    fchain->SetBranchAddress("CrystalsHit",&CrystalsHit, &bCrystalsHit);
    fchain->SetBranchAddress("NumbOfInteractions",&NumbOfInteractions, &bNumbOfInteractions);
    // fchain->SetBranchAddress("TotalCryEnergy",&TotalCryEnergy, &bTotalCryEnergy);

  }
  for(int i=0; i<adcChannels; i++)
  {
    std::stringstream sname;
    sname << "ch" << i;
    fchain->SetBranchAddress(sname.str().c_str(), &ChainAdcChannel[i], &bChainAdcChannel[i]);
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
    stype << "ch" << detector[i].digitizerChannel << "/S";
    ftree->Branch(sname.str().c_str(),&TreeAdcChannel[i],stype.str().c_str());
  }
  if(usingTaggingBench) ftree->Branch("Tagging",&TreeTagging,"Tagging/S");
  ftree->Branch("TriggerChannel",&TreeTriggerChannel,"TriggerChannel/I");
  ftree->Branch("FloodX",&TreeFloodX,"FloodX/F");
  ftree->Branch("FloodY",&TreeFloodY,"FloodY/F");
  ftree->Branch("FloodZ",&TreeFloodZ,"FloodZ/F");
  if(usingTaggingBench) ftree->Branch("ZPosition",&TreeZPosition,"TreeZPosition/F");
  ftree->Branch("Theta",&TreeTheta,"Theta/F");
  ftree->Branch("Phi",&TreePhi,"Phi/F");
  ftree->Branch("BadEvent",&TreeBadevent,"BadEvent/O");
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

// Runs on the input TChain elements
// and fills the analysis TTree
void InputFile::FillTree()
{
  //creates the TTree from the input Tchain
  std::cout << "Filling the TTree for the analysis... " << std::endl;
  Int_t nevent = fchain->GetEntries();
  long long int GoodCounter = 0;
  long long int badEvents = 0;
  long long int counter = 0;

  Point point;
  ofstream output_file;
  if(binary)
    output_file.open(BinaryOutputFileName.c_str(), std::ios::binary);

  for(int i = 0; i < adcChannels; i++)
  {
    DigitizerChannelOn[i] = false;
  }
  for(int i = 0; i < inputChannels; i++)
  {
    DigitizerChannelOn[detector[i].digitizerChannel] = true;
  }


  for (Int_t i=0;i<nevent;i++)
    //   for(Int_t i=0;i<2;i++)
  {
    //loop on all the entries of tchain
    fchain->GetEvent(i);              //read complete accepted event in memory
    double maxCharge = 0;
    double secondCharge = 0;
    float columnsum= 0;
    float rowsum= 0;
    float total=  0;
    float totalForFloodZ = 0;
    TreeBadevent = false;
    TreeExtendedTimeTag = ChainExtendedTimeTag;
    TreeDeltaTimeTag = ChainDeltaTimeTag;

    // int TreeEntryCounter = 0;
    int detectorCounter = 0;
    int TriggerID = 0;
    //loop to fill the channel
    for (int j = 0 ; j < adcChannels ; j++)
    {
      for(int iDet = 0 ; iDet < detector.size(); iDet++)
      {
        if(j == detector[iDet].digitizerChannel)
        {
          if(correctingSaturation)
          {
            if( ChainAdcChannel[j] > detector[iDet].saturation)
            {
              TreeBadevent = true;
            }
            TreeAdcChannel[iDet] = (Short_t)round(-detector[iDet].saturation * TMath::Log(1.0 - ( ChainAdcChannel[j]/detector[iDet].saturation )));
          }
          else
            TreeAdcChannel[iDet] = ChainAdcChannel[j];
              //find the max charge and therefore the TriggerChannel
          if (TreeAdcChannel[iDet] > maxCharge)
          {
            maxCharge = TreeAdcChannel[iDet];
            TreeTriggerChannel = detector[iDet].digitizerChannel;
            TriggerID = iDet;
          }
          if(usingTaggingBench)
          {
            if( j == taggingCrystalChannel)
            {
              TreeTagging = ChainAdcChannel[j]; // no saturation correction for the tagging crystal..
                //this is the tagging crystal data
            }
          }
        }
      }
    }

    //localize trigger channel detector


    //loop to calculate u,v
    // int counterFill = 0;

    for(int iDet = 0 ; iDet < detector.size(); iDet++)
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
          for(int iNeighbour = 0; iNeighbour < detector[TriggerID].neighbourChannels.size(); iNeighbour++)
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
          for(int iNeighbour = 0; iNeighbour < detector[TriggerID].neighbourChannels.size(); iNeighbour++)
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

    if(usingTaggingBench) TreeZPosition = taggingPosition;

    TreeTheta = std::acos(TreeFloodZ /( std::sqrt( std::pow(TreeFloodX - 2.8,2) + std::pow(TreeFloodY - (-1.0),2) + std::pow(TreeFloodZ,2)) ));
    TreePhi =  std::atan (TreeFloodY / TreeFloodX);

    point.x = TreeFloodX;
    point.y = TreeFloodY;
    point.z = TreeFloodZ;

    if(usingRealSimData)
    {
      TreeRealX = RealX;
      TreeRealY = RealY;
      TreeRealZ = RealZ;
      TreeCrystalsHit = CrystalsHit;
      TreeNumbOfInteractions = NumbOfInteractions;
      // for(int k=0; k<(TotalCryEnergy->size()); k++)
        //std::cout<< TotalCryEnergy->at(k) <<std::endl;
        // TreeTotalCryEnergy.push_back(TotalCryEnergy->at(k));

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
    if(binary)
      output_file.write((char*)&point,sizeof(point));

    //counter to give a feedback to the user
    counter++;
    // TreeTotalCryEnergy.clear();

    int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
  }

  //some feedback...
  std::cout << std::endl;
  std::cout << "Tot events = \t" << counter << std::endl;
  std::cout << "Accepted events = \t" << GoodCounter << std::endl;
  std::cout << "Bad events = \t" << badEvents << std::endl;

  //close the binary file if it was opened
  if(binary)
    output_file.close();
  //   std::cout << "Accepted events = \t" << GoodCounter << std::endl;
  //std::cout << "Bad events = \t" << badEvents << std::endl;

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

      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
        for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
        {
          //calc I and J
          int mppcI = (iModule*nmppcx)+iMppc;
          int mppcJ = (jModule*nmppcy)+jMppc;

          //find corresponding detector
          int detID;
          for(int iDet = 0; iDet < detector.size(); iDet++)
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
          mppc[mppcI][mppcJ]->SetClusterLevelPrecision(global_div);
          mppc[mppcI][mppcJ]->SetClusterVolumeCut(global_clusterVolumeCut);
          //now override
          for(int modCounter = 0; modCounter < specificMPPC.size(); modCounter++)
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
              for(int modCounter = 0; modCounter < crystalOFF.size(); modCounter++)
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
