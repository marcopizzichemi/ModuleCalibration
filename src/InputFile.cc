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
InputFile::InputFile (ConfigFile& config)
{
  gErrorIgnoreLevel = kError;     // set ROOT error level
  ReadConfig(config);             //read config file and fill the variables

  //prepare and fill detector struct - all info except isNeighbour and isCross (because those will change from trigger to trigger and will be used temporarily)
  for(unsigned int i = 0 ; i < digitizer.size() ; i++)
  {
    detector_t det;    // create a detector_t structure
    det.digitizerChannel = digitizer[i];    // set channel
    det.label            = mppc_label[i];   // set label
    if(correctingSaturation)                // if saturation parameters are given
    {
      if(electronics == 0)                  // if electronics is CAEN x740
      {
        //then saturation is only one parameter
        det.saturation.push_back(saturation[i]);
      }
      else if(electronics == 1)             // if electronics is NINO
      {
        //then saturation is ninoSaturationParameters parameters per channel
        for(int pSat = 0 ; pSat < ninoSaturationParameters ; pSat++)
          det.saturation.push_back(saturation[i*ninoSaturationParameters + pSat]);
      }
    }
    det.plotPosition     = plotPositions[i];//canvas position
    det.xPosition        = xPositions[i];   //mppc x
    det.yPosition        = yPositions[i];   //mppc y
    detector.push_back(det);  // add this struct to a std::vector of detector_t
  }

  inputChannels       = detector.size();
  //set channels on, i.e. channel biased
  DigitizerChannelOn  = new bool[adcChannels];
  for(int i = 0; i < adcChannels; i++) // first set all the possible inputs from ADC to OFF
  {
    DigitizerChannelOn[i] = false;
  }
  for(int i = 0; i < inputChannels; i++) // then turn on only the ones were there is bias (i.e. the ones mentioned in the config digitizer key)
  {
    DigitizerChannelOn[detector[i].digitizerChannel] = true;
  }
  std::vector<int> all;
  for(unsigned int i = 0 ; i < detector.size() ; i++)
  {
    all.push_back(detector[i].digitizerChannel);
  }
  //find neighbour channels for each detector, write the numbers in the detector structs
  for(unsigned int i = 0 ; i < detector.size() ; i++)
  {
    detector[i].neighbours = FindNeighbours(detector[i].digitizerChannel);
    detector[i].cross      = FindCross(detector[i].digitizerChannel);
    detector[i].all        = all;

    //now assign the vectors
    //u-v
    // std::cout << allChannels << " " << neighbourChannels << " " << crossChannels << std::endl;
    // std::cout << relevantForUV << std::endl;
    switch(relevantForUV)
    {
      case allChannels:
        // std::cout << "allChannels" << std::endl;
        detector[i].relevantForUV = detector[i].all;
        break;
      case neighbourChannels:
        // std::cout << "neighbourChannels" << std::endl;
        detector[i].relevantForUV = detector[i].neighbours;
        break;
      case crossChannels:
        // std::cout << "crossChannels" << std::endl;
        detector[i].relevantForUV = detector[i].cross;
        break;
    }
    // std::cout << relevantForW<< std::endl;
    switch(relevantForW)
    {
      case allChannels:
        detector[i].relevantForW = detector[i].all;
        break;
      case neighbourChannels:
        detector[i].relevantForW = detector[i].neighbours;
        break;
      case crossChannels:
        detector[i].relevantForW = detector[i].cross;
        break;
    }
    // std::cout << relevantForE << std::endl;
    switch(relevantForE)
    {
      case allChannels:
        detector[i].relevantForE = detector[i].all;
        break;
      case neighbourChannels:
        detector[i].relevantForE = detector[i].neighbours;
        break;
      case crossChannels:
        detector[i].relevantForE = detector[i].cross;
        break;
    }
  }

  //crystals
  //set all crystals to off, to start
  crystalIsOn = new bool*[ncrystalsx*nmppcx*nmodulex];
  for(int i = 0 ; i < ncrystalsx*nmppcx*nmodulex ; i++)
  {
    crystalIsOn[i] = new bool[ncrystalsy*nmppcy*nmoduley];
    for(int j = 0 ; j < ncrystalsy*nmppcy*nmoduley ; j++)
    {
      crystalIsOn[i][j] = false;
    }
  }


  //feedback to the user
  std::cout << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << " Channels configuration " << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << "ADC input\tMPPC ch\tCanvas\tx[mm]\ty[mm]\tSaturation Parameters" << std::endl;
  std::cout << "------------------------" << std::endl;
  for(unsigned int i = 0 ; i < digitizer.size() ; i++)
  {
    std::cout   << "Channel["
                << detector[i].digitizerChannel << "] = \t"
                << detector[i].label << "\t"
                << detector[i].plotPosition << "\t"
                << detector[i].xPosition << "\t"
                << detector[i].yPosition << "\t";
    for(unsigned int j = 0 ; j < detector[0].saturation.size() ; j++)
      std::cout << detector[i].saturation[j] << "\t";
    std::cout << std::endl;
  }
  std::cout << "------------------------" << std::endl;
  std::cout << std::endl;






  //------------------------------------------------------------------------------------------//
  //  opens the Tchain, set its branches, create the TTree that will be used for the analysis //
  //------------------------------------------------------------------------------------------//

  fchain                  = new TChain(fname.c_str());  // create the input tchain and the analysis ttree
  ftree                   = new TTree(fname.c_str(),fname.c_str());
  // first, create the adc channels variables and branches
  CAENx740AdcChannel      = new Short_t [adcChannels]; // input from ADC
  bCAENx740AdcChannel     = new TBranch* [adcChannels];
  // analysis variables
  TreeCAENx740AdcChannel  = new Short_t [inputChannels]; // channels analyzed
  TreeNINOChargeChannel   = new Float_t [inputChannels]; // channels analyzed
  TreeNINOTimeChannel     = new Float_t [inputChannels]; // channels analyzed
}


//Imports TChain from CAEN x740 output or from g4matrix simulation
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
    fchain->SetBranchAddress("TotalEnergyDeposited",&TotalEnergyDeposited, &bTotalEnergyDeposited);
  }
  for(int i=0; i<adcChannels; i++)
  {
    std::stringstream sname;
    sname << "ch" << i;
    fchain->SetBranchAddress(sname.str().c_str(), &CAENx740AdcChannel[i], &bCAENx740AdcChannel[i]);
  }
}


void InputFile::FillTreeNINO(int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>");
  int startFiles = 1;

  // first argument is -c, then the config file name is passed by command line
  if(std::string(argv[1]) == std::string("-c")) startFiles = 3;

  for (int m = startFiles; m < argc ; m++) // run on the remaining arguments to add all the input files
  {
    std::ifstream inFile;
    inFile.open( argv[m],std::ios::in);
    std::cout << "Adding file " << argv[m] << std::endl;

    while(!inFile.eof())
    // for(int aaa = 0 ; aaa < 2 ; aaa++) //TEMP for debugging
    {
      //take all data
      Float_t tempT[adcChannels];
      Float_t tempW[adcChannels];
      int TreeEntryCounter = 0;
      double maxCharge = 0;
      float columnsum= 0;
      float rowsum= 0;
      float totalUV=  0;
      float totalW=  0;

      for (int i = 0; i < adcChannels; i++)
      {
        inFile >> tempT[i];
      }
      for (int i = 0; i < adcChannels; i++)
      {
        inFile >> tempW[i];
      }

      for (int j = 0 ; j < adcChannels ; j++) // runs on all the possible inputs of ADC/electronics, but then uses only the biased ones
      {
        if(DigitizerChannelOn[j])
        {
          // fill tree with data from the channels
          // also correcting for saturation if it's set in the config file

          TreeNINOTimeChannel[TreeEntryCounter] = (Float_t) tempT[j];
          detector[TreeEntryCounter].EventTime = TreeNINOTimeChannel[TreeEntryCounter];
          if(correctingSaturation)
          {
            TreeNINOChargeChannel[TreeEntryCounter] = (Float_t) (-1./detector[TreeEntryCounter].saturation[2])*(TMath::Log(1-((1./detector[TreeEntryCounter].saturation[1])*(TMath::Exp((tempW[j])/detector[TreeEntryCounter].saturation[0])-detector[TreeEntryCounter].saturation[3]))));

          }
          else
            TreeNINOChargeChannel[TreeEntryCounter] = (Float_t) tempW[j];

          detector[TreeEntryCounter].EventCharge = TreeNINOChargeChannel[TreeEntryCounter];
          //find the max charge and therefore the TriggerChannel
          if (TreeNINOChargeChannel[TreeEntryCounter] > maxCharge)
          {
            maxCharge = TreeNINOChargeChannel[TreeEntryCounter];
            TreeTriggerChannel = digitizer[TreeEntryCounter];
          }

          TreeEntryCounter++;
        }

        if(usingTaggingBench)
        {
          if( j == taggingCrystalChannel)
          {
            TreeNINOtaggingCharge = (Float_t) tempW[j]; // no saturation correction for the tagging crystal..
            TreeNINOtaggingTime   = (Float_t) tempT[j];
          }
        }
      }

      // more general implementation. We already have the all, neighbours and cross digitizer channels for each detector. We also know the TreeTriggerChannel. Let's get that channel and the relevant vectors
      for(unsigned int i = 0 ; i < detector.size() ; i++) // run on all detectors
      {
        for(unsigned int a = 0; a < detector[TreeTriggerChannel].relevantForUV.size() ; a++) // for as many detectors as the ones set to relevantForUV
        {
          if(detector[i].digitizerChannel == detector[TreeTriggerChannel].relevantForUV[a]) // if the detector is among the ones that are relevantForUV, use them for U-V calculation
          {
            rowsum    += detector[i].EventCharge*detector[i].xPosition;
            columnsum += detector[i].EventCharge*detector[i].yPosition;
            totalUV   += detector[i].EventCharge;
          }
        }
        for(unsigned int a = 0; a < detector[TreeTriggerChannel].relevantForW.size() ; a++) // for as many detectors as the ones set to relevantForW
        {
          if(detector[i].digitizerChannel == detector[TreeTriggerChannel].relevantForW[a]) // if the detector is among the ones that are relevantForUV, use them for W calculation
          {
            totalW    += detector[i].EventCharge;
          }
        }
      }

      // finally, compute u,v,w
      // near channels vs. total channels depending on what decided before
      TreeFloodX = rowsum/totalUV;
      TreeFloodY = columnsum/totalUV;
      TreeFloodZ =  maxCharge/totalW;
      if(usingTaggingBench) TreeZPosition = taggingPosition;
      if(usingRealSimData)
      {
        TreeRealX = RealX;
        TreeRealY = RealY;
        TreeRealZ = RealZ;
        TreeCrystalsHit = CrystalsHit;
        TreeNumbOfInteractions = NumbOfInteractions;
        TreeTotalEnergyDeposited = TotalEnergyDeposited;
      }
      ftree->Fill();
    }
    inFile.close();
  }
}





void InputFile::PrepareTTree()
{
  //set branches also for the analysis ttree
  //branches of the channels data
  if(electronics == 0) //CAEN x740
  {
    ftree->Branch("ExtendedTimeTag",&TreeExtendedTimeTag,"ExtendedTimeTag/l");
    ftree->Branch("DeltaTimeTag",&TreeDeltaTimeTag,"DeltaTimeTag/l");
    for (int i = 0 ; i < inputChannels ; i++)
    {
      std::stringstream sname,stype;
      sname << "ch" << detector[i].digitizerChannel;
      stype << "ch" << detector[i].digitizerChannel << "/S";
      ftree->Branch(sname.str().c_str(),&TreeCAENx740AdcChannel[i],stype.str().c_str());
    }
    ftree->Branch("TriggerChannel",&TreeTriggerChannel,"TriggerChannel/I");
    ftree->Branch("FloodX",&TreeFloodX,"FloodX/F");
    ftree->Branch("FloodY",&TreeFloodY,"FloodY/F");
    ftree->Branch("FloodZ",&TreeFloodZ,"FloodZ/F");
    if(usingTaggingBench)
    {
      ftree->Branch("Tagging",&TreeTagging,"Tagging/S");
      ftree->Branch("ZPosition",&TreeZPosition,"TreeZPosition/F");
    }
    ftree->Branch("BadEvent",&TreeBadevent,"BadEvent/O");
  }
  else if(electronics == 1) // NINO
  {
    for (int i = 0 ; i < inputChannels ; i++)
    {
      std::stringstream sname,stype;
      sname << "ch" << detector[i].digitizerChannel;
      stype << "ch" << detector[i].digitizerChannel << "/F";
      ftree->Branch(sname.str().c_str(),&TreeNINOChargeChannel[i],stype.str().c_str());
      // sname.str("");
      // stype.str("");
    }
    for (int i = 0 ; i < inputChannels ; i++)
    {
      std::stringstream sname,stype;
      sname << "t" << detector[i].digitizerChannel;
      stype << "t" << detector[i].digitizerChannel << "/F";
      ftree->Branch(sname.str().c_str(),&TreeNINOTimeChannel[i],stype.str().c_str());
    }
    ftree->Branch("TriggerChannel",&TreeTriggerChannel,"TriggerChannel/I");
    ftree->Branch("FloodX",&TreeFloodX,"FloodX/F");
    ftree->Branch("FloodY",&TreeFloodY,"FloodY/F");
    ftree->Branch("FloodZ",&TreeFloodZ,"FloodZ/F");
    if(usingTaggingBench)
    {
      ftree->Branch("Tagging",&TreeNINOtaggingCharge,"Tagging/F");
      ftree->Branch("t_Tagging",&TreeNINOtaggingTime,"t_Tagging/F");
      ftree->Branch("ZPosition",&TreeZPosition,"TreeZPosition/F");
    }
  }
  if(usingRealSimData) // sim data have same type of CAEN x740, plus other data //TODO decide how to integrate time. in principle we should use the NINO format
  {
    // pTreeTotalCryEnergy = &TreeTotalCryEnergy;
    ftree->Branch("RealX",&TreeRealX,"RealX/F");
    ftree->Branch("RealY",&TreeRealY,"RealY/F");
    ftree->Branch("RealZ",&TreeRealZ,"RealZ/F");
    ftree->Branch("CrystalsHit",&TreeCrystalsHit,"CrystalsHit/S");
    ftree->Branch("NumbOfInteractions",&TreeNumbOfInteractions,"NumbOfInteractions/S");
    ftree->Branch("TotalEnergyDeposited",&TreeTotalEnergyDeposited,"TotalEnergyDeposited/F");
    // ftree->Branch("TotalCryEnergy","std::vector<float>",&pTreeTotalCryEnergy);
  }
}






// Runs on the input TChain elements
// and fills the analysis TTree
void InputFile::FillTreeCAENx740()
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



  for (Int_t i=0;i<nevent;i++)
      // for(Int_t i=0;i<2;i++)
  {
    //loop on all the entries of tchain
    fchain->GetEvent(i);              //read complete accepted event in memory
    double maxCharge = 0;
    // double secondCharge = 0;
    float columnsum= 0;
    float rowsum= 0;
    float totalUV=  0;
    float totalW=  0;
    // float totalForFloodZ = 0;
    TreeBadevent = false;
    TreeExtendedTimeTag = ChainExtendedTimeTag;
    TreeDeltaTimeTag = ChainDeltaTimeTag;

    int TreeEntryCounter = 0;
    //loop to fill the channel
    for (int j = 0 ; j < adcChannels ; j++)
    {
      if(DigitizerChannelOn[j])
      {
        // 	std::cout << ChainAdcChannel[j] << " " ;
        // fill tree with data from the channels
        // also correcting for saturation if it's set in the config file
        if(correctingSaturation)
        {
          if(TreeCAENx740AdcChannel[TreeEntryCounter] > detector[TreeEntryCounter].saturation[0])
          {
            TreeBadevent = true;
            // 	    std::cout << "BadCharge " << (int)round(-Input[i].param0 * TMath::Log(1.0 - ( charge[i]/Input[i].param0 ))) << std::endl;
          }
          TreeCAENx740AdcChannel[TreeEntryCounter] = (Short_t)round(-detector[TreeEntryCounter].saturation[0] * TMath::Log(1.0 - ( CAENx740AdcChannel[j]/detector[TreeEntryCounter].saturation[0] )));

        }
        else
          TreeCAENx740AdcChannel[TreeEntryCounter] = CAENx740AdcChannel[j];

        detector[TreeEntryCounter].EventCharge = TreeCAENx740AdcChannel[TreeEntryCounter];
        //find the max charge and therefore the TriggerChannel
        if (TreeCAENx740AdcChannel[TreeEntryCounter] > maxCharge)
        {
          maxCharge = TreeCAENx740AdcChannel[TreeEntryCounter];
          TreeTriggerChannel = digitizer[TreeEntryCounter];
        }
        TreeEntryCounter++;
      }

      if(usingTaggingBench)
      {
        if( j == taggingCrystalChannel)
        {
          TreeTagging = CAENx740AdcChannel[j]; // no saturation correction for the tagging crystal..
          //this is the tagging crystal data
        }
      }
    }
    //     std::cout << std::endl;

    for(unsigned int i = 0 ; i < detector.size() ; i++)
    {
      for(unsigned int a = 0; a < detector[TreeTriggerChannel].relevantForUV.size() ; a++)
      {
        if(detector[i].digitizerChannel == detector[TreeTriggerChannel].relevantForUV[a])
        {
          rowsum    += detector[i].EventCharge*detector[i].xPosition;
          columnsum += detector[i].EventCharge*detector[i].yPosition;
          totalUV   += detector[i].EventCharge;
        }
      }
      for(unsigned int a = 0; a < detector[TreeTriggerChannel].relevantForW.size() ; a++)
      {
        if(detector[i].digitizerChannel == detector[TreeTriggerChannel].relevantForW[a])
        {
          totalW    += detector[i].EventCharge;
        }
      }
    }

    //compute u,v,w
    // near channels vs. total channels depending on what decided before
    TreeFloodX = rowsum/totalUV;
    TreeFloodY = columnsum/totalUV;
    TreeFloodZ =  maxCharge/totalW;
    //     TreeFloodZ =  maxCharge/totalForFloodZ;

    if(usingTaggingBench) TreeZPosition = taggingPosition;

    // TreeTheta = std::acos(TreeFloodZ /( std::sqrt( std::pow(TreeFloodX - 2.8,2) + std::pow(TreeFloodY - (-1.0),2) + std::pow(TreeFloodZ,2)) ));
    // TreePhi =  std::atan (TreeFloodY / TreeFloodX);

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


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to take the Elements created in the main file and fill them with information (hierarchy, positions, names, etc..) //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void InputFile::FillElements(Module*** module,Mppc*** mppc,Crystal*** crystal)
{

  // changing i and j, no more from top left and y x, but x y from bottom left
  std::stringstream sname;
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      sname << "Module " << iModule << "." << jModule;
      // std::cout << sname.str() << std::endl;
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
          int mppcI = (iModule*nmppcx)+iMppc;
          int mppcJ = (jModule*nmppcy)+jMppc;
          // compute the position of this mppc element from the plotPositions taken from config file
          //
          //  	  int pos = (mppcI*nmppcx + mppcJ) +1 ; // +1 is because of root multicanvas starting from 1...
          int pos = 1 + (nmppcx * nmppcy) - (mppcJ * nmppcx) - ((nmppcy-mppcI));
          // 	  std::cout << mppcI << "." << mppcJ << " " << pos << std::endl;
          int posID = -1;
          for(unsigned int arrayCounter = 0 ; arrayCounter < digitizer.size() ; arrayCounter++) // get which input mppc goes here...
          {
            if(plotPositions[arrayCounter] == pos) posID = arrayCounter;
          }

          if(posID != -1)
          {
            sname << "MPPC " << mppcI << "." << mppcJ;
            // std::cout << sname.str() << std::endl;
            mppc[mppcI][mppcJ] = new Mppc();
            mppc[mppcI][mppcJ]->SetName(sname.str().c_str());
            mppc[mppcI][mppcJ]->SetLabel( mppc_label[posID] ); // ----
            sname.str("");
            sname << iModule << "." << jModule << "-" << iMppc << "." << jMppc;

            mppc[mppcI][mppcJ]->SetIsOnForModular(true);
            for(unsigned int modCounter = 0; modCounter < mppcOFF.size(); modCounter++)
            {
              if(mppcOFF[modCounter].compare(mppc_label[posID]) == 0) mppc[mppcI][mppcJ]->SetIsOnForModular(false);
            }


            //set the global mppc variables for 3d plots, then override if specified
            mppc[mppcI][mppcJ]->SetHisto3DchannelBin(global_histo3DchannelBin);
            mppc[mppcI][mppcJ]->SetClusterLevelPrecision(global_div);
            mppc[mppcI][mppcJ]->SetClusterVolumeCut(global_clusterVolumeCut);
            //now override
            for(unsigned int modCounter = 0; modCounter < specificMPPC.size(); modCounter++)
            {
              if(specificMPPC[modCounter].compare(mppc_label[posID]) == 0)
              {
                mppc[mppcI][mppcJ]->SetHisto3DchannelBin(specificBin[modCounter]);
                mppc[mppcI][mppcJ]->SetClusterLevelPrecision(specificPrecision[modCounter]);
                mppc[mppcI][mppcJ]->SetClusterVolumeCut(specificCut[modCounter]);
              }
            }

            mppc[mppcI][mppcJ]->SetExtendedID(sname.str().c_str());
            mppc[mppcI][mppcJ]->SetID(posID);                  // ----
            mppc[mppcI][mppcJ]->SetI(mppcI);
            mppc[mppcI][mppcJ]->SetJ(mppcJ);
            mppc[mppcI][mppcJ]->SetChildrenI(ncrystalsx);
            mppc[mppcI][mppcJ]->SetChildrenJ(ncrystalsy);
            mppc[mppcI][mppcJ]->SetPosition(xPositions[posID],yPositions[posID],0);
            mppc[mppcI][mppcJ]->SetDigitizerChannel(digitizer[posID]);
            mppc[mppcI][mppcJ]->SetCanvasPosition(pos);
            mppc[mppcI][mppcJ]->SetParentName(module[iModule][jModule]->GetName());
            for(unsigned int iDet = 0 ; iDet < detector.size() ; iDet++)
            {
              if(detector[iDet].digitizerChannel == digitizer[posID])
              {
                mppc[mppcI][mppcJ]->SetRelevantForUV(detector[iDet].relevantForUV);
                mppc[mppcI][mppcJ]->SetRelevantForW(detector[iDet].relevantForW);
                mppc[mppcI][mppcJ]->SetRelevantForE(detector[iDet].relevantForE);
              }

            }

            mppc[mppcI][mppcJ]->SetIsOnForDoi(false);
            for(unsigned int iDoi = 0 ; iDoi < digitizerDoi.size(); iDoi++)
            {
              if(mppc[mppcI][mppcJ]->GetDigitizerChannel() == digitizerDoi[iDoi])
              mppc[mppcI][mppcJ]->SetIsOnForDoi(true);
            }

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
                // std::cout << sname.str() << std::endl;
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
                  xPositions[posID] + iCry*(crystalx+esrThickness) - (ncrystalsx-1)*((crystalx+esrThickness)/2.0)
                  , yPositions[posID] + jCry*(crystaly+esrThickness) - (ncrystalsy-1)*((crystaly+esrThickness)/2.0)
                  ,  0  ); //FIXME z not useful so set to 0, but maybe we should put the real one..
                  crystal[cryI][cryJ]->SetDimension(crystalx,crystaly,crystalz);
                  // 	      double u  = crystaldata[cryI][cryJ][0];
                  // 	      double v  = crystaldata[cryI][cryJ][1];
                  // 	      double wu = crystaldata[cryI][cryJ][2];
                  // 	      double wv = crystaldata[cryI][cryJ][3];
                  // 	      double t  = crystaldata[cryI][cryJ][4];
                  // 	      crystal[cryI][cryJ]->SetCrystalData(u,v,wu,wv,t);
                  // 	      std::cout << u << " " << v << " " << wu << " " << wv << " " << t << std::endl;
                  // 	      crystal[cryI][cryJ]->SetEllipses(u,v,wu,wv,t);
                  // 	      TEllipse *ellipse = new TEllipse(u,v,wu,wv,0,360,t);
                  // 	      crystal[cryI][cryJ]->SetGraphicalCut(*ellipse);
                  // 	      crystal[cryI][cryJ]->Print();

                  mppc[mppcI][mppcJ]->AddChild( crystal[cryI][cryJ]->GetName() );
                }
              }
            }
            else //set them off for modular
            {
              mppc[mppcI][mppcJ] = new Mppc();
              mppc[mppcI][mppcJ]->SetIsOnForModular(false);
              for(int iCry = 0; iCry < ncrystalsx ; iCry++)
              {
                for(int jCry = 0; jCry < ncrystalsy ; jCry++)
                {
                  int cryI = (iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry);
                  int cryJ = (jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry);
                  crystal[cryI][cryJ] = new Crystal();
                  crystal[cryI][cryJ]->SetIsOnForModular(false);
                }
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
      if(mppc[iMppc][jMppc]->GetIsOnForModular()) std::cout <<  mppc[iMppc][jMppc]->GetLabel() << "\t";
      else std::cout << "\t";
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
      if(mppc[iMppc][jMppc]->GetIsOnForModular()) std::cout << mppc[iMppc][jMppc]->GetDigitizerChannel() << "\t";
      else std::cout << "\t";
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
      if(crystal[iCrystal][jCrystal]->GetIsOnForModular()) std::cout << crystal[iCrystal][jCrystal]->GetID() << "\t";
      else std::cout << "\t";
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
      if(crystal[iCrystal][jCrystal]->GetIsOnForModular()) std::cout << "(" << crystal[iCrystal][jCrystal]->GetX() << "," << crystal[iCrystal][jCrystal]->GetY() << ")"  << "\t";
      else std::cout << "\t\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to read the config file and set the keys                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void InputFile::ReadConfig(ConfigFile& config)
{
   //read configuration file
   fname                       = config.read<std::string>("chainName","adc");
   ncrystalsx                  = config.read<int>("ncrystalsx",1);
   ncrystalsy                  = config.read<int>("ncrystalsy",1);
   nmppcx                      = config.read<int>("nmppcx",1);
   nmppcy                      = config.read<int>("nmppcy",1);
   nmodulex                    = config.read<int>("nmodulex",1);
   nmoduley                    = config.read<int>("nmoduley",1);
   taggingPosition             = config.read<float>("taggingPosition",0);
   usingTaggingBench           = config.read<bool>("usingTaggingBench",0);
   taggingCrystalChannel       = config.read<int>("taggingCrystalChannel",16);
   usingRealSimData            = config.read<bool>("usingRealSimData",0);
   binary                      = config.read<bool>("binary",0);
   correctingSaturation        = config.read<bool>("correctingSaturation",0);
   BinaryOutputFileName        = config.read<std::string>("output","binOutput");
   BinaryOutputFileName       += ".bin";
   //read the strings that describe the input channels
   digitizer_s                 = config.read<std::string>("digitizer");     //MANDATORY, so no default
   mppc_s                      = config.read<std::string>("mppc");          //MANDATORY, so no default
   plotPositions_s             = config.read<std::string>("plotPositions"); //MANDATORY, so no default
   xPositions_s                = config.read<std::string>("xPositions");    //MANDATORY, so no default
   yPositions_s                = config.read<std::string>("yPositions");    //MANDATORY, so no default
   saturation_s                = config.read<std::string>("saturation","0");
   adcChannels                 = config.read<int>("digitizerTotalCh",32);
   nclock                      = config.read<double>("nclock",0);
   crystalx                    = config.read<double>("crystalx",1.53);
   crystaly                    = config.read<double>("crystaly",1.53);
   crystalz                    = config.read<double>("crystalz",15);
   esrThickness                = config.read<double>("esrThickness",0.07);
   usingAllChannels            = config.read<bool>("usingAllChannels",0);   // whether to use the sum of all channels to compute u,v or just the neighbours. Deafult to false = 0.
   wAllChannels                = config.read<bool>("wAllChannels",0);       // whether we use the sum of all channels to compute w ( of just the neighbours. Deafult to false = 0.
   electronics                 = config.read<int>("electronics",0) ;
   ninoSaturationParameters    = config.read<int>("ninoSaturationParameters",4) ;
   digitizerDoi_s              = config.read<std::string>("digiChannelsForDoi","8,9,10,11");
   mppcOFF_s                   = config.read<std::string>("mppcOFF","");
   crystalOFF_s                = config.read<std::string>("crystalOFF","-1");
   // global 3d plots variables for single mppcs
   global_histo3DchannelBin    = config.read<int>("histo3DchannelBin",100);
   global_div                  = config.read<int>("clusterLevelPrecision",10);
   global_clusterVolumeCut     = config.read<double>("clusterVolumeCut",0.001);
   specificMPPC_s              = config.read<std::string>("specificMPPCname","");
   specificBin_s               = config.read<std::string>("specificBin","");
   specificPrecision_s         = config.read<std::string>("specificPrecision","");
   specificCut_s               = config.read<std::string>("specificCut","");
   // u-v-w and e computation
   // choice of channels involved in calculation of u-v-w and total energy
   // Possibilities:
   // 0 = all channels in the mppc array
   // 1 = only the trigger channel + neighbour channels (DEFAULT for all 3)
   // 2 = only the trigger channel + cross channels
   relevantForUV               = config.read<int>("relevantForUV",1);
   relevantForW                = config.read<int>("relevantForW",1);
   relevantForE                = config.read<int>("relevantForE",1);

   //PARSE STRINGS
   //split them using the config file class
   config.split( digitizer_f, digitizer_s, "," );
   config.split( mppc_f, mppc_s, "," );
   config.split( plotPositions_f, plotPositions_s, "," );
   config.split( xPositions_f, xPositions_s, "," );
   config.split( yPositions_f, yPositions_s, "," );
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
   //check if the vectors just built have the same size
   assert( (digitizer.size() == mppc_label.size() ) && (digitizer.size() == plotPositions.size()) && (digitizer.size() == xPositions.size()) && (digitizer.size() == yPositions.size()) );


   //saturation part. different for caen and NINO, but only difference is that CAEN takes 1 saturation param, NINO 4. so the asset is the only difference
   if(saturation_s.compare("0") != 0)
   {
     correctingSaturation = true;
     config.split( saturation_f, saturation_s, "," );
     for(unsigned int i = 0 ; i < saturation_f.size() ; i++)
     {
       config.trim(saturation_f[i]);
       saturation.push_back(atof(saturation_f[i].c_str()));
     }
     if(electronics == 0) assert( (digitizer.size() == saturation.size()) );   //CAEN
     else if(electronics == 1) assert( (digitizer.size()*ninoSaturationParameters == saturation.size()) ); //NINO
   }

   //read string for doi analysis channels
   // digitizerDoi_s = config.read<std::string>("digiChannelsForDoi","8,9,10,11");
   config.split( digitizerDoi_f, digitizerDoi_s, "," );
   for(unsigned int i = 0 ; i < digitizerDoi_f.size() ; i++)
   {
     config.trim(digitizerDoi_f[i]);
     digitizerDoi.push_back(atoi(digitizerDoi_f[i].c_str()));
   }
   //on or off for modular analysis
   config.split( mppcOFF_f, mppcOFF_s, "," );
   for(unsigned int i = 0 ; i < mppcOFF_f.size() ; i++)
   {
     config.trim(mppcOFF_f[i]);
     mppcOFF.push_back(mppcOFF_f[i]);
   }

   config.split( crystalOFF_f, crystalOFF_s, "," );
   for(unsigned int i = 0 ; i < crystalOFF_f.size() ; i++)
   {
     config.trim(crystalOFF_f[i]);
     crystalOFF.push_back(atoi(crystalOFF_f[i].c_str()));
   }
   //specific variables for mppcs
   config.split( specificMPPC_f, specificMPPC_s, "," );
   for(unsigned int i = 0 ; i < specificMPPC_f.size() ; i++)
   {
     config.trim(specificMPPC_f[i]);
     specificMPPC.push_back(specificMPPC_f[i]);
   }

   config.split( specificBin_f, specificBin_s, "," );
   for(unsigned int i = 0 ; i < specificBin_f.size() ; i++)
   {
     config.trim(specificBin_f[i]);
     specificBin.push_back(atoi(specificBin_f[i].c_str()));
   }

   config.split( specificPrecision_f, specificPrecision_s, "," );
   for(unsigned int i = 0 ; i < specificPrecision_f.size() ; i++)
   {
     config.trim(specificPrecision_f[i]);
     specificPrecision.push_back(atoi(specificPrecision_f[i].c_str()));
   }
   config.split( specificCut_f, specificCut_s, "," );
   for(unsigned int i = 0 ; i < specificCut_f.size() ; i++)
   {
     config.trim(specificCut_f[i]);
     specificCut.push_back(atof(specificCut_f[i].c_str()));
   }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to find the neighbouring channels of the channel specified by TreeTriggerChannel                               //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int> InputFile::FindNeighbours(int TreeTriggerChannel)
{
  std::vector<int> channels;

  double xTrigger = 0;
  double yTrigger = 0;
  std::vector<float> allowedX,allowedY;
  for(int iFill = 0; iFill < inputChannels; iFill++)
  {
    detector[iFill].isNeighbour = false;
  }

  std::vector<float> xCopyTemp;
  std::vector<float> yCopyTemp;
  for(int iFill = 0; iFill < inputChannels; iFill++)
  {
    if(detector[iFill].digitizerChannel == TreeTriggerChannel)
    {
      xTrigger = detector[iFill].xPosition;
      yTrigger = detector[iFill].yPosition;
    }
    xCopyTemp.push_back(detector[iFill].xPosition);
    yCopyTemp.push_back(detector[iFill].yPosition);
  }

  //sort the position vectors
  std::sort (xCopyTemp.begin(),xCopyTemp.end());
  std::sort (yCopyTemp.begin(),yCopyTemp.end());
  //filter repetitions
  std::vector<float> xCopy;
  std::vector<float> yCopy;
  xCopy.push_back(xCopyTemp[0]);
  yCopy.push_back(yCopyTemp[0]);
  int xCounter = 0;
  int yCounter = 0;
  for (unsigned int j = 0 ; j < xCopyTemp.size() ; j++)
  {
    if(xCopy[xCounter] != xCopyTemp[j])// this will work only because xCopyTemp is already sorted...
    {
      xCopy.push_back(xCopyTemp[j]);
      xCounter++;
    }
  }
  for (unsigned int j = 0 ; j < yCopyTemp.size() ; j++)
  {
    if(yCopy[yCounter] != yCopyTemp[j])// this will work only because yCopyTemp is already sorted...
    {
      yCopy.push_back(yCopyTemp[j]);
      yCounter++;
    }
  }

  // by default trigger row and column are allowed
  allowedX.push_back(xTrigger);
  allowedY.push_back(yTrigger);

  int TrigI = -1;
  int TrigJ = -1;
  //locate trigger in the new vectors
  for(unsigned int xPos = 0 ; xPos < xCopy.size() ; xPos++)
  {
    if(xCopy[xPos] == xTrigger) TrigI = xPos;
  }
  for(unsigned int yPos = 0 ; yPos < yCopy.size() ; yPos++)
  {
    if(yCopy[yPos] == yTrigger) TrigJ = yPos;
  }

  // take x and y of neighbours
  if(TrigI != 0)
  {
    allowedX.push_back(xCopy[TrigI-1]);
  }
  if(TrigI != nmppcx-1)
  {
    allowedX.push_back(xCopy[TrigI+1]);
  }
  if(TrigJ != 0)
  {
    allowedY.push_back(yCopy[TrigJ-1]);
  }
  if(TrigJ != nmppcy-1)
  {
    allowedY.push_back(yCopy[TrigJ+1]);
  }
  int counterNeighbour = 0;
  for (int j = 0 ; j < adcChannels ; j++)
  {
    if(DigitizerChannelOn[j])
    {
      for(unsigned int iCheck = 0; iCheck < allowedX.size(); iCheck++)
      {
        if(detector[counterNeighbour].xPosition == allowedX[iCheck]) //check if x is allowed
        {
          for(unsigned int jCheck = 0; jCheck < allowedY.size(); jCheck++)
          {
            if(detector[counterNeighbour].yPosition == allowedY[jCheck]) //check if y is allowed
            {
              detector[counterNeighbour].isNeighbour = true;
              channels.push_back(detector[counterNeighbour].digitizerChannel);
            }
          }
        }
      }
      counterNeighbour++;
    }
  }
  return channels;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to find the channels above, below, on the right and on the left of the channel specified by TreeTriggerChannel //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int> InputFile::FindCross(int TreeTriggerChannel)
{
  std::vector<int> channels;

  double xTrigger = 0;
  double yTrigger = 0;
  std::vector<float> allowedX,allowedY;
  for(int iFill = 0; iFill < inputChannels; iFill++)
  {
    detector[iFill].isCross = false;
  }

  std::vector<float> xCopyTemp;
  std::vector<float> yCopyTemp;
  for(int iFill = 0; iFill < inputChannels; iFill++)
  {
    if(detector[iFill].digitizerChannel == TreeTriggerChannel)
    {
      xTrigger = detector[iFill].xPosition;
      yTrigger = detector[iFill].yPosition;
    }
    xCopyTemp.push_back(detector[iFill].xPosition);
    yCopyTemp.push_back(detector[iFill].yPosition);
  }

  //sort the position vectors
  std::sort (xCopyTemp.begin(),xCopyTemp.end());
  std::sort (yCopyTemp.begin(),yCopyTemp.end());
  //filter repetitions
  std::vector<float> xCopy;
  std::vector<float> yCopy;
  xCopy.push_back(xCopyTemp[0]);
  yCopy.push_back(yCopyTemp[0]);
  int xCounter = 0;
  int yCounter = 0;
  for (unsigned int j = 0 ; j < xCopyTemp.size() ; j++)
  {
    if(xCopy[xCounter] != xCopyTemp[j])// this will work only because xCopyTemp is already sorted...
    {
      xCopy.push_back(xCopyTemp[j]);
      xCounter++;
    }
  }
  for (unsigned int j = 0 ; j < yCopyTemp.size() ; j++)
  {
    if(yCopy[yCounter] != yCopyTemp[j])// this will work only because yCopyTemp is already sorted...
    {
      yCopy.push_back(yCopyTemp[j]);
      yCounter++;
    }
  }

  // by default trigger row and column are allowed
  allowedX.push_back(xTrigger);
  allowedY.push_back(yTrigger);

  int TrigI = -1;
  int TrigJ = -1;
  //locate trigger in the new vectors
  for(unsigned int xPos = 0 ; xPos < xCopy.size() ; xPos++)
  {
    if(xCopy[xPos] == xTrigger) TrigI = xPos;
  }
  for(unsigned int yPos = 0 ; yPos < yCopy.size() ; yPos++)
  {
    if(yCopy[yPos] == yTrigger) TrigJ = yPos;
  }

  // take x and y of neighbours
  if(TrigI != 0)
  {
    allowedX.push_back(xCopy[TrigI-1]);
  }
  if(TrigI != nmppcx-1)
  {
    allowedX.push_back(xCopy[TrigI+1]);
  }
  if(TrigJ != 0)
  {
    allowedY.push_back(yCopy[TrigJ-1]);
  }
  if(TrigJ != nmppcy-1)
  {
    allowedY.push_back(yCopy[TrigJ+1]);
  }
  int counterNeighbour = 0;
  for (int j = 0 ; j < adcChannels ; j++)
  {
    if(DigitizerChannelOn[j])
    {
      if(detector[counterNeighbour].xPosition == xTrigger && detector[counterNeighbour].yPosition == yTrigger) //trigger is part of the cross
      {
        detector[counterNeighbour].isCross = true;
        channels.push_back(detector[counterNeighbour].digitizerChannel);
      }
      for(unsigned int iCheck = 0; iCheck < allowedX.size(); iCheck++)  // allowed X, only if y is trigger y
      {
        if(detector[counterNeighbour].xPosition == allowedX[iCheck] && detector[counterNeighbour].xPosition != xTrigger) //check if x is allowed
        {
          if(detector[counterNeighbour].yPosition == yTrigger)
          {
            detector[counterNeighbour].isCross = true;
            channels.push_back(detector[counterNeighbour].digitizerChannel);
          }
        }
      }
      for(unsigned int jCheck = 0; jCheck < allowedY.size(); jCheck++) // allowed Y, only if x is trigger x
      {
        if(detector[counterNeighbour].yPosition == allowedY[jCheck] && detector[counterNeighbour].yPosition != yTrigger) //check if y is allowed
        {
          if(detector[counterNeighbour].xPosition == xTrigger )
          {
            detector[counterNeighbour].isCross = true;
            channels.push_back(detector[counterNeighbour].digitizerChannel);
          }
        }
      }
      counterNeighbour++;
    }
  }
  return channels;
}
