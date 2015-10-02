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

struct Point 
{
  float x;
  float y;
  float z;
} __attribute__((__packed__));     


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
  taggingPosition             = config.read<float>("taggingPosition");
  usingTaggingBench           = config.read<bool>("usingTaggingBench");
  taggingCrystalChannel       = config.read<int>("taggingCrystalChannel");
  usingRealSimData            = config.read<bool>("usingRealSimData");
  binary                      = config.read<bool>("binary");
  correctingSaturation        = config.read<bool>("correctingSaturation");
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
  
  //read strings that describes crystals
  
  
  //a string for input for each crystal
  crystal_s = new std::string*[ncrystalsx*nmppcx*nmodulex];
  for(int i = 0 ; i < ncrystalsx*nmppcx*nmodulex ; i++) crystal_s[i] = new std::string[ncrystalsy*nmppcy*nmoduley];
  //crystal on
  crystalIsOn = new bool*[ncrystalsx*nmppcx*nmodulex];
  for(int i = 0 ; i < ncrystalsx*nmppcx*nmodulex ; i++) 
  {
    crystalIsOn[i] = new bool[ncrystalsy*nmppcy*nmoduley];
    for(int j = 0 ; j < ncrystalsy*nmppcy*nmoduley ; j++) 
    {
      crystalIsOn[i][j] = false;
    }
  }
  //a float vector for each crystal
  crystaldata = new float**[ncrystalsx*nmppcx*nmodulex];
  for(int i = 0 ; i < ncrystalsx*nmppcx*nmodulex ; i++) 
  {
    crystaldata[i] = new float*[ncrystalsy*nmppcy*nmoduley];
    for(int j = 0 ; j < ncrystalsy*nmppcy*nmoduley ; j++) 
    { 
      crystaldata[i][j] = new float[5];
      for (int k = 0 ; k < 5 ; k++)
      {
	crystaldata[i][j][k] = 0;
      }
    }
  }
  
  int crystalCounter = 0;
  for(int ii = 0; ii < ncrystalsx*nmppcx*nmodulex*ncrystalsy*nmppcy*nmoduley ; ii++)
  {
    std::stringstream crystalstring;
    crystalstring << "crystal" << crystalCounter;
    std::string tempString;
    std::vector<std::string> tempStringVector;
    tempString = config.read<std::string>(crystalstring.str().c_str(),"0,0,0,0,0,0,0"); //FIXME i and j here are not the i and j set on the crystals!!! actually, in the loop for the crystals, below, it's wrong as well, as we don't consider the multiple modules case!! it will work for one module, it has to be fixed for multiple ones.
    
//     std::cout << crystalCounter << " " << crystalstring.str() << " " << tempString /*<< std::endl*/;
    
    config.split( tempStringVector, tempString, "," );
    for(int i = 0 ; i < tempStringVector.size() ; i++)
    {
      config.trim(tempStringVector[i]);
      // 	crystaldata.push_back(atof(crystal_f[i].c_str()));
    }
    int CryIDi = atoi(tempStringVector[0].c_str());
    int CryIDj = atoi(tempStringVector[1].c_str());
    
    for(int i = 0 ; i < 5 ; i++)
    {
      crystaldata[CryIDi][CryIDj][i] = atof(tempStringVector[i+2].c_str());
//       std::cout << " " <<  crystaldata[CryIDi][CryIDj][i] ;
    }
    
    if(tempString != "0,0,0,0,0,0,0")
      crystalIsOn[CryIDi][CryIDj] = true;
    
//     std::cout << " " << crystalIsOn[CryIDi][CryIDj] << std::endl;
    
    crystalCounter++;
  }
  
  //   crystalu_s              = config.read<std::string>("crystalu");
  //   crystalv_s              = config.read<std::string>("crystalv");
  //   crystalwu_s             = config.read<std::string>("crystalwu");
  //   crystalwv_s             = config.read<std::string>("crystalwv");
  //   crystalt_s              = config.read<std::string>("crystalt");
  //   config.split( crystalu_f, crystalu_s, "," );
  //   config.split( crystalv_f, crystalv_s, "," );
  //   config.split( crystalwu_f, crystalwu_s, "," );
  //   config.split( crystalwv_f, crystalwv_s, "," );
  //   config.split( crystalt_f, crystalt_s, "," );
  //   for(int i = 0 ; i < crystalu_f.size() ; i++)
  //   {
  //     config.trim(crystalu_f[i]);
  //     crystalu.push_back(atof(crystalu_f[i].c_str()));
  //   }
  //   for(int i = 0 ; i < crystalv_f.size() ; i++)
  //   {
  //     config.trim(crystalv_f[i]);
  //     crystalv.push_back(atof(crystalv_f[i].c_str()));
  //   }
  //   for(int i = 0 ; i < crystalwu_f.size() ; i++)
  //   {
  //     config.trim(crystalwu_f[i]);
  //     crystalwu.push_back(atof(crystalwu_f[i].c_str()));
  //   }
  //   for(int i = 0 ; i < crystalwv_f.size() ; i++)
  //   {
  //     config.trim(crystalwv_f[i]);
  //     crystalwv.push_back(atof(crystalwv_f[i].c_str()));
  //   }
  //   for(int i = 0 ; i < crystalt_f.size() ; i++)
  //   {
  //     config.trim(crystalt_f[i]);
  //     crystalt.push_back(atof(crystalt_f[i].c_str()));
  //   }
  
  
  
  
  if(digitizer.size() > 16) 
  {
    std::cout << "ERROR: Only one module can be analyzed at a time! Set 16 or less input channels in the config file!" << std::endl;
  }
  //feedback to the user
  std::cout << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << " Channels configuration " << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << "ADC input\tMPPC ch\tCanvas\tx[mm]\ty[mm]" << std::endl;
  std::cout << "------------------------" << std::endl;
  for(int i = 0 ; i < digitizer.size() ; i++)
  {
    std::cout << "Channel[" << digitizer[i] << "] = \t" <<  mppc_label[i] << "\t" << plotPositions[i] << "\t" << xPositions[i] << "\t" << yPositions[i] << std::endl;
  }
  std::cout << "------------------------" << std::endl;
  std::cout << std::endl;
  
  //crystal ON for analysis
  //   crystalIsOn = new bool* [nmodulex*nmppcx*ncrystalsx]; 
  //   for(int j = 0; j < nmodulex*nmppcx*ncrystalsx ; j++) crystalIsOn[j] = new bool [nmoduley*nmppcy*ncrystalsy];
  //   for(int i = 0; i < nmodulex*nmppcx*ncrystalsx; i++)
  //     for(int j = 0; j < nmoduley*nmppcy*ncrystalsy; j++)
  //       crystalIsOn = false;
  
  
  
  
  
  //------------------------------------------------------------------------------------------//
  //  opens the Tchain, set its branches, create the TTree that will be used for the analysis //
  //------------------------------------------------------------------------------------------//  
  inputChannels = digitizer.size();
  fchain              = new TChain(fname.c_str());  // create the input tchain and the analysis ttree
  ftree               = new TTree(fname.c_str(),fname.c_str());
  // first, create the adc channels variables and branches
  ChainAdcChannel     = new Short_t [adcChannels]; // input from ADC is always 32 ch for CAEN - FIXME at some point if necessary...
  bChainAdcChannel    = new TBranch* [adcChannels];
  TreeAdcChannel      = new Short_t [inputChannels]; // channels analyzed instead can be up to 16
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
  }
  for(int i=0; i<adcChannels; i++)
  {
    std::stringstream sname;
    sname << "ch" << i;
    fchain->SetBranchAddress(sname.str().c_str(), &ChainAdcChannel[i], &bChainAdcChannel[i]);
  }
  //set branches also for the analysis ttree
  ftree->Branch("ExtendedTimeTag",&TreeExtendedTimeTag,"ExtendedTimeTag/l"); 
  ftree->Branch("DeltaTimeTag",&TreeDeltaTimeTag,"DeltaTimeTag/l");
  //branches of the 32 channels data
  for (int i = 0 ; i < inputChannels ; i++)
  {
    //empty the stringstreams
    std::stringstream sname,stype;
    sname << "ch" << i;
    stype << "ch" << i << "/S";  
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
    ftree->Branch("RealX",&TreeRealX,"RealX/F"); 
    ftree->Branch("RealY",&TreeRealY,"RealY/F"); 
    ftree->Branch("RealZ",&TreeRealZ,"RealZ/F"); 
  }
}

void InputFile::CreateTree()
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
    DigitizerChannelOn[digitizer[i]] = true;
  }
  
  for (Int_t i=0;i<nevent;i++) 
  { 
    //loop on all the entries of tchain
    fchain->GetEvent(i);              //read complete accepted event in memory
    double maxCharge = 0;
    double secondCharge = 0;
    float columnsum= 0;
    float rowsum= 0;
    float total=  0;
    TreeBadevent = false;
    TreeExtendedTimeTag = ChainExtendedTimeTag;
    TreeDeltaTimeTag = ChainDeltaTimeTag;
    
    int TreeEntryCounter = 0;
    
    for (int j = 0 ; j < 32 ; j++) //loop on input chain channels
    {
      if(DigitizerChannelOn[j])
      {
	// fill tree with data from the channels
	// also correcting for saturation if it's set in the config file
	if(correctingSaturation)
	{
	  if(TreeAdcChannel[TreeEntryCounter] > saturation[TreeEntryCounter])
	  {
	    TreeBadevent = true;
	    //std::cout << "BadCharge " << (int)round(-Input[i].param0 * TMath::Log(1.0 - ( charge[i]/Input[i].param0 ))) << std::endl;
	  }
	  TreeAdcChannel[TreeEntryCounter] = (Short_t)round(-saturation[TreeEntryCounter] * TMath::Log(1.0 - ( ChainAdcChannel[j]/saturation[TreeEntryCounter] )));
	}
	else
	  TreeAdcChannel[TreeEntryCounter] = ChainAdcChannel[j];
	
	//find the max charge and therefore the TriggerChannel
	if (TreeAdcChannel[TreeEntryCounter] > maxCharge)
	{
	  maxCharge = TreeAdcChannel[TreeEntryCounter];
	  TreeTriggerChannel = TreeEntryCounter;
	}
	//update total, row and columnsum
	total += TreeAdcChannel[TreeEntryCounter];
	rowsum += TreeAdcChannel[TreeEntryCounter]*xPositions[TreeEntryCounter];
	columnsum += TreeAdcChannel[TreeEntryCounter]*yPositions[TreeEntryCounter];
	TreeEntryCounter++;
      }
      
      if(usingTaggingBench)
      {
	if( j == taggingCrystalChannel)
	{
	  TreeTagging = ChainAdcChannel[j];
	  //this is the tagging crystal data
	}
      }
    }
    
    TreeFloodX = rowsum/total;
    TreeFloodY = columnsum/total;
    TreeFloodZ =  maxCharge/total;
    
    if(usingTaggingBench) TreeZPosition = taggingPosition;
    
    TreeTheta = std::acos(TreeFloodZ /( std::sqrt( std::pow(TreeFloodX,2) + std::pow(TreeFloodY,2) + std::pow(TreeFloodZ,2)) ));
    TreePhi =  std::atan (TreeFloodY / TreeFloodX);
    
    point.x = TreeFloodX;
    point.y = TreeFloodY;
    point.z = TreeFloodZ;
    
    if(usingRealSimData)
    {
      TreeRealX = RealX;
      TreeRealY = RealY;
      TreeRealZ = RealZ;
    }
    
    if(!TreeBadevent)
    {
      ftree->Fill();
      GoodCounter++;
    }
    else
    {
      badEvents++;
    }
    
    if(binary)
      output_file.write((char*)&point,sizeof(point));
    
    //counter to give a feedback to the user
    counter++;
    
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

void InputFile::FillElements(Module*** module,Mppc*** mppc,Crystal*** crystal)
{
  //temp
  int translateCh[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  int moduleCounter = 0;
  
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      std::stringstream sname;
      sname << "Module " << iModule << "." << jModule;
      module[iModule][jModule] = new Module(); // creates a default module   
      module[iModule][jModule]->SetName(sname.str().c_str());          // assign a name
      module[iModule][jModule]->SetID(moduleCounter);                  // assign an ID number
      module[iModule][jModule]->SetI(iModule); 
      module[iModule][jModule]->SetJ(jModule);
      module[iModule][jModule]->SetChildrenI(nmppcx);
      module[iModule][jModule]->SetChildrenJ(nmppcy);
      moduleCounter++;
    }
  }
  
  int mppcCounter = 0;
  
  for(int iMppc = 0; iMppc < nmppcx*nmodulex ; iMppc++)
  {
    for(int jMppc = 0; jMppc < nmppcy*nmoduley ; jMppc++)
    {
      std::stringstream sname;
      sname << "Mppc " << mppc_label[mppcCounter];
      mppc[iMppc][jMppc] = new Mppc();
      mppc[iMppc][jMppc]->SetName(sname.str().c_str());   // assign a name
      mppc[iMppc][jMppc]->SetLabel(mppc_label[mppcCounter]);
      mppc[iMppc][jMppc]->SetID(mppcCounter);             // assign an ID
      mppc[iMppc][jMppc]->SetI(iMppc); 
      mppc[iMppc][jMppc]->SetJ(jMppc);
      mppc[iMppc][jMppc]->SetChildrenI(ncrystalsx);
      mppc[iMppc][jMppc]->SetChildrenJ(ncrystalsy);
      mppc[iMppc][jMppc]->SetPosition(xPositions[mppcCounter],yPositions[mppcCounter],0);
      mppc[iMppc][jMppc]->SetDigitizerChannel(translateCh[mppcCounter]);
      mppc[iMppc][jMppc]->SetCanvasPosition(plotPositions[mppcCounter]);
      mppc[iMppc][jMppc]->SetParentName(module[iMppc/nmppcx][jMppc/nmppcy]->GetName());
      mppcCounter++;
    }
  }
  
  int crystalCounter = 0;
  
  
//   for(int iCrystal = 0; iCrystal < ncrystalsx*nmppcx*nmodulex ; iCrystal++)
//   {
//     for(int jCrystal = 0; jCrystal < ncrystalsy*nmppcy*nmoduley ; jCrystal++)
//     {
//       for (int zz = 0 ; zz < 5 ; zz++)
//       {
// 	std::cout << crystaldata[iCrystal][jCrystal][zz] << " ";
//       }
//       std::cout << std::endl;
//     }
//     
//   }
  
  
  for(int iCrystal = 0; iCrystal < ncrystalsx*nmppcx*nmodulex ; iCrystal++)
  {
    for(int jCrystal = 0; jCrystal < ncrystalsy*nmppcy*nmoduley ; jCrystal++)
    {
      std::stringstream stream;
      stream << "Crystal " << crystalCounter;
      crystal[iCrystal][jCrystal] = new Crystal();
      crystal[iCrystal][jCrystal]->SetName(stream.str().c_str());   // assign a name
      crystal[iCrystal][jCrystal]->SetID(crystalCounter);          // assign an ID
      crystal[iCrystal][jCrystal]->SetI(iCrystal); 
      crystal[iCrystal][jCrystal]->SetJ(jCrystal);
      crystal[iCrystal][jCrystal]->SetParentName(mppc[iCrystal/ncrystalsx][jCrystal/ncrystalsy]->GetName());
      crystal[iCrystal][jCrystal]->SetCrystalOn(crystalIsOn[iCrystal][jCrystal]);
      double u  = crystaldata[iCrystal][jCrystal][0];
      double v  = crystaldata[iCrystal][jCrystal][1];
      double wu = crystaldata[iCrystal][jCrystal][2];
      double wv = crystaldata[iCrystal][jCrystal][3];
      double t  = crystaldata[iCrystal][jCrystal][4];
//       std::cout << u << " " << v << " " << wu << " " << wv << " " << t << std::endl;
      crystal[iCrystal][jCrystal]->SetEllipses(u,v,wu,wv,t);
      TEllipse *ellipse = new TEllipse(u,v,wu,wv,0,360,t);
      crystal[iCrystal][jCrystal]->SetGraphicalCut(*ellipse);
//       crystal[iCrystal][jCrystal]->Print();
      crystalCounter++;
    }
  }
  
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
	{
	  module[iModule][jModule]->AddChild( mppc[(iModule * nmppcx) + iMppc][(jModule * nmppcy) + jMppc]->GetName() );
	  for(int iCrystal = 0; iCrystal < ncrystalsx ; iCrystal++)
	  {
	    for(int jCrystal = 0; jCrystal < ncrystalsy ; jCrystal++)
	    {
	      mppc[iMppc][jMppc]->AddChild( crystal[(iMppc * ncrystalsx) + iCrystal][(jMppc * ncrystalsy) + jCrystal]->GetName() );
	    }
	  }
	}
      }
    } 
  }
  
  
  
  //   std::cout << "--------------------------------------------" << std::endl;
  //   module[0][0]->Print();
  //   std::cout <<module[0][0]->GetMppcsNumber() << std::endl;
  //   for(int iModule = 0; iModule < nmodulex ; iModule++)
  //   {
  //     for(int jModule = 0; jModule < nmoduley ; jModule++)
  //     {
  //       module[iModule][jModule]->Print();
  //       for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
  //       {
  // 	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
  // 	{
  // 	  mppc[(iModule * nmppcx) + iMppc][(jModule * nmppcy) + jMppc]->Print();
  // 	  for(int iCrystal = 0; iCrystal < ncrystalsx ; iCrystal++)
  // 	  {
  // 	    for(int jCrystal = 0; jCrystal < ncrystalsy ; jCrystal++)
  // 	    {
  // 	      crystal[(iMppc * ncrystalsx) + iCrystal][(jMppc * ncrystalsy) + jCrystal]->Print();
  // 	    }
  // 	  }
  // 	}
  //       }
  //     }
  //   }
  // //   std::cout <<module[0][0]->GetMppcsNumber() << std::endl;
  // //   std::cout << "--------------------------------------------" << std::endl;
  // //   std::cout << std::endl;
  //   
  //   // draw a map
  //   
  //   // MPPC labels
  //   
  //   
  for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
  {
    //     std::cout << "|---"
    for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
    {
      
      std::cout << mppc[iMppc][jMppc]->GetLabel() << "\t";
    }  
    std::cout << std::endl;
  }
  
  for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
  {
    for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
    {
      std::cout << mppc[iMppc][jMppc]->GetDigitizerChannel() << "\t";
    }  
    std::cout << std::endl;
  }
  
  for(int iCrystal = 0; iCrystal < ncrystalsx*nmppcx*nmodulex ; iCrystal++)
  {
    for(int jCrystal = 0; jCrystal < ncrystalsy*nmppcy*nmoduley ; jCrystal++)
    {
      std::cout << crystal[iCrystal][jCrystal]->GetID() << "\t";
    } 
    std::cout << std::endl;
  }
  
  
  
}










