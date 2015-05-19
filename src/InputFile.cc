#include <iostream>
#include <sstream>

#include "TChain.h"

#include "InputFile.h"
#include <cmath> 

InputFile::InputFile (int argc, char** argv, std::string chainName,int nCh): fname(chainName), inputChannels(nCh)
{
  // Constructor 
  // opens the Tchain, set its branches, create the TTree that will be used for the analysis
  
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
  ftree->Branch("TriggerChannel",&TreeTriggerChannel,"TriggerChannel/I"); 
  ftree->Branch("FloodX",&TreeFloodX,"FloodX/F"); 
  ftree->Branch("FloodY",&TreeFloodY,"FloodY/F"); 
  ftree->Branch("FloodZ",&TreeFloodZ,"FloodZ/F");
  ftree->Branch("Theta",&TreeTheta,"Theta/F"); 
  ftree->Branch("Phi",&TreePhi,"Phi/F");
  ftree->Branch("BadEvent",&TreeBadevent,"BadEvent/O"); 
  
}

void InputFile::CreateTree(std::vector<int> digitizer , std::vector<float> xmppc , std::vector<float> ymppc)
{
  //creates the TTree from the input Tchain
  std::cout << "Filling the TTree for the analysis... " << std::endl;
  Int_t nevent = fchain->GetEntries();
  long long int GoodCounter = 0;
  long long int BadEvent = 0;
  long long int counter = 0;
  
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
    int badCharge = 0;
    bool badevent = false;
    
    TreeExtendedTimeTag = ChainExtendedTimeTag;
    TreeDeltaTimeTag = ChainDeltaTimeTag;
    
    int TreeEntryCounter = 0;
    
    
    for (int j = 0 ; j < 32 ; j++) //loop on input chain channels
    {
      if(DigitizerChannelOn[j])
      {
//       for(int k= 0 ; k < digitizer.size() ; k++) // loop on the channels of the digitizer for this module
//       {
// 	if(k == j) //so if this channel was enabled for input of this module, save the data in the analysis ttree
// 	{
	  //TODO implement the saturation correction here
	  
	  // fill tree with data from the channels
	  TreeAdcChannel[TreeEntryCounter] = ChainAdcChannel[j];
	  	  
	  //find the max charge and therefore the TriggerChannel
	  if (TreeAdcChannel[TreeEntryCounter] > maxCharge)
	  {
	    maxCharge = TreeAdcChannel[TreeEntryCounter];
	    TreeTriggerChannel = TreeEntryCounter;
	  }
	  
	  //update total, row and columnsum
	  total += TreeAdcChannel[TreeEntryCounter];
	  rowsum += TreeAdcChannel[TreeEntryCounter]*xmppc[TreeEntryCounter];
	  columnsum += TreeAdcChannel[TreeEntryCounter]*ymppc[TreeEntryCounter];
	  
	  
	  TreeEntryCounter++;
// 	}
//       }
      }
    }
    
    TreeFloodX = rowsum/total;
    TreeFloodY = columnsum/total;
    TreeFloodZ =  maxCharge/total;
    
    TreeTheta = std::acos(TreeFloodZ /( std::sqrt( std::pow(TreeFloodX-0,2) + std::pow(TreeFloodY-0,2) + std::pow(TreeFloodZ+0,2)) ));
    TreePhi =  std::atan (TreeFloodY / TreeFloodX);
    
    TreeBadevent = false;
    ftree->Fill();
    
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
  std::cout << std::endl;
  
  std::cout << "Tot events = \t" << counter << std::endl;
//   std::cout << "Accepted events = \t" << GoodCounter << std::endl;
  //std::cout << "Bad events = \t" << badEvents << std::endl;
  
}