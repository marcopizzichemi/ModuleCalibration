#ifndef INPUTFILE_H
#define INPUTFILE_H

#include "TChain.h"
#include "TTree.h"
#include <string>
#include "ConfigFile.h"
#include "Element.h"
#include "Crystal.h"
#include "Module.h"
#include "Mppc.h"

class InputFile
{
  
  
private:

  TChain*                        fchain;
  std::string                    fname;
  TTree*                         ftree;
  int                            inputChannels;
  int                            adcChannels;
  bool                           DigitizerChannelOn[32];
  
  //variables read from the config file 
  std::string                    ConfigFileName;
  std::string                    chainName; 
  std::string                    digitizer_s;    
  std::string                    mppc_s;         
  std::string                    plotPositions_s;
  std::string                    xPositions_s;   
  std::string                    yPositions_s;  
  std::string                    saturation_s;
  std::vector <std::string>      digitizer_f;
  std::vector <std::string>      mppc_f;
  std::vector <std::string>      plotPositions_f;
  std::vector <std::string>      xPositions_f;
  std::vector <std::string>      yPositions_f;
  std::vector <std::string>      saturation_f;
  std::vector <int>              digitizer;
  std::vector <std::string>      mppc_label;
  std::vector <int>              plotPositions;
  std::vector <float>            xPositions;
  std::vector <float>            yPositions;
  std::vector <float>            saturation;
  
  int                            ncrystalsx;
  int                            ncrystalsy;
  int                            nmppcx;
  int                            nmppcy;
  int                            nmodulex;
  int                            nmoduley;
  std::string                    BinaryOutputFileName;
  bool                           binary;
  bool                           correctingSaturation;
  
  
  
  
  
  //variables for the input TChain
  ULong64_t     ChainExtendedTimeTag;
  ULong64_t     ChainDeltaTimeTag;
  Short_t      *ChainAdcChannel;
  //branches for the input TChain
  TBranch      *bChainExtendedTimeTag;
  TBranch      *bChainDeltaTimeTag;
  TBranch     **bChainAdcChannel;
  //variables for the analysis TTree
  ULong64_t     TreeExtendedTimeTag;
  ULong64_t     TreeDeltaTimeTag;
  Short_t      *TreeAdcChannel;
  int           TreeTriggerChannel;
  Float_t       TreeFloodX;
  Float_t       TreeFloodY;
  Float_t       TreeFloodZ;
  Float_t       TreeTheta;
  Float_t       TreePhi;
  Bool_t        TreeBadevent;
  
  
public:
  
  InputFile     (int argc, char** argv, ConfigFile& config);
  TChain*       GetChain() const { return fchain; };
  TTree*        GetTree() const { return ftree; };
  void          CreateTree();
  void          FillElements(Module*** module,Mppc*** mppc,Crystal*** crystal);
  
};




#endif  // INPUTFILE_H