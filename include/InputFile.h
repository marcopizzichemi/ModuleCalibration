#ifndef INPUTFILE_H
#define INPUTFILE_H

#include "TChain.h"
#include "TTree.h"
#include <string>

class InputFile
{
  
  
private:

  TChain*       fchain;
  std::string   fname;
  TTree*        ftree;
  int           inputChannels;
  
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
  Bool_t        TreeBadevent;
  
  
public:
  
  InputFile     (int argc, char** argv, std::string chainName,int nCh);
  TChain*       GetChain() const { return fchain; };
  TTree*        GetTree() const { return ftree; };
  
  
};




#endif  // INPUTFILE_H