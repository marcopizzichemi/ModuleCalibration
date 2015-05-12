#ifndef INPUTFILE_H
#define INPUTFILE_H

#include "TChain.h"
#include <string>

class InputFile
{
  
private:

  TChain*       fchain;
  std::string   fname;
  
public:
  
  InputFile     (int argc, char** argv, std::string chainName);
  TChain*       GetChain() const { return fchain; };
  
  
  
  
};

#endif  // INPUTFILE_H