#include <iostream>
#include "InputFile.h"
#include "TChain.h"


InputFile::InputFile (int argc, char** argv, std::string chainName): fname(chainName)
{
  // create braches for reading the input files
  
  fchain = new TChain(fname.c_str());
  
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
  
  
  
}