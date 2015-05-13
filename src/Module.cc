#include <iostream>

#include "Module.h"

Module::Module()
{
  //default contructor
  name = "Default Module";
}

Module::Module(const Module &obj) 
{
  //copy contructor
}

Module::~Module()
{
  //destructor
}

// Module::MakeMppcPointers(int i, int j)
// {
//   childrenMppc = new Mppc** [i];
//   for(int k = 0 ; k < i ; k++)
//   {
//     childrenMppc[k] = new Mppc* [j]
//   }
// }

void Module::PrintSpecific()
{
  std::cout << "Element Type \t = module"  << std::endl;
}