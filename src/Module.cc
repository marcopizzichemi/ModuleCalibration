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

void Module::PrintSpecific()
{
  std::cout << "Element Type \t = module"  << std::endl;
}

Mppc* Module::GetMppc()//TODO
{
  
}

void Module::MakeMppcPointers(int i , int j)
{
  iChildren = i;   // this sets the number of iChildren in this element
  jChildren = j;   // this sets the number of jChildren in this element
  mppc = new Mppc** [i];
  for(int k = 0 ; k < i ; k++)
  {
    mppc[k] = new Mppc* ;
  }
}
