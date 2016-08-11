#include <iostream>

#include "Module.h"

Module::Module()
{
  //default contructor
  name = "Default Module";
}

Module::Module(const Module &obj) 
{
  std::cout << "Copy ctor"  << std::endl;
}

Module::~Module()
{
  std::cout << "destructing module " << name << std::endl;
}

void Module::PrintSpecific()
{
  
}

void Module::SetMppc(Mppc *pMppc)
{
  Element *aElement = new Element;
  *aElement = *((Element*)pMppc);
  vMppc.push_back(aElement);
}

Mppc* Module::GetMppc(int pi, int pj)
{
  return (Mppc*)vMppc[pi*jChildren + pj];
}
