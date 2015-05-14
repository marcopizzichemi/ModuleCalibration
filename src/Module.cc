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

Element* Module::GetMppc(int digiCh) //get mppc from the digitizer channel that has been set
{
  for(int i = 0 ; i < iChildren ; i ++)
  {
    for(int j = 0 ; j < jChildren ; j++)
    {
      if(digiCh == GetChild(iChildren,jChildren)->GetDigitizerChannel())
      {
	return GetChild(iChildren,jChildren);
      }
    }
  }
  
}

Element* Module::GetMppc(std::string alabel) //get mppc from the digitizer channel that has been set
{
  for(int i = 0 ; i < iChildren ; i ++)
  {
    for(int j = 0 ; j < jChildren ; j++)
    {
      if(/*digiCh == GetChild(iChildren,jChildren)->GetDigitizerChannel()*/)
      {
	return GetChild(iChildren,jChildren);
      }
    }
  }
  
}