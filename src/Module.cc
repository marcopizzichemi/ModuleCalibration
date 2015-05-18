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
//   std::cout << "copy constructor" << std::endl;
//   mppc = new Element**;
//   
//   ***mppc = ***obj.mppc;
//   std::cout << "qui" << std::endl;
  
  /*
  FloodMap2D = new TH2F();
  *FloodMap2D = *obj.FloodMap2D;*/
}

Module::~Module()
{
   std::cout << "destructing module " << name << std::endl;
//   std::cout << "Destructing " << name << std::endl;
//   if(mppc) delete mppc;
//   if(FloodMap2D) delete FloodMap2D;
}

void Module::PrintSpecific()
{
//   std::cout << "Element Type \t = module"  << std::endl;
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

// void Module::SetMppc(int i, int j, Element* pMppc)
// {
//   mppc[i][j] = new Element;
// //   std::cout << pMppc << std::endl;
// //   std::cout << mppc[i][j] << std::endl;
//   *mppc[i][j] = *pMppc;
// //   std::cout << mppc[i][j] << std::endl;
// //   std::cout << std::endl;
// }


// Mppc* Module::GetMppc()//TODO
// {
//   
// }

// void Module::MakeMppcPointers(int i , int j)
// {
//   iChildren = i;   // this sets the number of iChildren in this element
//   jChildren = j;   // this sets the number of jChildren in this element
//   mppc = new Element** [i];
//   for(int k = 0 ; k < i ; k++)
//   {
//     mppc[k] = new Element* ;
//   }
// }
