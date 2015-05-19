#include <iostream>

#include "Mppc.h"

Mppc::Mppc()
{
  //default contructor
  name = "Default MPPC";
//   parentModule = NULL;
  digitizerChannel = -1;
//   label = "VOID";
  canvasPosition = -1;
//   RawSpectrum = new TH1F("Spectrum","Spectrum",);
}

// Mppc::Mppc(int nCrystalPerMppcX , int nCrystalPerMppcY)
// {
//   iChildren = nCrystalPerMppcX;
//   jChildren = nCrystalPerMppcY;
//   crystal = new Crystal**[nCrystalPerMppcX];
//   for(int i = 0 ; i < nCrystalPerMppcX; i++)
//   {
//     crystal[i] = new Crystal*[nCrystalPerMppcY];
//   }
// }

Mppc::Mppc(const Mppc &obj) 
{
  std::cout << "Copy ctor"  << std::endl;
//   parentModule = new Element;
//   *parentModule = *obj.parentModule;
//   
//   crystal = new Element**;
//   ***crystal = ***obj.crystal;
  
//   RawSpectrum = new TH1F();
//   *RawSpectrum = *obj.RawSpectrum;
  
}

Mppc::~Mppc()
{
  std::cout << "destructing mppc " << name << std::endl;
//   if(parentModule) delete parentModule;
//   if(crystal) delete crystal;
//   if(RawSpectrum) delete RawSpectrum;
}

void Mppc::SetModule(Module *amodule)
{
//   std::cout << "--------" << std::endl;
//   std::cout << amodule << std::endl;
  parentModule = new Element;
//   std::cout << parentModule << std::endl;
  *parentModule = *((Element*)amodule);
//   std::cout << parentModule << std::endl;
//   std::cout << "--------" << std::endl;
}

// void Mppc::SetCrystal(int i, int j, Element *pCrystal)
// {
//   crystal[i][j] = new Element;
//   *crystal[i][j] = *pCrystal;
// }

void Mppc::SetCrystal(Crystal *pCrystal)
{
  Element *aElement = new Element;
  *aElement = *((Element*)pCrystal);
  vCrystal.push_back(aElement);
}

Crystal* Mppc::GetCrystal(int pi, int pj)
{
  return (Crystal*)vCrystal[pi*jChildren + pj];
}


// void Mppc::MakeCrystalPointers(int i , int j)
// {
//   iChildren = i;   // this sets the number of iChildren in this element
//   jChildren = j;   // this sets the number of jChildren in this element
//   crystal = new Element** [i];
//   for(int k = 0 ; k < i ; k++)
//   {
//     crystal[k] = new Element* ;
//   }
//   
// }

void Mppc::PrintSpecific()
{
//   std::cout << "Element Type \t = mppc"  << std::endl;
//   std::cout << "Digitizer Ch \t = "  << digitizerChannel <<  std::endl;
//   std::cout << "Label \t\t = " << label  << std::endl;
//   std::cout << "Canvas Posit. \t = "<< canvasPosition  << std::endl;
}