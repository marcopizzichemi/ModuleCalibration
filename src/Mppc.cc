#include <iostream>

#include "Mppc.h"

Mppc::Mppc()
{
  //default contructor
  name = "Default MPPC";
  parentModule = NULL;
  digitizerChannel = -1;
//   label = "VOID";
  canvasPosition = -1;
  RawSpectrum = NULL;
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
  //copy contructor
}

Mppc::~Mppc()
{
  //destructor
}

Crystal* Mppc::GetCrystal() //TODO
{
  
}

void Mppc::MakeCrystalPointers(int i , int j)
{
  iChildren = i;   // this sets the number of iChildren in this element
  jChildren = j;   // this sets the number of jChildren in this element
  crystal = new Crystal** [i];
  for(int k = 0 ; k < i ; k++)
  {
    crystal[k] = new Crystal* ;
  }
  
}

void Mppc::PrintSpecific()
{
  std::cout << "Element Type \t = mppc"  << std::endl;
  std::cout << "Digitizer Ch \t = "  << digitizerChannel <<  std::endl;
//   std::cout << "Label \t\t = " << label  << std::endl;
  std::cout << "Canvas Posit. \t = "<< canvasPosition  << std::endl;
}