#include <iostream>

#include "Mppc.h"

Mppc::Mppc()
{
  //default contructor
  name = "Default MPPC";
  parentModule = NULL;
  digitizerChannel = -1;
  label = "VOID";
  canvasPosition = -1;
  RawSpectrum = NULL;
}
Mppc::Mppc(const Mppc &obj) 
{
  //copy contructor
}

Mppc::~Mppc()
{
  //destructor
}


void Mppc::PrintSpecific()
{
  std::cout << "Element Type \t = mppc"  << std::endl;
  std::cout << "Digitizer Ch \t = "  << digitizerChannel <<  std::endl;
  std::cout << "Label \t\t = " << label  << std::endl;
  std::cout << "Canvas Posit. \t = "<< canvasPosition  << std::endl;
}