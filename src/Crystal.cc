#include <iostream>
#include "Crystal.h"
#include "TCut.h"

Crystal::Crystal()
{
  //default contructor
  name = "Default Crystal";
  averageDoiResolution = 0;
}

Crystal::Crystal(const Crystal &obj)
{
  std::cout << "Copy ctor"  << std::endl;
}

Crystal::~Crystal()
{
  std::cout << "destructing crystal " << name << std::endl;
  //   if(parentMppc) delete parentMppc;
}

void Crystal::SetMppc(Mppc *amppc)
{
  parentMppc = new Element;
  *parentMppc = *((Element *)amppc);
}

float Crystal::GetPhotopeakEnergyResolution()
{
  if(peakPosition != 0) return ((peakSigma*2.355)/peakPosition);
  else return 0;
}

float Crystal::GetPhotopeakPositionCorrected()
{
  return peakPositionCorrected;
}

float Crystal::GetPhotopeakSigmaCorrected()
{
  return peakSigmaCorrected;
}

float Crystal::GetPhotopeakEnergyResolutionCorrected()
{
  if(peakPositionCorrected != 0) return ((peakSigmaCorrected*2.355)/peakPositionCorrected);
  else return 0;
}


void Crystal::PrintSpecific()
{
  std::cout << "TCut \t \t= " << Ellipses  << std::endl;
  std::cout << "On \t\t= " << isOn  << std::endl;
}
