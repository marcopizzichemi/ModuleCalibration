#include <iostream>

#include "Crystal.h"
#include "TCut.h"

Crystal::Crystal()
{
  //default contructor
  name = "Default Crystal";
}

Crystal::Crystal(const Crystal &obj) 
{
  std::cout << "Copy ctor"  << std::endl;
//   if(parentMppc)
//   {
//     parentMppc = new Element;
//     *parentMppc = *obj.parentMppc;
//   }
//   else
//   {
//     parentMppc = NULL;
//   }
//   parentMppc = new Mppc();
//   *parentMppc = *obj.parentMppc;
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


void Crystal::SetEllipses(double u,double v,double a,double b,double t)
{
  //takes the ellipses center (u,v), width (a,b) and inclination (t) 
  //and writes the TCut expression used for selecting the events assigned to the crystal
  //units are in the u,v scale, t in degrees
  
  //create the string for the elliptic cut around each crystal
  //FIXME for the moment, the cuts are input by config file
  //convoluted way using stringstream, but TCutG and friends don't seem to work for me...
  std::stringstream crystalFloodCut;
  crystalFloodCut << "TMath::Power(((FloodX - " 
                  << u 
                  << ")*TMath::Cos( ("
		  << t
		  << "/180.0) * TMath::Pi() ) + (FloodY -"
		  << v
		  << ")*TMath::Sin( ("
		  << t
		  << "/180.0) * TMath::Pi() ) ) / "
		  << a
		  << ",2) + "
		  << "TMath::Power(((FloodX - " 
                  << u 
                  << ")*TMath::Sin( ("
		  << t
		  << "/180.0) * TMath::Pi() ) - (FloodY -"
		  << v
		  << ")*TMath::Cos( ("
		  << t
		  << "/180.0) * TMath::Pi() ) ) / "
		  << b
		  << ",2) < 1";
   Ellipses = crystalFloodCut.str().c_str();
};



void Crystal::PrintSpecific()
{
  
  std::cout << "TCut \t \t= " << Ellipses  << std::endl;
  std::cout << "On \t \t= " << isOn  << std::endl;
  
}