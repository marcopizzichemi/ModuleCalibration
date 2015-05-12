#include <iostream>

#include "Mppc.h"

Mppc::Mppc()
{
  //default contructor
  name = "Default MPPC Name";
  x = y = z = 0;
  //BEGIN of debug output
  std::cout << "Constructed Default MPPC" << std::endl;
  //END of debug output
}

Mppc::Mppc(std::string aname, float px, float py, float pz)
{
  //contructor
  name = aname;
  x = px;
  y = py;
  z = pz;
  
  //BEGIN of debug output
  std::cout << "Constructed MPPC " << name << " in " << x << "," << y << "," << z << std::endl;
  //END of debug output
}

Mppc::Mppc(const Mppc &obj) 
{
  //copy contructor
}

Mppc::~Mppc()
{
  //destructor
}
