#include <iostream>
#include "Element.h"


Element::Element()
{
  //default constructor
  name = "Default Element";
  id = -1;
  x = y = z = 0;
  GlobalTag = "x-x.x-x.x";
  
}

// Element::Element(std::string aname, int pid, float px, float py, float pz)
// {
//   //constructor
//   name = aname;
//   id = pid;
//   x = px;
//   y = py;
//   z = pz;
// }

Element::Element(const Element &obj) 
{
  //copy constructor
}

Element::~Element()
{
  //destructor
}

void Element::SetGlobalTag(int module, int mppcx, int mppcy, int cryx , int cryy)
{
  std::stringstream stream;
  stream << module << "-" << mppcx << "." << mppcy << "-" << cryx << "." << cryy;
  GlobalTag = stream.str();
}


void Element::PrintGlobal()
{
  std::cout << "Name \t\t = "   << name << std::endl;
  std::cout << "ID \t\t = "     << id << std::endl;
  std::cout << "Position \t = (" << x << "," << y << "," << z << ")" << std::endl;
  std::cout << "GlobalTag \t = " << GlobalTag << std::endl;
}

void Element::PrintSpecific()
{
  std::cout << "Abstract Element, nothing specific to print" << std::endl;
}
