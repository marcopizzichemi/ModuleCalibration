#include <iostream>
#include "Element.h"

Element::Element()
{
  //default constructor
  name = "Default Element";
  id = i = j = 0;
  x = y = z = 0;
//   GlobalTag = "x.x-x.x-x.x";
//   parentID = 0;
  
}

Element::Element(const Element &obj) 
{
  //copy constructor
}

void Element::PrintGlobal()
{
  std::cout << "Name \t\t = "    << name << std::endl;
  std::cout << "ID \t\t = "      << id << std::endl;
  std::cout << "i = \t\t = "      << i << std::endl;
  std::cout << "j = \t\t = "      << j << std::endl;
  std::cout << "Position \t = (" << x << "," << y << "," << z << ")" << std::endl;
//   std::cout << "GlobalTag \t = " << GlobalTag << std::endl;
  std::cout << "Parent \t\t = " << parentName << std::endl;
  for(int k = 0; k < childrenName.size() ; k++)
  {
    std::cout << "Child \t\t = " << childrenName[k] << std::endl;
  }
}

void Element::PrintSpecific()
{
//   std::cout << "Abstract Element, nothing specific to print" << std::endl;
}