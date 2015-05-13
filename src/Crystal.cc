#include <iostream>

#include "Crystal.h"

Crystal::Crystal()
{
  //default contructor
  name = "Default Name";
  id = -1;
  x = y = z = 0;
}

Crystal::Crystal(const Crystal &obj) 
{
  //copy contructor
}

Crystal::~Crystal()
{
  //destructor
}


// void Crystal::Print()
// {
//   std::cout << std::endl;
//   std::cout << "--------------------------------------------------"  << std::endl;
//   std::cout << "Element Type \t = crystal"  << std::endl;
//   std::cout << "Name \t\t = "   << name << std::endl;
//   std::cout << "ID \t\t = "     << id << std::endl;
//   std::cout << "Position \t = (" << x << "," << y << "," << z << ")" << std::endl;
//   std::cout << "--------------------------------------------------"  << std::endl;
//   std::cout << std::endl;
// }

void Crystal::PrintSpecific()
{
  std::cout << "Element Type \t = crystal"  << std::endl;
}