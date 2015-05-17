#include <iostream>

#include "Crystal.h"

Crystal::Crystal()
{
  //default contructor
  name = "Default Crystal";
}

Crystal::Crystal(const Crystal &obj) 
{
  //copy contructor
}

Crystal::~Crystal()
{
  //destructor
}

void Crystal::PrintSpecific()
{
  
  std::cout << "Element Type \t = crystal"  << std::endl;
  
}