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
  if(parentMppc != NULL)
  {
    std::cout << "Parent Mppc \t = "  << parentMppc->GetI() << "." << parentMppc->GetJ() <<  std::endl;
    if(parentMppc->GetModule() != NULL)
    {
      std::cout << "Parent Module \t = " << parentMppc->GetModule()->GetI() << "." << parentMppc->GetModule()->GetJ() <<  std::endl;
    }
    else
    {
      std::cout << "Crystal not assigned to any Module."  <<  std::endl;
    }
  }
  else
  {
    std::cout << "Crystal not assigned to any Mppc."  <<  std::endl;
  }
  
}