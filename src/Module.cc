#include <iostream>

#include "Module.h"

Module::Module()
{
  //default contructor
  name = "Default Module Name";
  id = -1;
  x = y = z = 0;
}

// Module::Module(std::string aname, float px, float py, float pz)
// {
//   //contructor
//   name = aname;
//   x = px;
//   y = py;
//   z = pz;
//   
//   //BEGIN of debug output
//   std::cout << "Constructed Module " << name << " in " << x << "," << y << "," << z << std::endl;
//   //END of debug output
// }

Module::Module(const Module &obj) 
{
  //copy contructor
}

Module::~Module()
{
  //destructor
}

void Module::PrintSpecific()
{
  std::cout << "Element Type \t = module"  << std::endl;
}