#include <iostream>

#include "Module.h"

Module::Module()
{
  //default contructor
  name = "Default Module Name";
  x = y = z = 0;
  //BEGIN of debug output
  std::cout << "Constructed Default Module" << std::endl;
  //END of debug output
}

Module::Module(std::string aname, float px, float py, float pz)
{
  //contructor
  name = aname;
  x = px;
  y = py;
  z = pz;
  
  //BEGIN of debug output
  std::cout << "Constructed Module " << name << " in " << x << "," << y << "," << z << std::endl;
  //END of debug output
}

Module::Module(const Module &obj) 
{
  //copy contructor
}

Module::~Module()
{
  //destructor
}

