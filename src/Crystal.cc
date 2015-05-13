#include <iostream>

#include "Crystal.h"

Crystal::Crystal()
{
  //default contructor
  name = "Default Name";
  x = y = z = 0;
}

// Crystal::Crystal(std::string aname, int pid, float px, float py, float pz)
// {
//   //contructor
//   name = aname;
//   id = pid;
//   x = px;
//   y = py;
//   z = pz;
//   
//   //BEGIN of debug output
//   std::cout << "Constructed Crystal " << name << ", ID " << id << ", in " << x << "," << y << "," << z << std::endl;
//   //END of debug output
// }

Crystal::Crystal(const Crystal &obj) 
{
  //copy contructor
}

Crystal::~Crystal()
{
  //destructor
}
