#include <iostream>

#include "Crystal.h"

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



void Crystal::PrintSpecific()
{
  
//   std::cout << "Element Type \t = crystal"  << std::endl;
  
}