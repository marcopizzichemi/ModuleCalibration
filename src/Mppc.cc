#include <iostream>

#include "Mppc.h"

Mppc::Mppc()
{
  //default contructor
  name = "Default MPPC";
  
}
Mppc::Mppc(const Mppc &obj) 
{
  //copy contructor
}

Mppc::~Mppc()
{
  //destructor
}

// Mppc::MakeCrystalsPointers(int i, int j)
// {
//   childrenCrystal = new Crystal** [i];
//   for(int k = 0 ; k < i ; k++)
//   {
//     childrenCrystal[k] = new Crystal* [j]
//   }
// }

void Mppc::PrintSpecific()
{
  std::cout << "Element Type \t = mppc"  << std::endl;
}