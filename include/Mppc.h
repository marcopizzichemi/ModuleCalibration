#ifndef MPPC_H
#define MPPC_H
#include "Element.h"

// #include <string>

class Mppc : public Element
{
// private:
//   std::string name;  // mppc name or label
//   float x,y,z;       // mppc position in module coordinate space
  
public:
  Mppc(); // default constructor
  //Mppc(std::string aname, int pid, float px, float py, float pz); //constructor
  Mppc(const Mppc &obj); // copy constructor
  ~Mppc(); // destructor
  
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // MPPC_H