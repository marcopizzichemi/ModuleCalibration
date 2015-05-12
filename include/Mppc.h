#ifndef MPPC_H
#define MPPC_H
#include "Module.h"

// #include <string>

class Mppc: public Module
{
private:
  std::string name;  // mppc name or label
  float x,y,z;       // mppc position in module coordinate space
  
public:
  Mppc(); // default constructor
  Mppc(std::string aname, float px, float py, float pz); //constructor
  Mppc(const Mppc &obj); // copy constructor
  ~Mppc(); // destructor
};




#endif  // MPPC_H