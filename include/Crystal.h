#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <string>

#include "Element.h"

class Crystal : public Element
{
// protected:             //inherited from Element
//   std::string name;  // module name or label
//   float x,y,z;       // module position in "world" coordinate space
//   
public:
  Crystal(); // default constructor
  //Crystal(std::string aname, int pid, float px, float py, float pz); // constructor
  Crystal(const Crystal &obj); // copy constructor
  ~Crystal(); // destructor
    
  
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // CRYSTAL_H