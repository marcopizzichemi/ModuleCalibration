#ifndef MODULE_H
#define MODULE_H

#include <string>

#include "Element.h"

class Module : public Element
{
// private:
//   std::string name;  // module name or label
//   float x,y,z;       // module position in "world" coordinate space
//   
public:
  Module(); // default constructor
  //Module(std::string aname, int pid, float px, float py, float pz); //constructor
  Module(const Module &obj); // copy constructor
  ~Module(); // destructor
    
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // MODULE_H