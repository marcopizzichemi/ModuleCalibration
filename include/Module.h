#ifndef MODULE_H
#define MODULE_H

#include <string>

class Module
{
private:
  std::string name;  // module name or label
  float x,y,z;       // module position in "world" coordinate space
  
public:
  Module(); // default constructor
  Module(std::string aname, float px, float py, float pz); //constructor
  Module(const Module &obj); // copy constructor
  ~Module(); // destructor
};




#endif  // MODULE_H