#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <string>

#include "Element.h"

class Crystal : public Element
{
private:
  Mppc* parentMppc = NULL;
//   std::string mppcLabel;
//   int         mppcID;
  
public:
  Crystal(); // default constructor
  //Crystal(std::string aname, int pid, float px, float py, float pz); // constructor
  Crystal(const Crystal &obj); // copy constructor
  ~Crystal(); // destructor
  
  Mppc*                  GetMppc(){return parentMppc;};
  void                   SetMppc(Mppc* amppc){parentMppc = amppc;}; 
  
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // CRYSTAL_H