#ifndef MODULE_H
#define MODULE_H

#include <iostream>
#include <string>
#include <vector>
#include "Element.h"


class Module : public Element
{
private:
  std::vector<Element*>  vMppc;
  UInt_t seed;


public:
  Module(); // default constructor
  Module(const Module &obj); // copy constructor
  ~Module(); // destructor

  void                   SetMppc(Mppc *pMppc);
  int                    GetMppcsNumber(){return vMppc.size();};
  Mppc*                  GetMppc(int pi, int pj);
  void                   SetSeed(UInt_t aSeed){seed = aSeed;};
  UInt_t                 GetSeed(){return seed;};

  void PrintGlobal();
  void PrintSpecific();
};

#endif  // MODULE_H
