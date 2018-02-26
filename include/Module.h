#ifndef MODULE_H
#define MODULE_H

#include <iostream>
#include <string>
#include <vector>
#include "Element.h"
// #include "Detector.h"

class Module : public Element
{
private:
  std::vector<Element*>  vMppc;
  std::vector<int> channels;
  UInt_t seed;
  std::vector<detector_t> detector;

public:
  Module(); // default constructor
  Module(const Module &obj); // copy constructor
  ~Module(); // destructor

  void                   SetMppc(Mppc *pMppc);
  void                   SetDetector(std::vector<detector_t> aDetector){detector = aDetector;};
  std::vector<detector_t>   GetDetector(){return detector;};
  void                   SetChannels(std::vector<int> aVec){channels = aVec;};
  std::vector<int>       GetChannels(){return channels;};
  int                    GetMppcsNumber(){return vMppc.size();};
  Mppc*                  GetMppc(int pi, int pj);
  void                   SetSeed(UInt_t aSeed){seed = aSeed;};
  UInt_t                 GetSeed(){return seed;};

  void PrintGlobal();
  void PrintSpecific();
};

#endif  // MODULE_H
