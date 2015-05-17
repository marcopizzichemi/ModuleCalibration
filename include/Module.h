#ifndef MODULE_H
#define MODULE_H

#include <string>
#include <vector>
#include "Element.h"
#include "TH2F.h"


class Module : public Element
{
private: 
  Mppc*** mppc; // matrix of pointers to children mppcs
//   std::vector<Mppc*>  mppc; 
  TH2F* FloodMap2D;
  
public:
  Module(); // default constructor
  Module(const Module &obj); // copy constructor
  ~Module(); // destructor
    
  void MakeMppcPointers(int i , int j);
  Mppc* GetMppc();
  void SetMppc(int i, int j, Mppc* pMppc){mppc[i][j] = pMppc;};
  
  void SetFloodMap2D(TH2F* aHisto){FloodMap2D = aHisto;};
  TH2F* GetFloodMap2D(){return FloodMap2D;};
  
  void PrintGlobal();
  void PrintSpecific();
};

#endif  // MODULE_H