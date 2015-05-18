#ifndef MODULE_H
#define MODULE_H

#include <iostream>
#include <string>
#include <vector>
#include "Element.h"
#include "TH2F.h"


class Module : public Element
{
private: 
//   Element*** mppc; // matrix of pointers to children mppcs
  std::vector<Element*>  vMppc; 
  TH2F* FloodMap2D;
  
public:
  Module(); // default constructor
  Module(const Module &obj); // copy constructor
  ~Module(); // destructor
    
//   void  MakeMppcPointers(int i , int j);
  void                   SetMppc(Mppc *pMppc);
  int                    GetMppcsNumber(){return vMppc.size();};
  Mppc*                  GetMppc(int pi, int pj);
//   void  SetMppc(int i, int j, Element* pMppc);
  
//   void SetFloodMap2D(TH2F* aHisto){FloodMap2D = aHisto;};
//   TH2F* GetFloodMap2D(){return FloodMap2D;};
  
  void PrintGlobal();
  void PrintSpecific();
};

#endif  // MODULE_H