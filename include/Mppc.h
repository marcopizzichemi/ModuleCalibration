#ifndef MPPC_H
#define MPPC_H
#include "Element.h"
#include "TH1F.h"

#include <vector>
#include <iostream>
// #include <string>

class Mppc : public Element
{
private:
  Element*               parentModule;          // one pointer for its parent, it's only one
  std::vector<Element*>  vCrystal;               // vector of pointers for the children crystals
//   Element***             crystal;
  int                    digitizerChannel;      //which digitizer channel is assigned to this mppc
//   std::string            label;                 // label of the mppc
  int                    canvasPosition;        // position in the canvas of all channels
  
  //histograms
  //TH1F                  RawSpectrum;           // raw spectrum of all events seen by this mppc
  std::string            moduleName;
  
  
public:
  Mppc();                                       // default constructor
//   Mppc(int nCrystalPerMppcX , int nCrystalPerMppcY);
  Mppc(const Mppc &obj);                        // copy constructor
  ~Mppc();                                      // destructor
  
  void                   SetModule(Module *amodule); 
  Module*                GetModule(){return (Module*)parentModule;};
  
  void                   SetDigitizerChannel(int num){digitizerChannel = num;};
  int                    GetDigitizerChannel(){return digitizerChannel;};
  void                   SetCanvasPosition(int num){canvasPosition = num;};
  int                    GetCanvasPosition(){return canvasPosition;};
  
  void                   SetCrystal(Crystal *pCrystal);
  int                    GetCrystalsNumber(){return vCrystal.size();};
  Crystal*               GetCrystal(int pi, int pj);
  
  
//   void                   SetLabel(std::string string){label = string;};
  

//   void AddCrystal(Crystal* pCrystal){crystal.push_back(pCrystal);};
  
//   Module*                  GetModule(){return parentModule;};
  
//   Element* GetCrystal(int iCrystal , int jCrystal){return GetChild(iCrystal,jCrystal);};
//   void MakeCrystalPointers(int iCrystal , int jCrystal);
//   Crystal* GetCrystal(int iCrystal , int jCrystal){return crystal[iCrystal][jCrystal];};
//   void SetCrystal(int i, int j, Element *pCrystal);
  
  
  
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // MPPC_H