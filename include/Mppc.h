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
  Element*               parentModule;           // one pointer for its parent, it's only one
  std::vector<Element*>  vCrystal;               // vector of pointers for the children crystals
  int                    digitizerChannel;       //which digitizer channel is assigned to this mppc
  int                    canvasPosition;         // position in the canvas of all channels
  std::string            moduleName;             // name of the module
  
  //histograms
  TH1F                   RawSpectrum;            // raw spectrum of all events seen by this mppc
  TH1F                   TriggerSpectrum;        // raw spectrum of all events seen by this mppc
  
  
  
public:
  Mppc();                                        // default constructor
  Mppc(const Mppc &obj);                         // copy constructor
  ~Mppc();                                       // destructor
  
  // methods to get and set the private variables. Names should be self explanatory
  void                   SetModule(Module *amodule); 
  Module*                GetModule(){return (Module*)parentModule;};
  void                   SetDigitizerChannel(int num){digitizerChannel = num;};
  int                    GetDigitizerChannel(){return digitizerChannel;};
  void                   SetCanvasPosition(int num){canvasPosition = num;};
  int                    GetCanvasPosition(){return canvasPosition;};
  void                   SetCrystal(Crystal *pCrystal);
  int                    GetCrystalsNumber(){return vCrystal.size();};
  Crystal*               GetCrystal(int pi, int pj);
  TH1F*                  GetRawSpectrum(){return &RawSpectrum;};
  void                   SetRawSpectrum(TH1F aHisto){RawSpectrum = aHisto;};
  TH1F*                  GetTriggerSpectrum(){return &TriggerSpectrum;};
  void                   SetTriggerSpectrum(TH1F aHisto){TriggerSpectrum = aHisto;};
  
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // MPPC_H