#ifndef MPPC_H
#define MPPC_H
#include "Element.h"
#include "Module.h"
#include "TH1F.h"
// #include <string>

class Mppc : public Element
{
private:
  Module*                parentModule;          // one pointer for its parent, it's only one
  int                    digitizerChannel;      //which digitizer channel is assigned to this mppc
  std::string            label;                 // label of the mppc
  int                    canvasPosition;        // position in the canvas of all channels
  
  //histograms
  TH1F*                  RawSpectrum;           // raw spectrum of all events seen by this mppc
  
  
public:
  Mppc();                                       // default constructor
  Mppc(const Mppc &obj);                        // copy constructor
  ~Mppc();                                      // destructor
  
  Module*                GetModule(){return parentModule;};
  int                    GetDigitizerChannel(){return digitizerChannel;};
  std::string            GetLabel(){return label;};
  int                    GetCanvasPosition(){return canvasPosition;};
  
  void                   SetModule(Module* amodule){parentModule = amodule;}; 
  void                   SetDigitizerChannel(int num){digitizerChannel = num;};
  void                   SetLabel(std::string string){label = string;};
  void                   SetCanvasPosition(int num){canvasPosition = num;};

  
  Element* GetCrystal(int iCrystal , int jCrystal){return GetChild(iCrystal,jCrystal);};
  
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // MPPC_H