#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <string>
#include <iostream>

#include "Element.h"

class Crystal : public Element
{
private:
  Element* parentMppc = NULL;
  //spectra and co.
  TH1F                 Spectrum;
  TH1F                 HistoW;
//   int         mppcID;
  
public:
  Crystal(); // default constructor
  //Crystal(std::string aname, int pid, float px, float py, float pz); // constructor
  Crystal(const Crystal &obj); // copy constructor
  ~Crystal(); // destructor
  
  Mppc*                     GetMppc(){return (Mppc *)parentMppc;};
  void                      SetMppc(Mppc *amppc);
  
  
  TH1F*                GetSpectrum(){return &Spectrum;};
  TH1F*                GetHistoW(){return &HistoW;};
  
  void                 SetSpectrum(TH1F aHisto){Spectrum = aHisto;};
  void                 SetHistoW(TH1F aHisto){HistoW = aHisto;};
  
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // CRYSTAL_H