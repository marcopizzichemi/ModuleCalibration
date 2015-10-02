#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <string>
#include <iostream>

#include "Element.h"
#include "TCut.h"
#include "TEllipse.h"
#include "TF1.h"

class Crystal : public Element
{
private:
  Element* parentMppc = NULL;
  //spectra and co.
  TH1F                 Spectrum;
  TH1F                 HighlightedSpectrum;
  TH1F                 HistoW;
  
  TCut                 Ellipses;
  bool                 isOn;
  TEllipse             GraphicalCut;
  float                peakPosition;
  float                peakSigma;
  TF1                  Fit;
  double               w_fwhm;
//   int         mppcID;
  
  
public:
  Crystal(); // default constructor
  //Crystal(std::string aname, int pid, float px, float py, float pz); // constructor
  Crystal(const Crystal &obj); // copy constructor
  ~Crystal(); // destructor
  
  Mppc*                     GetMppc(){return (Mppc *)parentMppc;};
  void                      SetMppc(Mppc *amppc);
  
  
  TH1F*                GetSpectrum(){return &Spectrum;};
  TH1F*                GetHighlightedSpectrum(){return &HighlightedSpectrum;};
  TH1F*                GetHistoW(){return &HistoW;};
  TF1*                 GetFit(){return &Fit;};
  double               GetWfwhm(){return w_fwhm;};
  
  void                 SetSpectrum(TH1F aHisto){Spectrum = aHisto;};
  void                 SetHighlightedSpectrum(TH1F aHisto){HighlightedSpectrum = aHisto;};
  void                 SetHistoW(TH1F aHisto){HistoW = aHisto;};
  void                 SetFit(TF1 aFit){Fit = aFit;};
  void                 SetHistoWfwhm(double a){w_fwhm = a;};
  
  void                 SetEllipses(double u,double v,double a,double b,double t);
  TCut                 GetCrystalCut(){return Ellipses;};
  
  void                 SetGraphicalCut(TEllipse aEllipse){GraphicalCut = aEllipse;};
  TEllipse*            GetGraphicalCut(){return &GraphicalCut;};
  
  void                 SetCrystalOn(bool abool){isOn = abool;};
  bool                 CrystalIsOn(){return isOn;};
  
  void                 SetPhotopeak(float a, float b){peakPosition = a;peakSigma = b;};
  float                GetPhotopeakPosition(){return peakPosition;};
  float                GetPhotopeakSigma(){return peakSigma;};
  float                GetPhotopeakEnergyResolution(){return ((peakSigma*2.35)/peakPosition);};
  
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // CRYSTAL_H