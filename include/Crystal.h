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
  Element* parentMppc;                       // pointer for parent element
  //spectra and co.
  TH1F                 Spectrum;             // charge spectrum for this crystal. It's always the sum of all mppcs charges
  TH1F                 HighlightedSpectrum;  // same spectrum above, but in green and only for the photopeak
  TH1F                 HistoW;               // histogram of w values for this crystal
  
  TCut                 Ellipses;             // the elliptical TCut
  bool                 isOn;                 // if the crystal is on/off
  TEllipse             GraphicalCut;         // TEllipse to visualize the cut on the u,v global plot
  float                peakPosition;         // position of mean (after fitting) for the photopeak 
  float                peakSigma;            // sigma (after fitting) for the photopeak
  TF1                  Fit;                  // fit function (it's a gaussian)
  double               w_fwhm;               // width at half maximum for the w histogram
  
  
public:
  Crystal();                                 // default constructor
  Crystal(const Crystal &obj);               // copy constructor
  ~Crystal();                                // destructor
  
  // methods to get and set the private variables. Names should be self explanatory
  Mppc*                GetMppc(){return (Mppc *)parentMppc;};
  void                 SetMppc(Mppc *amppc);
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