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
  TH2F                 VersusTime;           // 2d histogram to plot the evolution of photopeak with time (in case of gain drift?)
  
  TCut                 Ellipses;             // the elliptical TCut
  bool                 isOn;                 // if the crystal is on/off
  TEllipse             GraphicalCut;         // TEllipse to visualize the cut on the u,v global plot
  float                peakPosition;         // position of mean (after fitting) for the photopeak 
  float                peakSigma;            // sigma (after fitting) for the photopeak
  TF1                  Fit;                  // fit function (it's a gaussian)
  double               w_fwhm;               // width at half maximum for the w histogram
  double               w_rms;                // rms of w histogram
  double               w_width20perc;        // width at 20% maximum for the w histogram
  
  
public:
  Crystal();                                 // default constructor
  Crystal(const Crystal &obj);               // copy constructor
  ~Crystal();                                // destructor
  
  // methods to get and set the private variables. Names should be self explanatory
  Mppc*                GetMppc(){return (Mppc *)parentMppc;};
  TH1F*                GetSpectrum(){return &Spectrum;};
  TH1F*                GetHighlightedSpectrum(){return &HighlightedSpectrum;};
  TH1F*                GetHistoW(){return &HistoW;};
  TF1*                 GetFit(){return &Fit;};
  double               GetWfwhm(){return w_fwhm;};
  double               GetWrms(){return w_rms;};
  double               GetWwidth20perc(){return w_width20perc;};
  TCut                 GetCrystalCut(){return Ellipses;};
  TEllipse*            GetGraphicalCut(){return &GraphicalCut;};
  float                GetPhotopeakPosition(){return peakPosition;};
  float                GetPhotopeakSigma(){return peakSigma;};
  float                GetPhotopeakEnergyResolution(){return ((peakSigma*2.35)/peakPosition);};
  bool                 CrystalIsOn(){return isOn;};
  TH2F*                GetVersusTime(){return &VersusTime;};
  
  
  void                 SetMppc(Mppc *amppc);
  void                 SetSpectrum(TH1F aHisto){Spectrum = aHisto;};
  void                 SetHighlightedSpectrum(TH1F aHisto){HighlightedSpectrum = aHisto;};
  void                 SetHistoW(TH1F aHisto){HistoW = aHisto;};
  void                 SetFit(TF1 aFit){Fit = aFit;};
  void                 SetHistoWfwhm(double a){w_fwhm = a;};
  void                 SetHistoWrms(double a){w_rms = a;};
  void                 SetHistoWwidth20perc(double a){w_width20perc = a;};
  void                 SetEllipses(double u,double v,double a,double b,double t);
  void                 SetCrystalOn(bool abool){isOn = abool;};
  void                 SetGraphicalCut(TEllipse aEllipse){GraphicalCut = aEllipse;};
  void                 SetPhotopeak(float a, float b){peakPosition = a;peakSigma = b;};
  void                 SetVersusTime(TH2F aHisto){VersusTime = aHisto;};
  
  void PrintGlobal();
  void PrintSpecific();
};

#endif  // CRYSTAL_H