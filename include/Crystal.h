#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <string>
#include <iostream>

#include "Element.h"
#include "TCut.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TGraph.h"

class Crystal : public Element
{
  
private:
  Element* parentMppc;                       // pointer for parent element
  //spectra and co.
  TH1F                 Spectrum;             // charge spectrum for this crystal. It's always the sum of all mppcs charges
  TH1F                 HighlightedSpectrum;  // same spectrum above, but in green and only for the photopeak
  TH1F                 HistoW;               // histogram of w values for this crystal
  TH2F                 VersusTime;           // 2d histogram to plot the evolution of photopeak with time (in case of gain drift?)
  TH2F                 SimDOIplot;           // 2d histogram for simulation, showing z versus w
  TGraph               SimGraph;             
  TCut                 Ellipses;             // the elliptical TCut
  TCut                 w20percCut;           ///< TCut from first bin above 20% to last bin above 20% for the w plot
  bool                 isOn;                 // if the crystal is on/off
  TEllipse             GraphicalCut;         // TEllipse to visualize the cut on the u,v global plot
  float                peakPosition;         // position of mean (after fitting) for the photopeak 
  float                peakSigma;            // sigma (after fitting) for the photopeak
  TF1                  Fit;                  // fit function (it's a gaussian)
  TF1                  SimFit;
  TF1                  Wfit;
  TF1                  ProfileXFit;
  double               w_fwhm;               // width at half maximum for the w histogram
  double               w_rms;                // rms of w histogram
  double               w_width20perc;        // width at 20% maximum for the w histogram
  double               u,v,wu,wv,t;
  
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
  TF1*                 GetSimFit(){return &SimFit;};
  TF1*                 GetHistoWfit(){return &Wfit;};
  TF1*                 GetProfileXFit(){return &ProfileXFit;};
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
  TH2F*                GetSimDOIplot(){return &SimDOIplot;};
  TGraph*              GetSimGraph(){return &SimGraph;};
  double               GetU(){return u;};
  double               GetV(){return v;};
  double               GetWU(){return wu;};
  double               GetWV(){return wv;};
  double               GetT(){return t;};
  TCut                 GetW20percCut(){return w20percCut;};
  
  void                 SetMppc(Mppc *amppc);
  void                 SetSpectrum(TH1F aHisto){Spectrum = aHisto;};
  void                 SetHighlightedSpectrum(TH1F aHisto){HighlightedSpectrum = aHisto;};
  void                 SetHistoW(TH1F aHisto){HistoW = aHisto;};
  void                 SetFit(TF1 aFit){Fit = aFit;};
  void                 SetHistoWfwhm(double a){w_fwhm = a;};
  void                 SetHistoWrms(double a){w_rms = a;};
  void                 SetHistoWwidth20perc(double a){w_width20perc = a;};
  void                 SetHistoWfit(TF1 aFit){Wfit = aFit;};
  void                 SetEllipses(std::string varX,std::string varY);
  void                 SetCrystalOn(bool abool){isOn = abool;};
  void                 SetCrystalData(double au,double av,double awu ,double awv, double at){u = au; v = av; wu = awu ; wv = awv ; t = at;};
  void                 SetGraphicalCut(TEllipse aEllipse){GraphicalCut = aEllipse;};
  void                 SetPhotopeak(float a, float b){peakPosition = a;peakSigma = b;};
  void                 SetVersusTime(TH2F aHisto){VersusTime = aHisto;};
  void                 SetSimDOIplot(TH2F aHisto){SimDOIplot = aHisto;};
  void                 SetSimGraph(TGraph aGraph){SimGraph = aGraph;};
  void                 SetSimFit(TF1 aFit){SimFit = aFit;};
  void                 SetW20percCut(TCut aCut){w20percCut = aCut;};
  void                 SetProfileXFit(TF1 aFit){ProfileXFit = aFit;};
  void                 Analyze();
  
  void PrintGlobal();
  void PrintSpecific();
};

#endif  // CRYSTAL_H