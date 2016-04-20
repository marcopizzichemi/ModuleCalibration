#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <string>
#include <iostream>

#include "Element.h"
#include "TCut.h"
#include "TCutG.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TGraph.h"

class Crystal : public Element
{
  
private:
  Element* parentMppc;                       ///< pointer for parent element
  //spectra and co.
  TH1F*                Spectrum;             ///< charge spectrum for this crystal. It's always the sum of all mppcs charges
  TH1F*                CorrectedSpectrum;    ///< charge spectrum for this crystal corrected by DOI
  TH1F*                HighlightedSpectrum;  ///< same spectrum above, but in green and only for the photopeak
  TH1F*                HighlightedSpectrumCorrected; ///<
  TH1F*                HistoWCorrectedSmooth;
  TH1F*                HistoW;               ///< histogram of w values for this crystal
  TH1F*                HistoWCorrected;
  TH1F*                DensityHisto;         ///< histogram of the entries per voxel in the crystal volume after separation
  TH1F*                pdfW;
  TH1F*                cumulativeW;
  TH1D*                SlicesMean;           ///< histogram of fitted mean values of profiles from TH2F ADCvsW distribution
  TH2F*                VersusTime;           ///< 2d histogram to plot the evolution of photopeak with time (in case of gain drift?)
  TH2F*                SimDOIplot;           ///< 2d histogram for simulation, showing z versus w
  TGraph*              SimGraph;             
  TCut                 Ellipses;             ///< the elliptical TCut
  TCut                 w20percCut;           ///< TCut from first bin above 20% to last bin above 20% for the w plot
  TCutG               *cutg[2];
  bool                 isOn;                 ///< if the crystal is on/off
  TEllipse             GraphicalCut;         ///< TEllipse to visualize the cut on the u,v global plot
  float                peakPosition;         ///< position of mean (after fitting) for the photopeak 
  float                peakSigma;            ///< sigma (after fitting) for the photopeak
  float                peakPositionCorrected;         ///< position of mean (after fitting) for the photopeak, corrected by DOI
  float                peakSigmaCorrected;            ///< sigma (after fitting) for the photopeak, corrected by DOI
  TGraph*              calibrationGraph;
  
  TF1*                 Fit;                  ///< fit function (it's a gaussian)
  TF1*                 SimFit;
  TF1*                 Wfit;
  TF1*                 FitCorrected;
  TF1*                 ThetaFit;
  TF1*                 deltaWfit;
  TF1*                 deltaWfit_2;
//   TF1                  ProfileXFit;
  TF1*                 SlicesMeanFit;
  double               w_fwhm;               ///< width at half maximum for the w histogram
  double               w_rms;                ///< rms of w histogram
  double               w_width20perc;        ///< width at 20% maximum for the w histogram
  double               u,v,wu,wv,t;
  double               wBegin;               ///< beginning of w histogram after fitting with theta function
  double               wEnd;                 ///< end of w histogram after fitting with theta function
  double               deltaW;               ///< delta of w for a fixed position, as calculated from the gaussian fit of rise in w plot
  
public:
  Crystal();                                 ///< default constructor
  Crystal(const Crystal &obj);               ///< copy constructor
  ~Crystal();                                ///< destructor
  
  // methods to get and set the private variables. Names should be self explanatory
  Mppc*                GetMppc(){return (Mppc *)parentMppc;};
  TH1F*                GetSpectrum(){return Spectrum;};
  TH1F*                GetCorrectedSpectrum(){return CorrectedSpectrum;};
  TH1F*                GetHighlightedSpectrum(){return HighlightedSpectrum;};
  TH1F*                GetHighlightedSpectrumCorrected(){return HighlightedSpectrumCorrected;};
  TH1F*                GetHistoW(){return HistoW;};
  TH1F*                GetHistoWCorrected(){return HistoWCorrected;};
  TH1F*                GetHistoWCorrectedSmooth(){return HistoWCorrectedSmooth;};
  TH1F*                GetDensityHisto(){return DensityHisto;};
  TH1F*                GetPdfW(){return pdfW;};
  TH1F*                GetCumulativeW(){return cumulativeW;};
  TGraph*              GetCalibrationGraph(){return calibrationGraph;};
  TF1*                 GetFit(){return Fit;};
  TF1*                 GetSimFit(){return SimFit;};
  TF1*                 GetHistoWfit(){return Wfit;};
//   TF1*                 GetProfileXFit(){return &ProfileXFit;};
  TF1*                 GetSlicesMeanFit(){return SlicesMeanFit;};
  double               GetWfwhm(){return w_fwhm;};
  double               GetWrms(){return w_rms;};
  double               GetWwidth20perc(){return w_width20perc;};
  TCut                 GetCrystalCut(){return Ellipses;};
  TEllipse*            GetGraphicalCut(){return &GraphicalCut;};
  float                GetPhotopeakPosition(){return peakPosition;};
  float                GetPhotopeakSigma(){return peakSigma;};
  float                GetPhotopeakEnergyResolution(){return ((peakSigma*2.355)/peakPosition);};
  float                GetPhotopeakPositionCorrected(){return peakPositionCorrected;};
  float                GetPhotopeakSigmaCorrected(){return peakSigmaCorrected;};
  float                GetPhotopeakEnergyResolutionCorrected(){return ((peakSigmaCorrected*2.355)/peakPositionCorrected);};
  double               GetWbegin(){return wBegin;};
  double               GetWend(){return wEnd;};
  double               GetDeltaW(){return std::abs(deltaW);};
   
  bool                 CrystalIsOn(){return isOn;};
  TH2F*                GetVersusTime(){return VersusTime;};
  TH2F*                GetSimDOIplot(){return SimDOIplot;};
  TGraph*              GetSimGraph(){return SimGraph;};
  double               GetU(){return u;};
  double               GetV(){return v;};
  double               GetWU(){return wu;};
  double               GetWV(){return wv;};
  double               GetT(){return t;};
  TCut                 GetW20percCut(){return w20percCut;};
  TH1D*                GetSlicesMean(){return SlicesMean;};
  TF1*                 GetFitCorrected(){return FitCorrected;};
  TCutG*               GetZXCut(){return cutg[0];};
  TCutG*               GetZYCut(){return cutg[1];};
  TF1*                 GetThetaFit(){return ThetaFit;};
  TF1*                 GetDeltaWfit(){return deltaWfit;};
  TF1*                 GetDeltaWfit_2(){return deltaWfit_2;};
  double               GetMcal(){return dz/(wBegin - wEnd);};
  double               GetQcal(){return -(dz*wEnd)/(wBegin-wEnd);};
  double               GetDoiResolutionFWHM(){return 2.355 * std::abs((dz/(wBegin - wEnd))) * std::abs(deltaW);};
  
  
  void                 SetZXCut(TCutG *aCut){cutg[0] = aCut;};
  void                 SetZYCut(TCutG *aCut){cutg[1] = aCut;};
  void                 SetMppc(Mppc *amppc);
  void                 SetSpectrum(TH1F* aHisto){Spectrum = aHisto;};
  void                 SetHighlightedSpectrum(TH1F* aHisto){HighlightedSpectrum = aHisto;};
  void                 SetHighlightedSpectrumCorrected(TH1F* aHisto){HighlightedSpectrumCorrected = aHisto;};
  void                 SetHistoW(TH1F* aHisto){HistoW = aHisto;};
  void                 SetHistoWCorrected(TH1F* aHisto){HistoWCorrected = aHisto;};
  void                 SetHistoWCorrectedSmooth(TH1F* aHisto){HistoWCorrectedSmooth = aHisto;};
  void                 SetDensityHisto(TH1F* aHisto){DensityHisto = aHisto;};
  void                 SetFit(TF1* aFit){Fit = aFit;};
  void                 SetHistoWfwhm(double a){w_fwhm = a;};
  void                 SetHistoWrms(double a){w_rms = a;};
  void                 SetHistoWwidth20perc(double a){w_width20perc = a;};
  void                 SetHistoWfit(TF1* aFit){Wfit = aFit;};
  void                 SetPdfW(TH1F* aHisto){pdfW = aHisto;};
  void                 SetCumulativeW(TH1F* aHisto){cumulativeW = aHisto;};
  void                 SetCalibrationGraph(TGraph* aGraph){calibrationGraph = aGraph;};
//   void                 SetEllipses(std::string varX,std::string varY);
  void                 SetCrystalOn(bool abool){isOn = abool;};
  void                 SetCrystalData(double au,double av,double awu ,double awv, double at){u = au; v = av; wu = awu ; wv = awv ; t = at;};
  void                 SetGraphicalCut(TEllipse aEllipse){GraphicalCut = aEllipse;};
  void                 SetPhotopeak(float a, float b){peakPosition = a;peakSigma = b;};
  void                 SetPhotopeakCorrected(float a, float b){peakPositionCorrected = a;peakSigmaCorrected = b;};
  void                 SetVersusTime(TH2F* aHisto){VersusTime = aHisto;};
  void                 SetSimDOIplot(TH2F* aHisto){SimDOIplot = aHisto;};
  void                 SetSimGraph(TGraph* aGraph){SimGraph = aGraph;};
  void                 SetSimFit(TF1* aFit){SimFit = aFit;};
  void                 SetW20percCut(TCut aCut){w20percCut = aCut;};
  void                 SetSlicesMean(TH1D* aHisto){SlicesMean = aHisto;};
  void                 SetCorrectedSpectrum(TH1F* aHisto){CorrectedSpectrum = aHisto;};
//   void                 SetProfileXFit(TF1 aFit){ProfileXFit = aFit;};
  void                 SetSlicesMeanFit(TF1* aFit){SlicesMeanFit = aFit;};
  void                 SetFitCorrected(TF1* aFit){FitCorrected = aFit;};
  void                 SetWbegin(double a){wBegin = a;};
  void                 SetWend(double a){wEnd = a;};
  void                 SetThetaFit(TF1* aFit){ThetaFit = aFit;};
  void                 SetDeltaW(double a){deltaW = a;};
  void                 SetDeltaWfit(TF1* aFit){deltaWfit = aFit;};
  void                 SetDeltaWfit_2(TF1* aFit){deltaWfit_2 = aFit;};
  
  void                 Analyze();
  
  void PrintGlobal();
  void PrintSpecific();
};

#endif  // CRYSTAL_H