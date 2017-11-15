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
#include "TGraphErrors.h"
#include "TGraphDelaunay.h"
#include "TGraph2DErrors.h"


class Crystal : public Element
{

private:
  Element* parentMppc;                       ///< pointer for parent element
  //spectra and co.
  TH1F*                Spectrum;             ///< charge spectrum for this crystal. It's always the sum of all mppcs charges
  TH1F*                singleChargeSpectrum; ///< charge spectrum for this crystal, but not sum of all mppcs, just the sipm in front of the crystal
  TH1F*                LYSpectrum;
  TH1F*                CorrectedSpectrum;    ///< charge spectrum for this crystal corrected by DOI
  TH1F*                HighlightedSpectrum;  ///< same spectrum above, but in green and only for the photopeak
  TH1F*                HighlightedSpectrumCorrected; ///<
  TH1F*                HistoWCorrectedSmooth;
  TH1F*                HistoW;               ///< histogram of w values for this crystal
  TH1F*                HistoWCorrected;
  TH1F*                DensityHisto;         ///< histogram of the entries per voxel in the crystal volume after separation
  TH1F*                pdfW;
  TH1F*                DerivativeDoiResolution;
  TH1F*                cumulativeW;
  TH1F*                resolutions;
  TH1F*                histoAltDoiRes;
  TH1F*                CorrectedSpectrumSearchArea;
  TH1F*                sumChargeSpectrum;
  TH1F*                ATMChargeSpectrum;
  TH1F*                investigatedSpectrum;
  TH1F*                doiResWithCalibration;
  TH1F*                CTRcentralCorrection;
  TH1F*                DeltaTimeWRTTagging;

  TH1D*                SlicesMean;           ///< histogram of fitted mean values of profiles from TH2F ADCvsW distribution
  //   TH1D*                FitSlicesSimDOIplot;  ///< FitSlicesY of the Real Z vs. W plot
  //   TH1D*                FitSlicesSimDOIplotSigma;
  TH2F*                VersusTime;           ///< 2d histogram to plot the evolution of photopeak with time (in case of gain drift?)
  TH2F*                WversusTime;
  TH2F*                SimDOIplot;           ///< 2d histogram for simulation, showing z versus w
  // TH1D*                slicesDelta;
  // TF1*                 slicesDeltaFit;
  TGraph2D***          ComptonCalibation;
  TGraph2D***          ConvertedComptonCalibation;
  TGraphDelaunay***    interpolationGraph;
  TH3I***              ComptonCalibationHistogram;
  std::vector<TH3I*>   ListOfComptonHisto;

  TH3I*                oneHisto;
  TGraph*              SimGraph;
  TCut                 Ellipses;             ///< the elliptical TCut
  TCut                 w20percCut;           ///< TCut from first bin above 20% to last bin above 20% for the w plot
  TCutG               *cutg[2];
  TCut                 crytalCut;
  TCut                 crytalCutWithoutCutG;
  TCut                 photopeakEnergyCut;
  bool                 isOn;                 ///< if the crystal is on/off
  TEllipse             GraphicalCut;         ///< TEllipse to visualize the cut on the u,v global plot
  float                peakPosition;         ///< position of mean (after fitting) for the photopeak
  float                peakSigma;            ///< sigma (after fitting) for the photopeak
  float                LY;
  float                LYSigma;
  float                peakPositionCorrected;         ///< position of mean (after fitting) for the photopeak, corrected by DOI
  float                peakSigmaCorrected;            ///< sigma (after fitting) for the photopeak, corrected by DOI
  TGraph*              calibrationGraph;
  TGraph*              doiResZ;
  TGraph*              wzgraph;
  TGraph*              graphDeltaW;
  TGraph*              graphDeltaRMS;
  TH1F*                simSigmaW;            ///< distribution of sigma W (sigma of gaussian fit of w histogram) for a simulation
  TGraphErrors*        simZvsW;              ///< Z energy deposition point versus peak position of W histogram for a simulation
  TF1*                 Fit;                  ///< fit function (it's a gaussian)
  TF1*                 LYFit;                  ///< fit function (it's a gaussian)
  TF1*                 histoAltDoiFit;


  TF1*                 SimFit;
  TF1*                 Wfit;
  TF1*                 FitCorrected;
  TF1*                 ThetaFit;
  TF1*                 deltaWfit;
  TF1*                 deltaWfit_2;
  //   TF1                  ProfileXFit;
  TF1*                 SlicesMeanFit;
  std::vector<TF1*>    saturationFits;
  double               w_fwhm;               ///< width at half maximum for the w histogram
  double               w_rms;                ///< rms of w histogram
  double               w_width20perc;        ///< width at 20% maximum for the w histogram
  double               u,v,wu,wv,t;
  double               wBegin;               ///< beginning of w histogram after fitting with theta function
  double               wEnd;                 ///< end of w histogram after fitting with theta function
  double               deltaW;               ///< delta of w for a fixed position, as calculated from the gaussian fit of rise in w plot
  double               averageDoiResolution;

  struct multiSpectrum_t
  {
    int canvasPosition;
    TH1F* spectrum;
  };
  struct multiScatter_t
  {
    int canvasPosition;
    TH2F* spectrum;
  };
  struct multiDeltaSlice_t
  {
    int canvasPosition;
    TH1D* spectrum;
    TF1* fit;
  };
  struct multiGraphDelayW_t
  {
    int canvasPosition;
    TGraph* spectrum;
  };
  struct multiGraphDelayRMS_t
  {
    int canvasPosition;
    TGraph* spectrum;
  };


  std::vector<multiSpectrum_t>   deltaTcryTneig;
  std::vector<multiScatter_t>   deltaT2vsW;
  std::vector<multiDeltaSlice_t> slicesDelta;
  std::vector<multiGraphDelayW_t> graphDelayW;
  std::vector<multiGraphDelayRMS_t> graphDelayRMS;

public:
  Crystal();                                 ///< default constructor
  Crystal(const Crystal &obj);               ///< copy constructor
  ~Crystal();                                ///< destructor

  // methods to get and set the private variables. Names should be self explanatory
  Mppc*                GetMppc(){return (Mppc *)parentMppc;};
  TH1F*                GetSpectrum(){return Spectrum;};
  TH1F*                GetLYSpectrum(){return LYSpectrum;};
  TH1F*                GetCorrectedSpectrum(){return CorrectedSpectrum;};
  TH1F*                GetHighlightedSpectrum(){return HighlightedSpectrum;};
  TH1F*                GetHighlightedSpectrumCorrected(){return HighlightedSpectrumCorrected;};
  TH1F*                GetHistoW(){return HistoW;};
  TH1F*                GetHistoWCorrected(){return HistoWCorrected;};
  TH1F*                GetHistoWCorrectedSmooth(){return HistoWCorrectedSmooth;};
  TH1F*                GetDensityHisto(){return DensityHisto;};
  TH1F*                GetPdfW(){return pdfW;};
  TH1F*                GetCumulativeW(){return cumulativeW;};
  TH1F*                GetDerivativeDoiResolution(){return DerivativeDoiResolution;};
  TH1F*                GetDoiResolutions(){return resolutions;};
  TH1F*                GetHistoAltDoiRes(){return histoAltDoiRes;};
  TH1F*                GetSingleChargeSpectrum(){return singleChargeSpectrum;};
  TH1F*                GetSumChargeSpectrum(){return sumChargeSpectrum;};
  TH1F*                GetAMTChargeSpectrum(){return ATMChargeSpectrum;};
  TH1F*                GetInvestigatedSpectrum(){return investigatedSpectrum;};
  TH1F*                GetDoiResWithCalibration(){return doiResWithCalibration;};
  TGraph*              GetCalibrationGraph(){return calibrationGraph;};
  TGraph*              GetDoiResZ(){return doiResZ;};
  TF1*                 GetFit(){return Fit;};
  TF1*                 GetLYFit(){return LYFit;};
  TF1*                 GetSimFit(){return SimFit;};
  TF1*                 GetHistoWfit(){return Wfit;};
  TF1*                 GetHistoAltDoiFit(){return histoAltDoiFit;};
  TGraph2D***          GetComptonCalibration(){return ComptonCalibation;};
  TGraph2D***          GetConvertedComptonCalibration(){return ConvertedComptonCalibation;};
  TGraphDelaunay***    GetInterpolationGraph(){return interpolationGraph;};
  TH3I***              GetComptonCalibrationHistogram(){return ComptonCalibationHistogram;};

  TH3I*                GetOneHisto(){return oneHisto;};
  TGraph*              GetWZgraph(){return wzgraph;};
  //   TF1*                 GetProfileXFit(){return &ProfileXFit;};
  TF1*                 GetSlicesMeanFit(){return SlicesMeanFit;};
  double               GetWfwhm(){return w_fwhm;};
  double               GetWrms(){return w_rms;};
  double               GetWwidth20perc(){return w_width20perc;};
  // TCut                 GetCrystalCut(){return Ellipses;};
  TEllipse*            GetGraphicalCut(){return &GraphicalCut;};
  float                GetPhotopeakPosition(){return peakPosition;};
  float                GetPhotopeakSigma(){return peakSigma;};
  float                GetLY(){return LY;};
  float                GetLYSigma(){return LYSigma;};
  float                GetPhotopeakEnergyResolution(){return ((peakSigma*2.355)/peakPosition);};
  float                GetPhotopeakPositionCorrected(){return peakPositionCorrected;};
  float                GetPhotopeakSigmaCorrected(){return peakSigmaCorrected;};
  float                GetPhotopeakEnergyResolutionCorrected(){return ((peakSigmaCorrected*2.355)/peakPositionCorrected);};
  double               GetWbegin(){return wBegin;};
  double               GetWend(){return wEnd;};
  double               GetDeltaW(){return std::abs(deltaW);};

  bool                 CrystalIsOn(){return isOn;};
  TH2F*                GetVersusTime(){return VersusTime;};
  TH2F*                GetWversusTime(){return WversusTime;};
  TH2F*                GetSimDOIplot(){return SimDOIplot;};
  //   TH1D*                GetFitSlicesSimDOIplot(){return FitSlicesSimDOIplot;};
  //   TH1D*                GetFitSlicesSimDOIplotSigma(){return FitSlicesSimDOIplotSigma;};
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
  TCut                 GetCrystalCut(){return crytalCut;};
  TCut                 GetCrystalCutWithoutCutG(){return crytalCutWithoutCutG;};
  TCut                 GetPhotopeakEnergyCut(){return photopeakEnergyCut;};

  TF1*                 GetThetaFit(){return ThetaFit;};
  TF1*                 GetDeltaWfit(){return deltaWfit;};
  TF1*                 GetDeltaWfit_2(){return deltaWfit_2;};
  double               GetMcal(){return dz/(wBegin - wEnd);};
  double               GetQcal(){return -(dz*wEnd)/(wBegin-wEnd);};
  double               GetDoiResolutionFWHM(){return 2.355 * std::abs((dz/(wBegin - wEnd))) * std::abs(deltaW);};
  double               GetAverageDoiResolution(){return averageDoiResolution;};
  TH1F*                GetSimSigmaW(){return simSigmaW;};
  TGraphErrors*        GetSimZvsW(){return simZvsW;};
  TH3I*                GetComptonHistogram(int i){return ListOfComptonHisto[i];};
  int                  GetNumOfComptonHisto(){return ListOfComptonHisto.size();};
  TH1F*                GetCorrectedSpectrumSearchArea(){return CorrectedSpectrumSearchArea;};
  TH1F*                GetCTRcentralCorrection(){return CTRcentralCorrection;};
  TGraph*              GetGraphDeltaW(){return graphDeltaW;};
  TGraph*              GetGraphDeltaRMS(){return graphDeltaRMS;};

  std::vector<TF1*>    GetSaturationFits(){return saturationFits;};
  std::vector<multiSpectrum_t>             GetDeltaTcryTneig(){return deltaTcryTneig;};
  std::vector<multiScatter_t>              GetDeltaT2vsW(){return deltaT2vsW;};
  std::vector<multiDeltaSlice_t>           GetDeltaSlice(){return slicesDelta;};
  std::vector<multiGraphDelayW_t>          GetGraphDelayW(){return graphDelayW;};
  std::vector<multiGraphDelayRMS_t>          GetGraphDelayRMS(){return graphDelayRMS;};

  TH1F*                GetDeltaTimeWRTTagging()                  {return DeltaTimeWRTTagging;};

  void                 SetDeltaTimeWRTTagging(TH1F* aHisto)       {DeltaTimeWRTTagging = aHisto;};
  void                 SetCorrectedSpectrumSearchArea(TH1F * aHisto){CorrectedSpectrumSearchArea = aHisto;};
  void                 SetZXCut(TCutG *aCut){cutg[0] = aCut;};
  void                 SetZYCut(TCutG *aCut){cutg[1] = aCut;};
  void                 SetCrystalCut(TCut aCut){crytalCut = aCut;};
  void                 SetCrystalCutWithoutCutG(TCut aCut){crytalCutWithoutCutG = aCut;};
  void                 SetPhotopeakEnergyCut(TCut aCut){photopeakEnergyCut = aCut;};
  void                 SetMppc(Mppc *amppc);
  void                 SetSpectrum(TH1F* aHisto){Spectrum = aHisto;};
  void                 SetLYSpectrum(TH1F* aHisto){LYSpectrum = aHisto;};
  void                 SetHighlightedSpectrum(TH1F* aHisto){HighlightedSpectrum = aHisto;};
  void                 SetHighlightedSpectrumCorrected(TH1F* aHisto){HighlightedSpectrumCorrected = aHisto;};
  void                 SetHistoW(TH1F* aHisto){HistoW = aHisto;};
  void                 SetHistoWCorrected(TH1F* aHisto){HistoWCorrected = aHisto;};
  void                 SetHistoWCorrectedSmooth(TH1F* aHisto){HistoWCorrectedSmooth = aHisto;};
  void                 SetDerivativeDoiResolution(TH1F* aHisto){DerivativeDoiResolution = aHisto;};
  void                 SetDoiResolutions(TH1F* aHisto){resolutions = aHisto;};
  void                 SetDensityHisto(TH1F* aHisto){DensityHisto = aHisto;};
  void                 SetFit(TF1* aFit){Fit = aFit;};
  void                 SetLYFit(TF1* aFit){LYFit = aFit;};
  void                 SetHistoWfwhm(double a){w_fwhm = a;};
  void                 SetHistoWrms(double a){w_rms = a;};
  void                 SetHistoWwidth20perc(double a){w_width20perc = a;};
  void                 SetHistoWfit(TF1* aFit){Wfit = aFit;};
  void                 SetHistoAltDoiFit(TF1* aFit){histoAltDoiFit = aFit;};
  void                 SetPdfW(TH1F* aHisto){pdfW = aHisto;};
  void                 SetDoiResWithCalibration(TH1F* aHisto){doiResWithCalibration = aHisto;};
  void                 SetCumulativeW(TH1F* aHisto){cumulativeW = aHisto;};
  void                 SetCalibrationGraph(TGraph* aGraph){calibrationGraph = aGraph;};
  void                 SetComptonCalibration(TGraph2D*** aGraph){ComptonCalibation = aGraph;};
  void                 SetConvertedComptonCalibration(TGraph2D*** aGraph){ConvertedComptonCalibation = aGraph;};
  void                 SetInterpolationGraph(TGraphDelaunay*** aGraph){interpolationGraph = aGraph;};
  void                 SetComptonCalibrationHistogram(TH3I*** aHisto){ComptonCalibationHistogram = aHisto;};
  void                 SetOne(TH3I* aHisto){oneHisto = aHisto;};
  void                 SetWZgraph(TGraph* aGraph){wzgraph = aGraph;};
  //   void                 SetEllipses(std::string varX,std::string varY);
  void                 SetCrystalOn(bool abool){isOn = abool;};
  void                 SetCrystalData(double au,double av,double awu ,double awv, double at){u = au; v = av; wu = awu ; wv = awv ; t = at;};
  void                 SetGraphicalCut(TEllipse aEllipse){GraphicalCut = aEllipse;};
  void                 SetPhotopeak(float a, float b){peakPosition = a;peakSigma = b;};
  void                 SetLY(float a, float b){LY = a;LYSigma = b;};
  void                 SetPhotopeakCorrected(float a, float b){peakPositionCorrected = a;peakSigmaCorrected = b;};
  void                 SetVersusTime(TH2F* aHisto){VersusTime = aHisto;};
  void                 SetWversusTime(TH2F* aHisto){WversusTime = aHisto;};
  void                 SetSimDOIplot(TH2F* aHisto){SimDOIplot = aHisto;};
  //   void                 SetFitSlicesSimDOIplot(TH1D* aHisto){FitSlicesSimDOIplot = aHisto;};
  //   void                 SetFitSlicesSimDOIplotSigma(TH1D* aHisto){FitSlicesSimDOIplotSigma = aHisto;};
  void                 SetSimZvsW(TGraphErrors* aGraph){simZvsW = aGraph;};
  void                 SetSimSigmaW(TH1F* aHisto){simSigmaW = aHisto;};
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
  void                 SetDoiResZ(TGraph* aGraph){doiResZ = aGraph;};
  void                 SetAverageDoiResolution(double a){averageDoiResolution = a;};
  void                 SetComptonHistogram(TH3I* aHisto){ListOfComptonHisto.push_back(aHisto);};
  void                 SetHistoAltDoiRes(TH1F *aHisto){histoAltDoiRes = aHisto;};
  void                 SetSingleChargeSpectrum(TH1F *aHisto){singleChargeSpectrum = aHisto;};
  void                 SetSumChargeSpectrum(TH1F *aHisto){sumChargeSpectrum = aHisto;};
  void                 SetAMTChargeSpectrum(TH1F *aHisto){ATMChargeSpectrum = aHisto;};
  void                 SetInvestigatedSpectrum(TH1F *aHisto){investigatedSpectrum = aHisto;};
  void                 SetSaturationFits(std::vector<TF1*> aVec){saturationFits = aVec;};
  void                 SetCTRcentralCorrection(TH1F* aHisto){CTRcentralCorrection = aHisto;};
  void                 Analyze();
  void                 SetGraphDeltaW(TGraph* aGraph){graphDeltaW = aGraph;};
  void                 SetGraphDeltaRMS(TGraph* aGraph){graphDeltaRMS = aGraph;};


  void                 AddDeltaTcryTneig(TH1F* aHisto,int aPos)
  {
    multiSpectrum_t tempSpectrum;
    tempSpectrum.canvasPosition = aPos;
    tempSpectrum.spectrum = aHisto;
    deltaTcryTneig.push_back(tempSpectrum);
  };
  void                 AddDeltaT2vsW(TH2F* aHisto,int aPos)
  {
    multiScatter_t tempSpectrum;
    tempSpectrum.canvasPosition = aPos;
    tempSpectrum.spectrum = aHisto;
    deltaT2vsW.push_back(tempSpectrum);
  };
  // void                 AddSlicesDelta(TH1D* aHisto,TF1* aFit,int aPos)
  void                 AddSlicesDelta(TH1D* aHisto,int aPos)
  {
    multiDeltaSlice_t tempDelta;
    tempDelta.canvasPosition = aPos;
    tempDelta.spectrum       = aHisto;
    // tempDelta.fit            = aFit;
    slicesDelta.push_back(tempDelta);
  };
  void                 AddGraphDelayW(TGraph* aGraph,int aPos)
  {
    multiGraphDelayW_t tempDelta;
    tempDelta.canvasPosition = aPos;
    tempDelta.spectrum       = aGraph;
    // tempDelta.fit            = aFit;
    graphDelayW.push_back(tempDelta);
  };

  void AddGraphDelayRMS(TGraph* aGraph,int aPos)
  {
    multiGraphDelayRMS_t tempDelta;
    tempDelta.canvasPosition = aPos;
    tempDelta.spectrum       = aGraph;
    // tempDelta.fit            = aFit;
    graphDelayRMS.push_back(tempDelta);
  }


  void PrintGlobal();
  void PrintSpecific();
};

#endif  // CRYSTAL_H
