#ifndef MODULE_H
#define MODULE_H

#include <iostream>
#include <string>
#include <vector>
#include "Element.h"
#include "TTreeFormula.h"
// #include "Detector.h"

class Module : public Element
{
private:
  std::vector<Element*>  vMppc;
  std::vector<int> channels;
  std::vector<int> detChData;
  std::vector<float> saturationData;
  std::vector<float> pedestalData;

  UInt_t seed;
  std::vector<detector_t> detector;
  int taggingCrystalTimingChannel;
  TH2F *TagVsExtTimeSpectrum;
  TH2F *TagVsDeltaTimeSpectrum;
  TTreeFormula* FormulaTag;

public:
  Module(); // default constructor
  Module(const Module &obj); // copy constructor
  ~Module(); // destructor

  void                   SetMppc(Mppc *pMppc);
  void                   SetDetector(std::vector<detector_t> aDetector){detector = aDetector;};
  std::vector<detector_t>   GetDetector(){return detector;};
  void                   SetChannels(std::vector<int> aVec){channels = aVec;};
  void                   SetTaggingTimingChannel (int aNum){taggingCrystalTimingChannel = aNum;};
  std::vector<int>       GetChannels(){return channels;};

  void                   SetDetChData(std::vector<int> aVec){detChData = aVec;};
  void                   SetSaturationData(std::vector<float> aVec){saturationData = aVec;};
  void                   SetPedestalData(std::vector<float> aVec){pedestalData = aVec;};

  std::vector<int>       GetDetChData(){return detChData;};
  std::vector<float>     GetSaturationData(){return saturationData;};
  std::vector<float>     GetPedestalData(){return pedestalData;};

  int                    GetMppcsNumber(){return vMppc.size();};
  Mppc*                  GetMppc(int pi, int pj);
  void                   SetSeed(UInt_t aSeed){seed = aSeed;};
  UInt_t                 GetSeed(){return seed;};
  int                    GetTaggingTimingChannel(){return taggingCrystalTimingChannel;};
  void                   SetTagVsExtTimeSpectrum(TH2F* aHisto){TagVsExtTimeSpectrum = aHisto;};
  TH2F*                  GetTagVsExtTimeSpectrum(){return TagVsExtTimeSpectrum;};

  void                   SetTagVsDeltaTimeSpectrum(TH2F* aHisto){TagVsDeltaTimeSpectrum = aHisto;};
  TH2F*                  GetTagVsDeltaTimeSpectrum(){return TagVsDeltaTimeSpectrum;};
  void                   SetFormulaTaggingPhotopeakCut(TTreeFormula* aFormula){FormulaTag = aFormula;};
  TTreeFormula*          GetFormulaTaggingPhotopeakCut(){return FormulaTag;};
  void PrintGlobal();
  void PrintSpecific();
};

#endif  // MODULE_H
