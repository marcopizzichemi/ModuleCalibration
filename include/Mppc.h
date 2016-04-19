#ifndef MPPC_H
#define MPPC_H
#include "Element.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH3I.h"
#include "TCutG.h"
#include "TProfile.h"
#include "TF1.h"

#include <vector>
#include <iostream>
// #include <string>

class Mppc : public Element
{
private:
  Element*               parentModule;           ///< one pointer for its parent, it's only one
  std::vector<Element*>  vCrystal;               ///< vector of pointers for the children crystals
  int                    digitizerChannel;       ///< which digitizer channel is assigned to this mppc
  int                    canvasPosition;         ///< position in the canvas of all channels
  std::string            moduleName;             ///< name of the module
  bool                   IsOnForDoi;             ///< if the mppc is on the DOI tagging line
  double                 ThetaWU;                ///< angle from w to u in radiants
  double                 ThetaWV;                ///< angle from w to v in radiants
  
  //histograms
  TH1F*                  RawSpectrum;            ///< raw spectrum of all events seen by this mppc
  TH1F*                  TriggerSpectrum;        ///< trigger spectrum of all events seen by this mppc
  TH1F*                  TriggerSpectrumHighlighted;///< trigger spectrum of all events seen by this mppc, highlighting the broad energy cut
  TH2D*                  projection_zy;          ///< Projection histogram of u,v,w points on the w,v plane
  TH2D*                  projection_zx;          ///< Projection histogram of u,v,w points on the w,u plane
//   TH1D*                  projection_x;           ///< 
//   TH1D*                  projection_y;           ///<
  TProfile*              profileX;               ///< Profile of w,u histogram (for each bin in w, mean u and sigma are plotted)
  TProfile*              profileY;               ///< Profile of w,v histogram (for each bin in w, mean v and sigma are plotted)
  TF1*                   lineX;                  ///< Line to fit the profileX plot -> u(w) = m*w + c
  TF1*                   lineY;                  ///< Line to fit the profileX plot -> v(w) = m*w + c
  
  std::vector<double>    fit2DmeanX ;            ///< arrays of mean and sigma for the 2d search of peaks in this mppc
  std::vector<double>    fit2DmeanY ;            ///< arrays of mean and sigma for the 2d search of peaks in this mppc
  std::vector<double>    fit2DsigmaX;            ///< arrays of mean and sigma for the 2d search of peaks in this mppc
  std::vector<double>    fit2DsigmaY;            ///< arrays of mean and sigma for the 2d search of peaks in this mppc
  std::vector<double>    fit2Dtheta;             ///< arrays of mean and sigma for the 2d search of peaks in this mppc
  struct point
  {
    int i;
    int j;
    int k;
  };
  
  struct masks_t
  {
    double meanx;
    double meany;
    double meanz;
    int maskID;
    int maskI;
    int maskJ;
    long int nBinsXMask;
//     bool operator<(const masks_t& rhs) const { meanx < rhs.meanx; }
  };
  
  
  struct compare_by_x
  {
    bool operator()(const masks_t& lhs, const masks_t& rhs) const
    {
      return lhs.meanx < rhs.meanx;
    }
  };

  struct compare_by_y
  {
    bool operator()(const masks_t& lhs, const masks_t& rhs) const
    {
      return lhs.meany < rhs.meany;
    }
  };
  
public:
  Mppc();                                        ///< default constructor
  Mppc(const Mppc &obj);                         ///< copy constructor
  ~Mppc();                                       ///< destructor
  
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
  TH1F*                  GetRawSpectrum(){return RawSpectrum;};
  void                   SetRawSpectrum(TH1F* aHisto){RawSpectrum = aHisto;};
  TH1F*                  GetTriggerSpectrum(){return TriggerSpectrum;};
  void                   SetTriggerSpectrum(TH1F* aHisto){TriggerSpectrum = aHisto;};
  TH1F*                  GetTriggerSpectrumHighlighted(){return TriggerSpectrumHighlighted;};
  void                   SetTriggerSpectrumHighlighted(TH1F* aHisto){TriggerSpectrumHighlighted = aHisto;};
  void                   SetIsOnForDoi(bool abool){IsOnForDoi = abool;};
  bool                   GetIsOnForDoi(){return IsOnForDoi;};
  std::vector<double>*   GetFit2DmeanX(){return &fit2DmeanX;};
  std::vector<double>*   GetFit2DmeanY(){return &fit2DmeanY;};
  std::vector<double>*   GetFit2DsigmaX(){return &fit2DsigmaX;};
  std::vector<double>*   GetFit2DsigmaY(){return &fit2DsigmaY;};
  std::vector<double>*   GetFit2Dtheta(){return &fit2Dtheta;};
  TProfile*              GetProfileX(){return profileX;};
  TProfile*              GetProfileY(){return profileY;};
  TH2D*                  GetProjectionZX(){return projection_zx;};
  TH2D*                  GetProjectionZY(){return projection_zy;};
  double                 GetThetaWU(){return ThetaWU;};
  double                 GetThetaWV(){return ThetaWV;};
  
  // methods to analyze the mppc
  int                    Find2Dpeaks(int nofcrystals,TH2F* histogram2d); ///< Finds the 2D peaks for crystals coupled to this module
  void                   FindProjectionPlane();                          ///< Finds the best projection plane for this mppc
  
  void                   MakeRotatedFlood();
  bool                   FindCrystalCuts(TCutG**** cutg,int histo3DchannelBin, int div,int nofcrystalsx,int nofcrystalsy);///< Finds the density volumes that represents the crystals
  
//   bool compare_by_x(const masks_t& lhs, const masks_t& rhs) { return lhs.meanx < rhs.meanx; };
//   bool compare_by_y(const masks_t& lhs, const masks_t& rhs) { return lhs.meany < rhs.meany; };
  
  
  
  
  //print methods
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // MPPC_H
