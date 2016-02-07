#ifndef MPPC_H
#define MPPC_H
#include "Element.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
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
  bool                   IsOnForDoi;             ///<
  double                 ThetaWU;                     ///< angle from w to u in radiants
  double                 ThetaWV;                     ///< angle from w to v in radiants
  
  //histograms
  TH1F                   RawSpectrum;            ///< raw spectrum of all events seen by this mppc
  TH1F                   TriggerSpectrum;        ///< raw spectrum of all events seen by this mppc
  TH2D*                  projection_zy;          ///<
  TH2D*                  projection_zx;          ///<
  TH1D*                  projection_x;           ///<
  TH1D*                  projection_y;           ///<
  TProfile*              profileX;               ///<
  TProfile*              profileY;               ///<
  TF1*                   lineX;                  ///<
  TF1*                   lineY;                  ///<
  
  std::vector<double>    fit2DmeanX ;            ///< arrays of mean and sigma for the 2d search of peaks in this mppc
  std::vector<double>    fit2DmeanY ;            ///< arrays of mean and sigma for the 2d search of peaks in this mppc
  std::vector<double>    fit2DsigmaX;            ///< arrays of mean and sigma for the 2d search of peaks in this mppc
  std::vector<double>    fit2DsigmaY;            ///< arrays of mean and sigma for the 2d search of peaks in this mppc
  std::vector<double>    fit2Dtheta;             ///< arrays of mean and sigma for the 2d search of peaks in this mppc
  
  
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
  TH1F*                  GetRawSpectrum(){return &RawSpectrum;};
  void                   SetRawSpectrum(TH1F aHisto){RawSpectrum = aHisto;};
  TH1F*                  GetTriggerSpectrum(){return &TriggerSpectrum;};
  void                   SetTriggerSpectrum(TH1F aHisto){TriggerSpectrum = aHisto;};
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
  
  //print methods
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // MPPC_H