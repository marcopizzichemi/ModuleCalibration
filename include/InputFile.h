// InputFile.h
// this class reads the input tchain and generates the TTree that will be used for the analysis /plot production
// Reading of config file and creation of TChain is performed directly in the class constructor, then CreateTree
// makes the TTree and FillElements is used to fill the modules, mppcs, and crystals with the information stored in
// the config file.

#ifndef INPUTFILE_H
#define INPUTFILE_H

#include "TChain.h"
#include "TTree.h"
#include <string>
#include "ConfigFile.h"
#include "Element.h"
#include "Crystal.h"
#include "Module.h"
#include "Mppc.h"
#include "TRandom3.h"
// #include "Detector.h"


/** @brief Class that controls the input files
 *
 *  This class reads the input tchain and generates the TTree that will be used for the analysis /plot production
 *  Reading of config file and creation of TChain is performed directly in the class constructor, then CreateTree
 *  makes the TTree and FillElements is used to fill the modules, mppcs, and crystals with the information stored in
 *  the config file.
 *
 *  @author M. Pizzichemi
 *  @date Nov 2015
 */


class InputFile
{

private:

  TChain*                        fchain;                             // pointer to input TChain
  std::string                    fname;                              // name of the input TChain
  TTree*                         ftree;                              // pointer to analysis TTree
  int                            inputChannels;                      // number of input channels
  int                            adcChannels;                        // number of output channels in the adc in use (32 for the CAEN DT5740)
  bool                           *DigitizerChannelOn;                // for each channel, if it's on or not. Meaning if it is specified in the config file "digitizer" key
  //variables read from the config file
  std::string                    ConfigFileName;                     // name of the config file
  //   std::string                    chainName;
  std::string                    digitizer_s;                        // input string of the digitizer key from config file
  std::string                    mppc_s;                             // input string of the mppc key from config file
  std::string                    plotPositions_s;                    // input string of the plotPositions key from config file
  std::string                    xPositions_s;                       // input string of the xPositions key from config file
  std::string                    yPositions_s;                       // input string of the yPositions key from config file
  std::string                    saturation_s;                       // input string of the saturation key from config file
  std::string                    pedestal_s;                       // input string of the saturation key from config file
  std::string                    noise_s;
  std::string                    gain_s;
  std::string                    mppcOFF_s;
  std::vector <std::string>      mppcOFF_f;
  std::vector <std::string>      mppcOFF;
  std::string                    crystalOFF_s;
  std::vector <std::string>      crystalOFF_f;
  std::vector <int>              crystalOFF;

  std::string                    specificMPPC_s;
  std::vector <std::string>      specificMPPC_f;
  std::vector <std::string>      specificMPPC;

  std::string                    specificBin_s;
  std::vector <std::string>      specificBin_f;
  std::vector <int>              specificBin;

  std::string                    specificPrecision_s;
  std::vector <std::string>      specificPrecision_f;
  std::vector <int>           specificPrecision;

  std::string                    specificCut_s;
  std::vector <std::string>      specificCut_f;
  std::vector <double>           specificCut;

  std::vector <std::string>      digitizer_f;                        // tokenization of above strings
  std::vector <std::string>      mppc_f;                             // tokenization of above strings
  std::vector <std::string>      plotPositions_f;                    // tokenization of above strings
  std::vector <std::string>      xPositions_f;                       // tokenization of above strings
  std::vector <std::string>      yPositions_f;                       // tokenization of above strings
  std::vector <std::string>      saturation_f;                       // tokenization of above strings
  std::vector <std::string>      pedestal_f;                         // tokenization of above strings
  std::vector <std::string>      noise_f;                         // tokenization of above strings
  std::vector <std::string>      gain_f;                         // tokenization of above strings
  std::vector <int>              digitizer;                          // above tokenized string transformed in the proper variable types


  std::vector <int>              pedestal;                          // above tokenized string transformed in the proper variable types
  std::vector <int>              noise;                          // above tokenized string transformed in the proper variable types
  std::vector <float>              gain;                          // above tokenized string transformed in the proper variable types
  std::vector <std::string>      mppc_label;                         // above tokenized string transformed in the proper variable types
  std::vector <int>              plotPositions;                      // above tokenized string transformed in the proper variable types
  std::vector <float>            xPositions;                         // above tokenized string transformed in the proper variable types
  std::vector <float>            yPositions;                         // above tokenized string transformed in the proper variable types
  std::vector <float>            saturation;                         // above tokenized string transformed in the proper variable types

  std::string                    timingCh_s;
  std::vector <std::string>      timingCh_f;
  std::vector <int>              timingCh;                          // above tokenized string

  std::string                    **crystal_s;                        // array of string, one for each crystal input by the user in the config file
  std::vector <std::string>      crystal_f;                          // tokenized version
  std::string                    digitizerDoi_s;                        // input string of the digitizer key from config file
  std::vector <std::string>      digitizerDoi_f;                        // tokenization of above strings
  std::vector <int>              digitizerDoi;                          // above tokenized string transformed in the proper variable types

  float                          ***crystaldata;                     // i,j matrix with data of the above string transformed in float
  bool                           **crystalIsOn;                      // i,j matrix with crystal ON/OFF
  // variable for the module(s) structure
  int                            ncrystalsx;                         // number of crystals in x direction per mppc
  int                            ncrystalsy;                         // number of crystals in y direction per mppc
  int                            nmppcx;                             // number of mppc in x direction per mppc
  int                            nmppcy;                             // number of mppc in y direction per mppc
  int                            nmodulex;                           // number of modules in x direction per mppc
  int                            nmoduley;                           // number of modules in y direction per mppc
  //
  std::string                    BinaryOutputFileName;               // name of the output binary file (if selected)
  bool                           binary;                             // true if the user wants a binary output, false otherwise
  bool                           correctingSaturation;               // true if the input dataset will be corrected for saturation
  bool                           saturationRun;                      // whether this is a saturation run or not
  float                          taggingPosition;                    // tagging DOI bench position
  bool                           usingTaggingBench;                  // true if the input datasets uses the DOI tagging
  bool                           usingRealSimData;                   // true if the "real" gamma interaction point is used (of course valid only for simulation datasets)
  bool                           usingAllChannels;
  bool                           wAllChannels;
  int                            taggingCrystalChannel;              // channel of the tagging crystal, only for DOI bench data
  int                            taggingCrystalTimingChannel;
  double                         nclock;                             // number of clock samples that will be ignored. clock is the digitizer clock (so 1 sample = 16ns for DT5740)
  int                            *translateCh;                       // translates adc channels to analysis channels
  double                         crystalx;                           // dimension of crystal in x [mm]
  double                         crystaly;                           // dimension of crystal in x [mm]
  double                         crystalz;                           // dimension of crystal in x [mm]
  double                         esrThickness;                       // thickness of esr separation foil [mm]
  double                         chargeBinningADC;                   // adc charge binning
  int                            saturationFormat;                   // format of saturation input. It can be in units of ADC_CHANNELS or CHARGE
  int digitizerType;
  float approximateTDCbinning;
  int TDCcalculationEntries;
  int                            global_histo3DchannelBin;
  float                          global_histo3Dmin;
  float                          global_histo3Dmax;
  int                            global_div;
  double                         global_clusterVolumeCut;

  long long int                  GoodCounter;
  long long int                  badEvents;
  long long int                  counter;
  bool calculateTDCbinning;
  float minDeltaForFT;
  float pedestalTag;
  bool taggingForTiming;

  //variables for the input TChain
  ULong64_t     ChainExtendedTimeTag;                                // extended time tag
  ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous event
  Int_t        *ChainAdcChannel;
  Short_t      *ChainDesktopAdcChannel;                              // input TChain data for desktop digitizers - data is int_16
  UShort_t     *ChainVMEadcChannel;                                  // input TChain data for VME digitizers - data is uint_16
  UShort_t     *ChainVMEAmplitudeChannel;                             // input TChain data for VME amplitude - data is also uint_16
  Float_t      *ChainTimeStamp;
  Float_t      *TDCBinning;
  // Short_t      *ChainPetirocChannel;                                 //FIXME temporary data type of petiroc charge input - ask
  Float_t       RealX;                                               // "real" gamma interaction positions (from simulation data)
  Float_t       RealY;                                               // "real" gamma interaction positions (from simulation data)
  Float_t       RealZ;                                               // "real" gamma interaction positions (from simulation data)
  Float_t       simTaggingCharge;
  Float_t       simTaggingTime;
  Short_t       CrystalsHit;                                         // "real" number of crystals hit in the event (from simulation data)
  Short_t       NumbOfInteractions;                                  // "real" number of interaction (energy depositions) in the event (from simulation data)
  std::vector <float>* TotalCryEnergy;

  //branches for the input TChain
  TBranch      *bChainExtendedTimeTag;                               // branches for above data
  TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  TBranch     **bChainAdcChannel;                                    // branches for above data
  TBranch     **bChainAmplChannel;                                    // branches for above data
  TBranch     **bChainTimeStamp;                                    // branches for above data
  TBranch      *bRealX;                                              // branches for above data
  TBranch      *bRealY;                                              // branches for above data
  TBranch      *bRealZ;                                              // branches for above data
  TBranch      *bsimTaggingCharge;                                              // branches for above data
  TBranch      *bsimTaggingTime;                                              // branches for above data
  TBranch      *bCrystalsHit;                                        // branches for above data
  TBranch      *bNumbOfInteractions;                                 // branches for above data
  TBranch      *bTotalCryEnergy;                                     //

  //variables for the analysis TTree
  ULong64_t     TreeExtendedTimeTag;                                 // extended time tag
  ULong64_t     TreeDeltaTimeTag;                                    // delta tag from previous event
  Float_t       *TreeAdcChannel;                                      // channels data for this event
  int           TreeTriggerChannel;                                  // trigger channel for this event
  Float_t       TreeTagging;                                         // tagging crystal data for this event
  Float_t       TaggingTimeStamp;
  Float_t       TreeFloodX;                                          // u position for this event
  Float_t       TreeFloodY;                                          // v position for this event
  Float_t       TreeFloodZ;                                          // w position for this event
  // Float_t       TreeTheta;                                           // u position for this event
  // Float_t       TreePhi;                                             // u position for this event
  Float_t      *TreeTimeStamp;                                       // timestamps data for each channel
  Bool_t        TreeBadevent;                                        // whether the event is "bad" --> too high to allow saturation correction with a logarithm
  Float_t       TreeZPosition;                                       // z position of the tagging bench
  Float_t       TreeRealX;                                           // "real" gamma interaction positions (from simulation data)
  Float_t       TreeRealY;                                           // "real" gamma interaction positions (from simulation data)
  Float_t       TreeRealZ;                                           // "real" gamma interaction positions (from simulation data)
  Short_t       TreeCrystalsHit;                                     // "real" number of crystals hit in the event (from simulation data)
  Short_t       TreeNumbOfInteractions;                              // "real" number of interaction (energy depositions) in the event (from simulation data)
  std::vector <float> TreeTotalCryEnergy;
  std::vector <float>* pTreeTotalCryEnergy;

  bool smearTaggingTime ;                                            // whether to smear the time stamp of external tagging. Needed for simulations, where the tagging time stamp is always 0 (i.e. the gamma emission time) - default = 0
  float sigmaTimeTag;                                                // sigma for the smearing of tagging time [ps]. it's the time resolution of an hypothetical external short crystal + fast sipm - default = 30.0, which corresponds to Hamamatsu MPPC + 2x2x3 m3 LSO-Ca codoped crystal (100ps FWHM CTR, see Mythra poster)
  TRandom3 *randGen;

  // struct detector_t
  // {
  //   int digitizerChannel;
  //   int timingChannel;
  //   std::string label;
  //   int i;
  //   int j;
  //   float saturation;
  //   int plotPosition;
  //   float xPosition;
  //   float yPosition;
  //   float pedestal;
  //   int OnForDOI;
  //   bool isNeighbour;
  //   std::vector<int> neighbourChannels;
  //   bool OnForModular;
  //   //     bool operator<(const masks_t& rhs) const { meanx < rhs.meanx; }
  // };

  std::vector<detector_t> detector;


  class inputPetirocFile_t
  {
  public:

    int Nchannels;
    std::vector<int> FineTime;
    std::vector<int> Charge;
    std::vector<int> CoarseTime;
    inputPetirocFile_t(int a){ Nchannels = a;};
    inputPetirocFile_t(){};
    void clear()
    {
      FineTime.clear();
      Charge.clear();
      CoarseTime.clear();
    };
    friend std::istream& operator>>(std::istream& input, inputPetirocFile_t& s)
    {
      for(int p = 0; p < s.Nchannels; p++)
      {
        int a,b;
        input >> a >> b;
        s.FineTime.push_back(a);
        s.Charge.push_back(b);
      }
      for(int p = 0; p < s.Nchannels; p++)
      {
        int a,b;
        input >> a >> b;
        s.CoarseTime.push_back(a);
      }
      return input;
    }
  };

public:

  InputFile(ConfigFile& config);              //ctor

  static InputFile*  Instance() { return fInstance; };               // not useful right now
  static InputFile*  fInstance;                                      // not useful right now

  TChain*       GetChain() const { return fchain; };                 // method to provide a pointer to the input TChain
  TTree*        GetTree()  const { return ftree; };                   // method to provide a pointer to the analysis TTree
  void          SetTree(TTree *aTree){ftree = aTree;};
  void          ImportTChain(int argc, char** argv);
  // void          ImportPetirocData(int argc, char** argv);
  void          PrepareTTree();
  void          FillTree(int argc, char** argv);                                        // method to run on the input and fill the analysis TTree
  void          FillElements(Module*** module,Mppc*** mppc,Crystal*** crystal);  // method to fill with info the elements (modules, mppcs, crystals)
  void          FillTreeCAEN();
  void          FillEvent();
  void          FillTreePetiroc(int argc, char** argv);
};




#endif  // INPUTFILE_H
