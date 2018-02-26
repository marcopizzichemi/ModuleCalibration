// class of input points from doi tag bench
class inputDoi_t
{
public:
  int i;
  int j;
  double m;
  double q;
  double doires;
  double avgs;
  std::vector<double> w;
  std::vector<double> sw;
  std::vector<double> sqrt_nentries;
  std::vector<double> z;
  std::vector<double> sz;
  int pointsFromDoi;
  inputDoi_t(int a){ pointsFromDoi = a;};
  inputDoi_t(){};
  void clear()
  {
    w.clear();
    sw.clear();
    sqrt_nentries.clear();
    z.clear();
    sz.clear();
  };
  void setPointsFromDoi(int a) {pointsFromDoi = a;};
  friend std::istream& operator>>(std::istream& input, inputDoi_t& s)
  {
    input >> s.i;
    input >> s.j;
    input >> s.m;
    input >> s.q;
    input >> s.doires;
    input >> s.avgs;
    for(int p = 0; p < s.pointsFromDoi; p++)
    {
      double wValue,swValue,sqrtValue;
      input >> wValue >> swValue >> sqrtValue;
      s.w.push_back(wValue);
      s.sw.push_back(swValue);
      s.sqrt_nentries.push_back(sqrtValue);
    }
    return input;
  }
};


struct SaturationPeak_t
{
  float energy;
  float peakMin;
  float peakMax;
};


//variables for the input TChain
ULong64_t     ChainExtendedTimeTag;                                // extended time tag
ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous event
Int_t        *ChainAdcChannel;
Short_t      *ChainDesktopAdcChannel;                              // input TChain data for desktop digitizers - data is int_16
UShort_t     *ChainVMEadcChannel;                                  // input TChain data for VME digitizers - data is uint_16
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

//branches for the input TChain
TBranch      *bChainExtendedTimeTag;                               // branches for above data
TBranch      *bChainDeltaTimeTag;                                  // branches for above data
TBranch     **bChainAdcChannel;                                    // branches for above data
TBranch     **bChainTimeStamp;
TBranch      *bRealX;                                              // branches for above data
TBranch      *bRealY;                                              // branches for above data
TBranch      *bRealZ;                                              // branches for above data
TBranch      *bsimTaggingCharge;                                              // branches for above data
TBranch      *bsimTaggingTime;                                              // branches for above data
TBranch      *bCrystalsHit;                                        // branches for above data
TBranch      *bNumbOfInteractions;                                 // branches for above data
TBranch      *bTotalCryEnergy;                                     //
