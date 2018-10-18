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


struct multi_channel_t
{
  int detectorIndex;
  std::string string;
};



//------------------------------------------//
// CREATION OF STRINGS FOR PLOTTING         //
//------------------------------------------//
// in order to compute all th plots on the fly, it is necessary to put together the proper strings to pass to the draw commands
// this has to include also the saturation correction, the pedestal position and the noise
// furthermore, we need to prepare some collection of elements for the channels that are relevant for computation of energy, u-v coordinates and w coordinates
// it all starts from the detector array, produced and filled by the InputFile.cc part of the program. Declaration is in ./include/Detector.h
// this vector of struct holds an entry for each mppc. each entry has data on the mppc. In particular, the digitizer (aka charge) channel number and
// the timing channel number, a vector with the digitizer channel numbers of all the neighbour channels, and the info on saturation, pedestal and noise
// The steps are the following (more explainations written just before each step below):

// 1. Get the index of this MPPC in the detector array                                                                 ---> thisChannelID
// 2. Get the indexes of the neighbour MPPCs in the detector array                                                     ---> neighbourChID[]
// 3. Translate indexes into strings, taking also into account saturation (if requested), pedestals and noise          ---> thisChannel
//                                                                                                                     ---> neighbourChannels[]
//                                                                                                                     ---> allChannels[]
// 4. Generate the TriggerChannel condition                                                                            ---> TriggerChannel
// 5. Get the relevant channels for this MPPC                                                                          ---> relevantForUV[]
//                                                                                                                     ---> relevantForW[]
//                                                                                                                     ---> relevantForEnRes[]
// 6. Compose the FloodX and FloodY strings                                                                            ---> FloodX
//                                                                                                                     ---> FloodY
// 7. Compose the FloodZ string                                                                                        ---> FloodZ
// 8. Compose an additional trigger condition, cutting the noise                                                       ---> relevantForNoiseCut[]
//                                                                                                                     ---> noiseCut
//                                                                                                                     ---> summed to CutTrigger


// 1. Get the index of this MPPC in the detector array
//    Simply, run on all the entries of detector array and get the one where detector[i].digitizerChannel == mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetDigitizerChannel()
//    ---> thisChannelID

// 2. Get the indexes of the neighbour MPPCs in the detector array
//    Same logic as in 1., run on all entries on detector[] and store in vector neighbourChID
//    ---> neighbourChID[]

// 3. Translate indexes into strings, taking also into account saturation (if requested), pedestals and noise
//    Two info need to be stored together:
//           a) the channel ID (i.e. the digitizer channel)
//           b) the complete string that will be passed to the Draw command
//    So it is necessary to create a multi_channel_t struct, with these two entries (see definition above)
//    Then, 3 arrays of this multi_channel_t elements are created:
//    ---> thisChannel        = holds the string and digitizerChannel for this MPPC, so it's not a vector but one struct
//    ---> neighbourChannels  = vector holding the string and digitizerChannel for each of the neighbour channels to this MPPC
//    ---> allChannels        = vector holding the string and digitizerChannel for all channels in this MPPC array
//    The 3 arrays will be used in step 4. and 5.

// 4. Generate the TriggerChannel condition
//    the idea is:        TriggerChannel     = "max(ch0,max(ch1,...max(chN))) == chI", for I == mppc channel and 0...N are all the channels in this module
//    in this case the MPPC channels involved are all the channels, so the allChannels array

// 5. Get the relevant channels for this MPPC
//    Relevant channels are split in 3
//    5.1. relevant for u-v coordinates
//    5.2. relevant for w coordinate
//    5.3. relevant for energy resolution (or better for sum spectrum)
//    The 3 conditions are controlled by the corresponding config file keys
//    5.1. usingAllChannels  -> 0 = just neighbour channels
//                         -> 1 = all channels in the module
//    5.2. wAllChannels      -> 0 = just neighbour channels
//                         -> 1 = all channels in the module
//    5.3. usingAllChannelsForEnergySpectra  -> 0 = just neighbour channels
//                                         -> 1 = all channels in the module
//    Remember that the string for the channels, already including the saturation correction if any, have been created in step 3.
//    multi_channel_t thisChannel;
//    std::vector<multi_channel_t> neighbourChannels;
//    std::vector<multi_channel_t> allChannels;

// 6. Compose the FloodX and FloodY strings
//    Standard Anger logic of charges (corrected by saturation and pedestal) with mppc position, for the channels defined as relevant for u-v

// 7. Compose the FloodZ string
//    Ratio of charge of this mppc versus charges of the relevant mppcs for w computation. all charges are corrected by saturation and pedestal

// 8. Compose an additional trigger condition, cutting the noise
//    the idea is to ignore events when the charge recorded by one of the detectors involved is == 0
//    ((ch12 > 0)+ (ch11>0) + (ch0>0)+ (ch4>0) + (ch7>0) +(ch10>0) +(ch13>0) +(ch14>0)+ (ch15>0)) > 8
//    the channels involved are defined as the biggest collection between relevantForUV and relevantForW
