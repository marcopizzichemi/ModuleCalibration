struct correction_t // for each "central" crystal, one of these for each relevant detector
{
  int digitizerChannel;    // number of digi channel for this detector (i.e. N in chN in the TTree files)
  int timingChannel;       // number of time channel for this detector (i.e. N in tN in the TTree files)
  TH1F* CTRhisto;          // histogram of timestamps[timingChannel] - timestamps[taggingCrystalTimingChannel] for this detector
  // TH1F* DelayHisto;        // histogram of delay of this channel wrt central detector. 0 if this is central detector
  TH2F*  CTRhistoVsW;
  // TH2F*  DelayHistoVsW;
};
