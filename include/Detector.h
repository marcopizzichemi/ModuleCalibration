struct detector_t
{
  int digitizerChannel;
  int timingChannel;
  std::string label;
  int i;
  int j;
  float saturation;
  int plotPosition;
  float xPosition;
  float yPosition;
  float pedestal;
  float noise;
  float gain;
  int OnForDOI;
  std::vector<int> neighbourChannels;
  bool OnForModular;
  //     bool operator<(const masks_t& rhs) const { meanx < rhs.meanx; }
};
