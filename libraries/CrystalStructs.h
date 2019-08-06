struct detector_t
{
  int digitizerChannel;
  float saturation;
  float pedestal;
};

struct slice_t
{
  Float_t wmin;
  Float_t wmax;
  Float_t wmean;
  Float_t werr;
  long long int entries;
  std::vector<int> tChannel;
  std::vector<Float_t> averageDelay;
  std::vector<Float_t> averageDeltaT;
  std::vector<Float_t> varianceDeltaT;
  std::vector<long long int> nVarianceDeltaT;
  std::vector<long long int> nDeltaT;
  Float_t **covariance;
  Float_t **inverse_covariance;
  long long int **entries_covariance;
  Float_t **normalized_covariance;

  // TGraphErrors* delay;
  // TF1* delay_line;
  // TH2F* deltaTscatter;
  // TGraphErrors* deltaTgraph;
  // TF1* deltaTline;
  // std::vector<TH1F*> deltaTslice;
  // std::vector<double> wmean;
  // std::vector<double> werr;
  //
  //
  // std::vector<double> dx;
  // std::vector<double> dy;
  // std::vector<double> dex;
  // std::vector<double> dey;
};


struct graphs_t
{
  int timingChannel;
  TGraphErrors* graph;
};

struct correction_graphs_t
{
  int timingChannel;
  bool isMainChannel;
  TGraphErrors* delay;
  TGraphErrors* rms;
};

struct polished_correction_t
{
  int timingChannel;
  double mean;
  double rms;
};

struct Crystal_t
{
  int number;
  int detectorChannel;
  int timingChannel;
  float length;
  std::vector<int> relevantForW;
  std::vector<int> delayTimingChannels;
  TCut *CrystalCut;
  TCut *CrystalCutWithoutCutG;
  TCut *PhotopeakEnergyCut;
  std::vector<TCutG*> cutg;
  TGraph *calibrationGraph;
  TGraph *wz;
  TTreeFormula *Formula;
  TTreeFormula *FormulaAnalysis;
  TTreeFormula *FormulaTagAnalysis;
  int taggingCrystalTimingChannel;
  int taggingCrystalChannel;
  int taggingPosition;
  TCut* taggingPhotopeakCut;
  std::vector<detector_t> detectorSaturation;
  float marginWZgraph;
  float WrangeBinsForTiming;
  float minAcceptedW;
  float maxAcceptedW;
  float wMinSlicing;
  float wMaxSlicing;
  float wStepSlicing;
  std::vector<double> z;
  TGraphErrors* tw_correction;
  TGraphErrors* rms_tw_correction;
  std::vector<TGraphErrors*> delay;
  std::vector<TGraphErrors*> rms_delay;
  TF1 *tw_correction_line;
  TF1 *rms_tw_correction_line;
  std::vector<TF1*> delay_line;
  std::vector<TF1*> rms_delay_line;
  const char* path;
  bool accepted;
  bool polishedCorrection;
  std::vector<int>    tChannelsForPolishedCorrectionMean;
  std::vector<int>    tChannelsForPolishedCorrectionFWHM;
  std::vector<double> meanForPolishedCorrection;
  std::vector<double> fwhmForPolishedCorrection;
  std::vector<polished_correction_t> polished_correction;
  std::vector<slice_t> slice;
  TGraph *** inverse_covariance_element; // matrix of TGraph, one for each element of the inverse covariance element s^{-1}_i,j(w) that is a function of doi...
  TF1 *** inverse_covariance_element_line;
  std::vector<graphs_t> delay_graphs;
  std::vector<graphs_t> rms_graphs;
  std::vector<correction_graphs_t> correction_graphs;
  TH1F* lightCentralHisto;
  TH1F* lightAllHisto;
  TH1F* basicCTRhisto;

  // CTR histograms. Here because it's easier like this, but they shouldn't belong...
  TH1F *simpleCTR;
  TH1F *centralCTR;
  TH1F *allCTR;
  TH1F *poliCorrCTR;
  TH1F *likeCTR;
  TH1F *hybridCTR;

  TH1F *simpleCTR_norm;
  TH1F *centralCTR_norm;
  TH1F *allCTR_norm;
  TH1F *poliCorrCTR_norm;
  TH1F *hybridCTR_norm;
  TH1F *likeCTR_norm;

  TH2F *ctrVSw;
  TH1D* aSlices;
  TGraph *crtVSw_gr;

  TCut* TriggerChannelCut ;
  TCut* broadCut          ;
  TCut* CutNoise          ;
  TCut* PhotopeakEnergyCutnergyCut;

  //CTR std::vectors for unbinned calculations. obsolete
  // std::vector<double> vSimple;
  // std::vector<double> vCentral;
  // std::vector<double> vAll;
  // std::vector<double> vPoli;
  // std::vector<double> vLike;
  // std::vector<double> vhybrid;
};



bool compareByNumber(const Crystal_t &a,const Crystal_t  &b)
{
  return a.number < b.number;
}
