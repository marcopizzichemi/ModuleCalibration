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
  std::vector<int> relevantForW;
  std::vector<int> delayTimingChannels;
  TCut *CrystalCut;
  TCut *CrystalCutWithoutCutG;
  TCut *PhotopeakEnergyCut;
  std::vector<TCutG*> cutg;
  TGraph *calibrationGraph;
  TGraph *wz;
  TH1F *simpleCTR;
  TH1F *centralCTR;
  TH1F *allCTR;
  TH1F *poliCorrCTR;
  TH1F *likeCTR;
  TH1F *hybridCTR;
  std::vector<double> vSimple;
  std::vector<double> vCentral;
  std::vector<double> vAll;
  std::vector<double> vPoli;
  std::vector<double> vLike;
  std::vector<double> vhybrid;
  TH1F *simpleCTR_norm;
  TH1F *centralCTR_norm;
  TH1F *allCTR_norm;
  TH1F *poliCorrCTR_norm;
  TH1F *hybridCTR_norm;
  TTreeFormula *Formula;
  TTreeFormula *FormulaAnalysis;
  TH1F *likeCTR_norm;

  float marginWZgraph;
  float WrangeBinsForTiming;

  float minAcceptedW;  // events w min e max for ->Eval operations
  float maxAcceptedW;
  float wMinSlicing;   // limits of w in scatter plot slicing
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
  // std::vector<ctr_aligned_t> ctr_aligned;


  // Float_t *variance;
  // long long int *entries_variance;
  // std::vector<Float_t**> covariance;
  // std::vector<long long int**> entries_covariance;
  // std::vector<Float_t**> inverse_covariance;
  //
  // TH2F** deltaTscatter;

  TH1F* lightCentralHisto;
  TH1F* lightAllHisto;
  TH1F* basicCTRhisto;

  // TCanvas
};

struct detector_t
{
  int digitizerChannel;
  float saturation;
  float pedestal;
};
