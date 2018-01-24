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
  float sigma;
  float fractionLow;
  float fractionHigh;
  // float peakMax;
};

struct Peaks_t
{
  float x;
  float y;
  float sx;
  float energy;
  TF1* fit;
};

bool compareByX(const Peaks_t &a,const Peaks_t  &b)
{
  return a.x < b.x;
}

bool compareByEnergy(const SaturationPeak_t &a,const SaturationPeak_t  &b)
{
  return a.energy < b.energy;
}
