// compile with
// g++ -o ../build/new_correlation new_correlation.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

// small program to extract timing calibration and data

#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"
#include "TEllipse.h"
#include "TFormula.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "TCutG.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TGraphDelaunay.h"
#include "TVector.h"
#include "TNamed.h"
#include "TPaveLabel.h"
#include "THStack.h"
#include "TFitResult.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <getopt.h>
#include <algorithm>    // std::sort

#include <sys/types.h>
#include <dirent.h>

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

// typedef std::vector<std::string> stringvec;
// list files in directory
// taken from
// http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
void read_directory(const std::string& name, std::vector<std::string> &v)
{
  DIR* dirp = opendir(name.c_str());
  struct dirent * dp;
  while ((dp = readdir(dirp)) != NULL) {
    v.push_back(dp->d_name);
  }
  closedir(dirp);
}


// void getCofactor(Double_t** mat, Double_t** temp, int p, int q, int n)
// {
//     int i = 0, j = 0;
//
//     // Looping for each element of the matrix
//     for (int row = 0; row < n; row++)
//     {
//         for (int col = 0; col < n; col++)
//         {
//             //  Copying into temporary matrix only those element
//             //  which are not in given row and column
//             if (row != p && col != q)
//             {
//                 temp[i][j++] = mat[row][col];
//
//                 // Row is filled, so increase row index and
//                 // reset col index
//                 if (j == n - 1)
//                 {
//                     j = 0;
//                     i++;
//                 }
//             }
//         }
//     }
// }
//
// /* Recursive function for finding determinant of matrix.
//    n is current dimension of mat[][]. */
// double determinantOfMatrix(Double_t mat[N][N], int n)
// {
//     double D = 0; // Initialize result
//
//     //  Base case : if matrix contains single element
//     if (n == 1)
//         return mat[0][0];
//
//     Double_t temp[N][N]; // To store cofactors
//
//     int sign = 1;  // To store sign multiplier
//
//      // Iterate for each element of first row
//     for (int f = 0; f < n; f++)
//     {
//         // Getting Cofactor of mat[0][f]
//         getCofactor(mat, temp, 0, f, n);
//         D += sign * mat[0][f] * determinantOfMatrix(temp, n - 1);
//
//         // terms are to be added with alternate sign
//         sign = -sign;
//     }
//
//     return D;
// }


gsl_matrix * invert_a_matrix(gsl_matrix *matrix,size_t size)
{
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;
}

void print_mat_contents(gsl_matrix *matrix,size_t size)
{
    size_t i, j;
    double element;

    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            element = gsl_matrix_get(matrix, i, j);
            printf("%.5e\t", element);
        }
        printf("\n");
    }
}


struct LikelihoodElement_t
{
  int t_channel;
  TH2F* mu_w;
};

struct Crystal_t
{
  int number;
  int detectorChannel;
  int timingChannel;
  std::vector<int> relevantForW;
  std::vector<int> delayTimingChannels;
  // std::vector<int> likelihoodTimingChannels;
  std::vector<LikelihoodElement_t> likelihood_mu_w;
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
  std::vector<double> vSimple;
  std::vector<double> vCentral;
  std::vector<double> vAll;
  std::vector<double> vPoli;
  TH1F *simpleCTR_norm;
  TH1F *centralCTR_norm;
  TH1F *allCTR_norm;
  TH1F *poliCorrCTR_norm;
  TH1F *likelihood_CTR;
  TTreeFormula *Formula;
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
  std::vector<int>    tChannelsForPolishedCorrection;
  std::vector<double> meanForPolishedCorrection;
  std::vector<double> fwhmForPolishedCorrection;
  // double *** corr_mat;
  TH2F*** correlation_matrix_element;
  TH1F**** h_correlation_matrix_element_w; //[iW][iLike][jLike]
  std::vector<double> w_limit_low;
  std::vector<double> w_limit_up;
  double** mu_i_w;
  TH1F*** mu_histo;
  long int** mu_i_w_counter;
  gsl_matrix **corr_mat; // covariance matrix gsl[iW][iLike][jLike]
  gsl_matrix **norm_corr_mat; // normalized covariance matrix gsl[iW][iLike][jLike]
  gsl_matrix **inverse;
  TH2F* scatter_likeT_w;

  long int*** corr_mat_counter;
  double *normalization_factor;

  // TCanvas
};

struct detector_t
{
  int digitizerChannel;
  float saturation;
  float pedestal;
};


/*** find width at half max ***/
// float ComputeFWHM(TH1F* histo)
// {
//    float max = histo->GetMaximum();
//    float halfMax = max / 2.0;
//    int binMin = histo->FindFirstBinAbove(halfMax);
//    int binMax = histo->FindLastBinAbove(halfMax);
//    float down = histo->GetBinCenter(binMin);
//    float up   = histo->GetBinCenter(binMax);
//    float ret = up -down;
//    return ret;
// }

void extractWithGaussAndExp(TH1F* histo,double fitPercMin,double fitPercMax, int divs, double tagFwhm, double* res)
{
  TCanvas *cTemp  = new TCanvas("temp","temp");
  TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  histo->Fit(gaussDummy,"QN");

  double f1min = histo->GetXaxis()->GetXmin();
  double f1max = histo->GetXaxis()->GetXmax();


  TF1* f1  = new TF1("f1","[0]/sqrt(2)*exp([2]^2/2/[3]^2-(x-[1])/[3])*(1-TMath::Erf(([1]-x+[2]^2/[3])/(sqrt(2*[2]^2))))");
  f1->SetLineColor(kBlack);
  f1->SetParName(0,"N");
  f1->SetParName(1,"Mean");
  f1->SetParName(2,"Sigma");
  f1->SetParName(3,"tau");
  // f1->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  f1->SetParameter(0,gaussDummy->GetParameter(0));
  f1->SetParameter(1,gaussDummy->GetParameter(1));
  f1->SetParameter(2,gaussDummy->GetParameter(2));
  f1->SetParameter(3,5e-9);
  double fitMin = gaussDummy->GetParameter(1) - fitPercMin*(gaussDummy->GetParameter(2));
  double fitMax = gaussDummy->GetParameter(1) + fitPercMax*(gaussDummy->GetParameter(2));
  if(fitMin < f1min)
  {
    fitMin = f1min;
  }
  if(fitMax > f1max)
  {
    fitMax = f1max;
  }

  histo->Fit(f1,"","",fitMin,fitMax);
  std::cout << "-----------------------------------------------------------------------" << std::endl;
  std::cout << histo->GetName() << std::endl;
  std::cout << f1->GetMaximum(fitMin,fitMax)    << "\t"
  // << fitMin              << "\t"
  // << fitMax              << "\t"
  << f1->GetParameter(0) << "\t"
  << f1->GetParameter(1) << "\t"
  << f1->GetParameter(2) << "\t"
  << f1->GetParameter(3) << "\t"
  << std::endl;
  double min,max,min10,max10;
  // int divs = 3000;
  double step = (fitMin-fitMax)/divs;
  double funcMax = f1->GetMaximum(fitMin,fitMax);
  for(int i = 0 ; i < divs ; i++)
  {
    if( (f1->Eval(f1min + i*step) < funcMax/2.0) && (f1->Eval(f1min + (i+1)*step) > funcMax/2.0) )
    {
      min = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) > funcMax/2.0) && (f1->Eval(f1min + (i+1)*step) < funcMax/2.0) )
    {
      max = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) < funcMax/10.0) && (f1->Eval(f1min + (i+1)*step) > funcMax/10.0) )
    {
      min10 = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) > funcMax/10.0) && (f1->Eval(f1min + (i+1)*step) < funcMax/10.0) )
    {
      max10 = f1min + (i+0.5)*step;
    }
  }
  // res[0] = f1->GetParameter(1);  // res[0] is mean
  // res[1] = max-min;              // res[1] is FWHM
  res[0] = sqrt(2)*sqrt(pow((max-min),2)-pow(tagFwhm,2));
  res[1] = sqrt(2)*sqrt(pow((max10-min10),2)-pow((tagFwhm/2.355)*4.29,2));
  std::cout <<"----> " << res[0] << " " << res[1]<< std::endl;
  delete cTemp;

}

void extractCTR(TH1D* histo,double fitPercMin,double fitPercMax, int divs, double* res)
{
  // TF1 *gexp;
  // TF1 *cb;

  // preliminary gauss fit
  TCanvas *cTemp  = new TCanvas("temp","temp");
  TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  gaussDummy->SetLineColor(kRed);
  histo->Fit(gaussDummy,"QN");
  double f1min = histo->GetXaxis()->GetXmin();
  double f1max = histo->GetXaxis()->GetXmax();
  double fitMin = gaussDummy->GetParameter(1) - fitPercMin*(gaussDummy->GetParameter(2));
  double fitMax = gaussDummy->GetParameter(1) + fitPercMax*(gaussDummy->GetParameter(2));
  if(fitMin < f1min)
  {
    fitMin = f1min;
  }
  if(fitMax > f1max)
  {
    fitMax = f1max;
  }

  //fit with crystalball
  TF1 *cb  = new TF1("cb","crystalball",f1min,f1max);
  cb->SetLineColor(kBlue);
  cb->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  TFitResultPtr rCb = histo->Fit(cb,"QNS","",fitMin,fitMax);

  //fit with gauss + exp
  TF1* gexp  = new TF1("gexp","[0]/sqrt(2)*exp([2]^2/2/[3]^2-(x-[1])/[3])*(1-TMath::Erf(([1]-x+[2]^2/[3])/(sqrt(2*[2]^2))))",f1min,f1max);
  gexp->SetLineColor(kGreen);
  gexp->SetParName(0,"N");
  gexp->SetParName(1,"Mean");
  gexp->SetParName(2,"Sigma");
  gexp->SetParName(3,"tau");
  // f1->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  gexp->SetParameter(0,gaussDummy->GetParameter(0));
  gexp->SetParameter(1,gaussDummy->GetParameter(1));
  gexp->SetParameter(2,gaussDummy->GetParameter(2));
  gexp->SetParameter(3,gaussDummy->GetParameter(2)); // ROOT really needs all parameters initialized, and a "good" guess for tau is the sigma of the previous fit...
  TFitResultPtr rGexp = histo->Fit(gexp,"QNS","",fitMin,fitMax);

  Int_t fitStatusCb = rCb;
  Int_t fitStatusGexp = rGexp;

  double chi2gexp;
  double chi2cb;

  if(!fitStatusGexp)
  {
    chi2gexp = rGexp->Chi2();
  }
  if(!fitStatusCb)
  {
    chi2cb   = rCb->Chi2();
  }


  //set function to measure
  TF1 *f1;
  if(fitStatusGexp && fitStatusCb) // both fit didn't work, just set everything to 0
  {
    res[0] = 0;
    res[1] = 0;
    res[2] = 0;
    res[3] = 0;
  }
  else
  {
    if(fitStatusGexp && !fitStatusCb) // only cb worked
    {
      f1 = cb;
      delete gexp;

    }
    else if(!fitStatusGexp && fitStatusCb) // only gexp worked
    {
      f1 = gexp;
      delete cb;
    }
    else // both worked
    {
      if(chi2gexp > chi2cb)
      {
        f1 = cb;
        delete gexp;
      }
      else
      {
        f1 = gexp;
        delete cb;
      }
    }

    f1->SetLineColor(kRed);
    // double f1min = histo->GetXaxis()->GetXmin();
    // double f1max = histo->GetXaxis()->GetXmax();
    // dummy re-fit just to draw the function and save it in the histogram...
    histo->Fit(f1,"Q","",fitMin,fitMax);

    // new version, get the histogram used to draw the f1 function and directly get mean and rms from it
    // idiotic way but easier and faster. thanks ROOT


    // double min,max;
    // double step = (f1max-f1min)/divs;
    // double funcMax = f1->GetMaximum(fitMin,fitMax);
    // // is [0] the max of the function???
    // for(int i = 0 ; i < divs ; i++)
    // {
    //   if( (f1->Eval(f1min + i*step) < funcMax/2.0) && (f1->Eval(f1min + (i+1)*step) > funcMax/2.0) )
    //   {
    //     min = f1min + (i+0.5)*step;
    //   }
    //   if( (f1->Eval(f1min + i*step) > funcMax/2.0) && (f1->Eval(f1min + (i+1)*step) < funcMax/2.0) )
    //   {
    //     max = f1min + (i+0.5)*step;
    //   }
    // }
    res[0] = f1->GetHistogram()->GetMean();      // res[0] is mean
    res[1] = f1->GetHistogram()->GetRMS();       // res[1] is now RMS
    res[2] = f1->GetHistogram()->GetMeanError(); // res[2] is mean error
    res[3] = f1->GetHistogram()->GetRMSError();  // res[3] is RMS error

    // //DEBUG
    // std::cout << f1->GetHistogram()->GetMean() << " "
    //           << f1->GetParameter(1) << " "
    //           << f1->GetHistogram()->GetRMS() << " "
    //           << max - min
    //           << std::endl;
  }
}



//**** per std::vector -- non binnata
double FindSmallestInterval(double& mean,
  double& meanErr,
  double& min,
  double& max,
  std::vector<double>& vals,
  const double& fraction,
  const bool& verbosity)
  {
    if( verbosity )
    std::cout << ">>>>>> FindSmallestInterval" << std::endl;


    std::sort(vals.begin(),vals.end());

    unsigned int nPoints = vals.size();
    unsigned int maxPoints = (unsigned int)(fraction * nPoints);

    unsigned int minPoint = 0;
    unsigned int maxPoint = 0;
    double delta = 999999.;
    for(unsigned int point = 0; point < nPoints-maxPoints; ++point)
    {
      double tmpMin = vals.at(point);
      double tmpMax = vals.at(point+maxPoints-1);
      if( tmpMax-tmpMin < delta )
      {
        delta = tmpMax - tmpMin;
        min = tmpMin;
        max = tmpMax;
        minPoint = point;
        maxPoint = point + maxPoints - 1;
      }
    }
    return delta;
  }


  /*** find effective sigma ***/
  void FindSmallestInterval(double* retValues, TH1F* histo, const float&
    fraction, const bool& verbosity, double tagFwhm)
    {
      float ret[4];
      float integralMax = fraction * histo->Integral();

      int N = histo -> GetNbinsX();
      std::vector<float> binCenters(N);
      std::vector<float> binContents(N);
      std::vector<float> binIntegrals(N);
      for(int bin1 = 0; bin1 < N; ++bin1)
      {
        binCenters[bin1] = histo->GetBinCenter(bin1+1);
        binContents[bin1] = histo->GetBinContent(bin1+1);

        for(int bin2 = 0; bin2 <= bin1; ++bin2)
        binIntegrals[bin1] += binContents[bin2];
      }

      float min = 0.;
      float max = 0.;
      float delta = 999999.;
      for(int bin1 = 0; bin1 < N; ++bin1)
      {
        for(int bin2 = bin1+1; bin2 < N; ++bin2)
        {
          if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;

          float tmpMin = histo -> GetBinCenter(bin1);
          float tmpMax = histo -> GetBinCenter(bin2);

          if( (tmpMax-tmpMin) < delta )
          {
            delta = (tmpMax - tmpMin);
            min = tmpMin;
            max = tmpMax;
          }

          break;
        }
      }

      TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
      for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
      {
        if( smallHisto->GetBinCenter(bin) < min )
        smallHisto -> SetBinContent(bin,0);

        if( smallHisto->GetBinCenter(bin) > max )
        smallHisto -> SetBinContent(bin,0);
      }
      smallHisto -> SetFillColor(kYellow);

      float mean = smallHisto -> GetMean();
      float meanErr = smallHisto -> GetMeanError();

      ret[0] = mean;
      ret[1] = meanErr;
      ret[2] = min;
      ret[3] = max;

      //mean is the smallest interval containing the 68% (fraction) of data. this would be from -1 sigma to +1 sigma, so 2 sigmas. therefore we get teh "fwhm" of this distro by
      double fwhm = 2.355* ((max-min) / 2.0);

      retValues[0] = sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
      retValues[1] = 0;
    }


    bool compareByNumber(const Crystal_t &a,const Crystal_t  &b)
    {
      return a.number < b.number;
    }

    void usage()
    {
      std::cout << "\t\t" << "[-i|--input] <file_prefix>  [-o|--output] <output.root> [-c|--calibration] calibration.root [--coincidence] coincidence.root [OPTIONS]" << std::endl
      << "\t\t" << "<file_prefix>                                      - prefix of TTree files to analyze"   << std::endl
      << "\t\t" << "<output.root>                                      - output file name"   << std::endl
      << "\t\t" << "<calibration.root>                                 - calibration file " << std::endl
      // << "\t\t" << "<coincidence.root>                                 - time calibration file " << std::endl
      // << "\t\t" << "--simulation                                       - the datast is from a simulation (therefore the tagging photopeak is ignored)" << std::endl
      // << "\t\t" << "--length <value>                                   - crystal length in mm, default = 15.0"  << std::endl
      // << "\t\t" << "--doiFraction <value>                              - fraction of DOI length towards which the time stamps are corrected (from 0 to 1)"  << std::endl
      // << "\t\t" << "                                                   - 0 = front of the crystal (DOI close to detector) "  << std::endl
      // << "\t\t" << "                                                   - 1 = back of the crystal (DOI far from detector) "  << std::endl
      << "\t\t" << "--tagFwhm <value>                                  - FWHM timing resolution of reference board, in sec - default = 70e-12"  << std::endl
      // << "\t\t" << "--rmsLow <value>                                   - lower bound of CTR fit -> mean - rmsLow*mean - default = 1.75"  << std::endl
      // << "\t\t" << "--rmsHigh <value>                                  - upper bound of CTR fit -> mean + rmsHigh*mean - default = 1.75"  << std::endl
      << "\t\t" << "--histoMin <value>                                 - lower limit of CTR spectra, in sec - default = -15e-9"  << std::endl
      << "\t\t" << "--histoMax <value>                                 - upper limit of CTR spectra, in sec - default = 15e-9"  << std::endl
      << "\t\t" << "--histoBins <value>                                - n of bins for CTR spectra - default = 500"  << std::endl
      << "\t\t" << "--corrMin <value>                                 - lower limit of y axis in correlation element spectra, in sec^2 - default = -5e-20"  << std::endl
      << "\t\t" << "--corrMax <value>                                 - upperlimit of y axis in correlation element spectra, in sec^2 - default = 5e-20"  << std::endl
      << "\t\t" << "--corrBins <value>                                - n of bins for y axis in correlation element spectra - default = 1000"  << std::endl
      // << "\t\t" << "--smooth <value>                                   - n of iteration in CTR histograms smoothing - default = 0 (no smoothing)"  << std::endl
      << "\t\t" << "--fitPercMin <value>                               - time fit min is set to ((gauss fit mean) - fitPercMin*(gauss fit sigma))  - default = 5"  << std::endl
      << "\t\t" << "--fitPercMax <value>                               - time fit max is set to ((gauus fit mean) - fitPercMax*(gauss fit sigma))  - default = 6" << std::endl
      << "\t\t" << "--divs <value>                                     - n of divisions when looking for FWHM - default = 10000"  << std::endl
      << "\t\t" << "--bins <value>                                     - n of bins in summary CTR histograms - deafult 40"  << std::endl
      << "\t\t" << "--func <value>                                     - function for fitting (default = 0)"  << std::endl
      << "\t\t" << "                                                   - 0 = crystalball "  << std::endl
      << "\t\t" << "                                                   - 1 = gauss+exp "  << std::endl
      // << "\t\t" << "--unbinned                                         - use also the unbinned method to calculate CTR - default = 0 (false)"  << std::endl
      // << "\t\t" << "--fitCorrection                                    - use line fit to perform correction   - default = 0 (false)"  << std::endl
      << "\t\t" << "--exclude-channels                                 - channels to exclude from time correction, comma separated - default = "" "  << std::endl
      << "\t\t" << "--sections                                 - number of doi slices - default = 10 "  << std::endl
      << "\t\t" << std::endl;
    }

    int main (int argc, char** argv)
    {
      if(argc < 2)
      {
        std::cout << argv[0] << std::endl;
        usage();
        return 1;
      }
      std::string inputFileName = "";
      std::string outputFileName = "";
      std::string calibrationFileName = "";
      // std::string coincidenceCalibrationFileName = "";

      std::string exclude_channels = "";
      bool exclude = false;
      int sections = 10;
      // bool simulation = false;
      // Float_t length = 15.0; //mm
      // Float_t doiFraction = 0.5;
      Float_t tagFwhm = 70.0e-12; //s
      // Float_t rmsLow = 1.75;
      // Float_t rmsHigh = 1.75;
      Float_t histoMin = -15e-9;//s
      Float_t histoMax = 15e-9;//s
      Float_t fitPercMin = 5;
      Float_t fitPercMax = 6;
      int divs       = 10000;
      int histoBins = 500;
      int smooth = 0; //
      int bins = 40;
      double minCTR = 100;
      double maxCTR = 500;
      int func = 0;
      bool unbinned = false;
      bool fitCorrection = false;
      Float_t corrMin = -5e-20;//s
      Float_t corrMax = 5e-20;//s
      int corrBins = 1000;
      // parse arguments
      static struct option longOptions[] =
      {
        { "input", required_argument, 0, 0 },
        { "output", required_argument, 0, 0 },
        { "calibration", required_argument, 0, 0 },
        // { "simulation", no_argument, 0, 0 },
        // { "length", required_argument, 0, 0 },
        // { "doiFraction", required_argument, 0, 0 },
        // { "coincidence", required_argument, 0, 0 },
        { "tagFwhm", required_argument, 0, 0 },
        // { "rmsLow", required_argument, 0, 0 },
        // { "rmsHigh", required_argument, 0, 0 },
        { "histoMin", required_argument, 0, 0 },
        { "histoMax", required_argument, 0, 0 },
        { "histoBins", required_argument, 0, 0 },
        // { "smooth", required_argument, 0, 0 },
        { "fitPercMin", required_argument, 0, 0 },
        { "fitPercMax", required_argument, 0, 0 },
        { "divs", required_argument, 0, 0 },
        { "bins", required_argument, 0, 0 },
        { "func", required_argument, 0, 0 },
        // { "unbinned", no_argument, 0, 0 },
        // { "fitCorrection", no_argument, 0, 0 },
        { "exclude-channels", required_argument, 0, 0 },
        { "corrMin", required_argument, 0, 0 },
        { "corrMax", required_argument, 0, 0 },
        { "corrBins", required_argument, 0, 0 },
        { "sections", required_argument, 0, 0 },
        { NULL, 0, 0, 0 }
      };

      while(1) {
        int optionIndex = 0;
        int c = getopt_long(argc, argv, "i:o:c:", longOptions, &optionIndex);
        if (c == -1) {
          break;
        }
        if (c == 'i'){
          inputFileName = (char *)optarg;
        }
        else if (c == 'o'){
          outputFileName = (char *)optarg;
        }
        else if (c == 'c'){
          calibrationFileName = (char *)optarg;
        }
        else if (c == 0 && optionIndex == 0){
          inputFileName = (char *)optarg;
        }
        else if (c == 0 && optionIndex == 1){
          outputFileName = (char *)optarg;
        }
        else if (c == 0 && optionIndex == 2){
          calibrationFileName = (char *)optarg;
        }
        // else if (c == 0 && optionIndex == 3){
        //   std::cout << "Dataset from simulation " << std::endl;
        //   simulation = true;
        // }
        // else if (c == 0 && optionIndex == 4){
        //   length = atof((char *)optarg);;
        // }
        // else if (c == 0 && optionIndex == 5){
        //   doiFraction = atof((char *)optarg);;
        // }
        // else if (c == 0 && optionIndex == 6){
        //   coincidenceCalibrationFileName = (char *)optarg;
        // }
        else if (c == 0 && optionIndex == 3){
          tagFwhm = atof((char *)optarg);
        }
        // else if (c == 0 && optionIndex == 7){
        //   rmsLow = atof((char *)optarg);
        // }
        // else if (c == 0 && optionIndex == 8){
        //   rmsHigh = atof((char *)optarg);
        // }
        else if (c == 0 && optionIndex == 4){
          histoMin = atof((char *)optarg);
        }
        else if (c == 0 && optionIndex == 5){
          histoMax = atof((char *)optarg);
        }
        else if (c == 0 && optionIndex == 6){
          histoBins = atoi((char *)optarg);
        }
        // else if (c == 0 && optionIndex == 12){
        //   smooth = atoi((char *)optarg);
        // }
        else if (c == 0 && optionIndex == 7){
          fitPercMin = atof((char *)optarg);
        }
        else if (c == 0 && optionIndex == 8){
          fitPercMax = atof((char *)optarg);
        }
        else if (c == 0 && optionIndex == 9){
          divs = atoi((char *)optarg);
        }
        else if (c == 0 && optionIndex == 10){
          bins = atoi((char *)optarg);
        }
        else if (c == 0 && optionIndex == 11){
          func = atoi((char *)optarg);
        }
        // else if (c == 0 && optionIndex == 12){
        //   unbinned = true;
        // }
        // else if (c == 0 && optionIndex == 13){
        //   fitCorrection = true;
        // }
        else if (c == 0 && optionIndex == 12){
          //take string
          // tokenize by ","
          // store channel numbers
          exclude = true;
          exclude_channels = (char *)optarg;

        }
        else if (c == 0 && optionIndex == 13){
          corrMin = atof((char *)optarg);
        }
        else if (c == 0 && optionIndex == 14){
          corrMax = atof((char *)optarg);
        }
        else if (c == 0 && optionIndex == 15){
          corrBins = atoi((char *)optarg);
        }
        else if (c == 0 && optionIndex == 16){
          sections = atoi((char *)optarg);
        }
        else {
          std::cout	<< "Usage: " << argv[0] << std::endl;
          usage();
          return 1;
        }
      }


      std::vector<int> forbidden_channels;
      if(exclude)
      {
        // std::vector<int> vect;
        std::stringstream ss(exclude_channels);
        int i;
        while (ss >> i)
        {
          forbidden_channels.push_back(i);
          if (ss.peek() == ',')
          ss.ignore();

        }
        std::cout << "Channels excluded from time correction (for depolished): " << std::endl;
        for (i=0; i< forbidden_channels.size(); i++)
        std::cout << forbidden_channels.at(i)<<std::endl;
      }

      // std::cout << "Chosen (length * doiFraction) = " << length * doiFraction << std::endl;
      // if(fitCorrection)
      // {
      //   std::cout << "Using linear fits to perform time correction" << std::endl;
      // }

      // read file in dir

      size_t foundDir = inputFileName.find_last_of("/");
      if(foundDir == std::string::npos)
      {
        foundDir = 0;
      }
      std::string dirName = inputFileName.substr(0,foundDir+1);
      std::string fileName = inputFileName.substr(foundDir+1,inputFileName.size()-(foundDir+1));
      std::cout << "Input Directory    = " << dirName  << std::endl;
      std::cout << "Input files prefix = " << fileName << std::endl;
      std::vector<std::string> v;
      read_directory(dirName.c_str(), v);
      // std::copy(v.begin(), v.end(),std::ostream_iterator<std::string>(std::cout, "\n"));
      // extract files with correct prefix
      std::vector<std::string> listInputFiles;

      for(unsigned int i = 0 ; i < v.size() ; i++)
      {
        if(!v[i].compare(0,fileName.size(),fileName))
        {
          std::string sumString = dirName + v[i];
          listInputFiles.push_back(sumString.c_str());
        }
      }


      //prepare output text file
      // std::string textFileName = outputFileName.substr(0,outputFileName.size()-5);
      // textFileName += ".txt";
      // std::cout << textFileName << std::endl;

      // std::ofstream textfile;
      // textfile.open (textFileName.c_str(),std::ofstream::out);




      // TFile *treeFile = new TFile(inputFileName.c_str());
      // TTree* tree = (TTree*) treeFile->Get("adc");

      //----------------------------------------------------------//
      //  Get TChain of input files                               //
      //----------------------------------------------------------//
      TChain* tree = new TChain("adc");  // create the input tchain and the analysis ttree
      for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
      {
        std::cout << "Adding file " << listInputFiles[i] << std::endl;
        tree->Add(listInputFiles[i].c_str());
      }

      // Add several TTreeFormula to the list;
      TList formulas;

      //variables for the input TChain
      // ULong64_t     ChainExtendedTimeTag;                                // extended time tag
      // ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous event
      // Int_t        *ChainAdcChannel;
      // Short_t      *ChainDesktopAdcChannel;                              // input TChain data for desktop digitizers - data is int_16
      // UShort_t     *ChainVMEadcChannel;                                  // input TChain data for VME digitizers - data is uint_16
      // Float_t      *ChainTimeStamp;
      // Float_t      *TDCBinning;
      // // Short_t      *ChainPetirocChannel;                                 //FIXME temporary data type of petiroc charge input - ask
      // Float_t       RealX;                                               // "real" gamma interaction positions (from simulation data)
      // Float_t       RealY;                                               // "real" gamma interaction positions (from simulation data)
      // Float_t       RealZ;                                               // "real" gamma interaction positions (from simulation data)
      // Float_t       simTaggingCharge;
      // Float_t       simTaggingTime;
      // Short_t       CrystalsHit;                                         // "real" number of crystals hit in the event (from simulation data)
      // Short_t       NumbOfInteractions;                                  // "real" number of interaction (energy depositions) in the event (from simulation data)
      //
      // //branches for the input TChain
      // TBranch      *bChainExtendedTimeTag;                               // branches for above data
      // TBranch      *bChainDeltaTimeTag;                                  // branches for above data
      // TBranch     **bChainAdcChannel;                                    // branches for above data
      // TBranch     **bChainTimeStamp;
      // TBranch      *bRealX;                                              // branches for above data
      // TBranch      *bRealY;                                              // branches for above data
      // TBranch      *bRealZ;                                              // branches for above data
      // TBranch      *bsimTaggingCharge;                                              // branches for above data
      // TBranch      *bsimTaggingTime;                                              // branches for above data
      // TBranch      *bCrystalsHit;                                        // branches for above data
      // TBranch      *bNumbOfInteractions;                                 // branches for above data
      // TBranch      *bTotalCryEnergy;                                     //
      //
      // if(std::string(argv[1]) == std::string("-c")) // first argument is -c, then the config file name is passed by command line
      // {
      //   for (int i = 3; i < argc ; i++) // run on the remaining arguments to add all the input files
      //   {
      //     std::cout << "Adding file " << argv[i] << std::endl;
      //     tree->Add(argv[i]);
      //   }
      // }
      // else // the config file was indeed the default one
      // {
      //   for (int i = 1; i < argc ; i++) // run on the remaining arguments to add all the input files
      //   {
      //     std::cout << "Adding file " << argv[i] << std::endl;
      //     tree->Add(argv[i]);
      //   }
      // }

      std::vector<int> detector_channels;

      TObjArray *leavescopy = tree->GetListOfLeaves();
      int nLeaves = leavescopy->GetEntries();
      std::vector<std::string> leavesName;
      // fill a vector with the leaves names
      for(int i = 0 ; i < nLeaves ; i++)
      {
        leavesName.push_back(leavescopy->At(i)->GetName());
      }
      // count the entries that start with "ch"
      int numOfCh = 0;
      // int numOfCry = 0;
      std::string ch_prefix("ch");
      std::string t_prefix("t");

      // std::string cry_prefix("cry");
      for(int i = 0 ; i < nLeaves ; i++)
      {
        //     leavesName.push_back(leavescopy->At(i)->GetName());
        if (!leavesName[i].compare(0, ch_prefix.size(), ch_prefix))
        {
          numOfCh++;
          detector_channels.push_back(atoi( (leavesName[i].substr(ch_prefix.size(),leavesName[i].size()-ch_prefix.size())).c_str() )) ;
        }
        // if (!leavesName[i].compare(0, cry_prefix.size(), cry_prefix))
        // numOfCry++;
      }
      //the string "cry" appears 4 times per crystal..
      // numOfCry = numOfCry / 4;
      std::cout << "Detector Channels \t= " << numOfCh << std::endl;
      // std::cout << "Number of Crystals \t= "<< numOfCry << std::endl;


      // first, create the adc channels variables and branches
      // ChainAdcChannel        = new Int_t [numOfCh];
      // ChainDesktopAdcChannel = new Short_t [numOfCh]; // input from ADC desktop
      // ChainVMEadcChannel     = new UShort_t [numOfCh]; // input from VME
      // ChainTimeStamp         = new Float_t[numOfCh];
      // // TDCBinning             = new Float_t[numOfCh];
      // // DigitizerChannelOn     = new bool[adcChannels];
      // bChainAdcChannel       = new TBranch* [numOfCh];
      // bChainTimeStamp        = new TBranch* [numOfCh];
      ULong64_t     ChainExtendedTimeTag;                                // extended time tag
      ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous
      UShort_t      *charge;
      Float_t      *timeStamp;
      TBranch      *bChainExtendedTimeTag;                               // branches for above data
      TBranch      *bChainDeltaTimeTag;                                  // branches for above data
      TBranch      **bCharge;
      TBranch      **btimeStamp;

      charge = new UShort_t[numOfCh];
      timeStamp = new Float_t[numOfCh];
      bCharge = new TBranch*[numOfCh];
      btimeStamp = new TBranch*[numOfCh];

      // set branches for reading the input files
      tree->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag, &bChainExtendedTimeTag);
      tree->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag, &bChainDeltaTimeTag);
      // if(usingRealSimData)
      // {
      //   tree->SetBranchAddress("RealX", &RealX, &bRealX);
      //   tree->SetBranchAddress("RealY", &RealY, &bRealY);
      //   tree->SetBranchAddress("RealZ", &RealZ, &bRealZ);
      //   // fchain->SetBranchAddress("Tagging", &simTaggingCharge, &bsimTaggingCharge);
      //   // fchain->SetBranchAddress("TaggingTimeStamp", &simTaggingTime, &bsimTaggingTime);
      //   tree->SetBranchAddress("CrystalsHit",&CrystalsHit, &bCrystalsHit);
      //   tree->SetBranchAddress("NumbOfInteractions",&NumbOfInteractions, &bNumbOfInteractions);
      //   // fchain->SetBranchAddress("TotalCryEnergy",&TotalCryEnergy, &bTotalCryEnergy);
      // }
      for (int i = 0 ; i < detector_channels.size() ; i++)
      {
        //empty the stringstreams
        std::stringstream sname;
        sname << "ch" << detector_channels[i];
        // stype << "ch" << i << "/F";
        tree->SetBranchAddress(sname.str().c_str(),&charge[detector_channels[i]],&bCharge[detector_channels[i]]);
        sname.str("");
        // stype.str("");

        sname << "t" << detector_channels[i];
        // stype << "t" << i << "/F";
        tree->SetBranchAddress(sname.str().c_str(),&timeStamp[detector_channels[i]],&btimeStamp[detector_channels[i]]);
        sname.str("");
        // stype.str("");
      }



      // find the crystals with complete calibration data
      // this means looking into TWO files:
      // 1. the standard calibration file, but produced with active timing part
      //    this is produced NOT in coincidence and will give all but two plot
      //    for the crystals calibration. the only two missing are the only two that HAVE
      //    to be measured in coincidence, i.e. the tw_correction plots (and it's RMS) and the tag photopeak cut
      // 2. a configuration run performed in coincidence, from which the tw_correction are taken
      //
      // afterwards the crystals will be accepted if they are present in both calibration files

      std::vector<Crystal_t> crystal;
      std::vector<detector_t> detectorSaturation;

      // STANDARD CALIBRATION FILE
      TFile *calibrationFile = new TFile(calibrationFileName.c_str());
      calibrationFile->cd("Module 0.0");
      TList *listModule = gDirectory->GetListOfKeys();
      int nKeysMod = listModule->GetEntries();
      std::vector<std::string> keysModName;
      // fill a vector with the leaves names
      std::string mppc_prefix("MPPC");
      for(int i = 0 ; i < nKeysMod ; i++){
        keysModName.push_back(listModule->At(i)->GetName());
      }

      TCut *taggingPhotopeakCut;
      int taggingCrystalTimingChannel;
      std::string taggingPhotopeakCut_prefix("taggingPhotopeakCut");
      std::string taggingCrystalTimingChannel_prefix("taggingCrystalTimingChannel");
      // std::string det_prefix("channels");
      // std::string saturation_prefix("saturation");

      std::vector<int> *pChannels;
      std::vector<float> *pSaturation;
      std::vector<float> *pPedestal;
      gDirectory->GetObject("channels",pChannels);
      gDirectory->GetObject("saturation",pSaturation);
      gDirectory->GetObject("pedestal",pPedestal);

      std::vector<int> DetChannels = pChannels[0];
      std::vector<float> saturation = pSaturation[0];
      std::vector<float> pedestal = pPedestal[0];

      for(unsigned int iSat = 0; iSat < DetChannels.size(); iSat++)
      {
        detector_t tempDetector;
        tempDetector.digitizerChannel = DetChannels[iSat];
        tempDetector.saturation = saturation[iSat];
        tempDetector.pedestal = pedestal[iSat];
        detectorSaturation.push_back(tempDetector);
      }

      std::vector<std::string> MPPCfolders;
      for(unsigned int i = 0 ; i < keysModName.size() ; i++)
      {
        if (!keysModName[i].compare(0, mppc_prefix.size(), mppc_prefix))
        {
          MPPCfolders.push_back(keysModName[i]);
        }
        if(!keysModName[i].compare(0,taggingPhotopeakCut_prefix.size(),taggingPhotopeakCut_prefix)) // find tcut
        {
          //       std::cout << keysCryName[i] << std::endl;
          taggingPhotopeakCut = (TCut*) gDirectory->Get( keysModName[i].c_str());
          //       if(cut)
          //         temp_crystal.CrystalCut = cut;
        }
        if(!keysModName[i].compare(0,taggingCrystalTimingChannel_prefix.size(),taggingCrystalTimingChannel_prefix)) // find tcut
        {
          std::stringstream snameCh;
          snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
          taggingCrystalTimingChannel = atoi(snameCh.str().c_str());
        }

        // if(!keysModName[i].compare(0,det_prefix.size(),det_prefix)) // find tcut
        // {
        //   std::stringstream snameCh;
        //   std::vector<int> *v;
        //   gDirectory->GetObject("channelsNumRelevantForW",v);
        //   temp_crystal.relevantForW = v[0];
        // }
        //
        // if(!keysModName[i].compare(0,saturation_prefix.size(),saturation_prefix)) // find tcut
        // {
        //   std::stringstream snameCh;
        //   std::vector<int> *v;
        //   gDirectory->GetObject("channelsNumRelevantForW",v);
        //   temp_crystal.relevantForW = v[0];
        // }


      }

      std::stringstream sformulaname;
      sformulaname << "FormulaTag";
      TCut taggingPhotopeakCutName;
      taggingPhotopeakCutName = taggingPhotopeakCut->GetTitle();
      // std::cout << "FormulaTag ------------- \n" << taggingPhotopeakCutName << std::endl;
      TTreeFormula* FormulaTag = new TTreeFormula("FormulaTag",taggingPhotopeakCutName,tree);
      formulas.Add(FormulaTag);



      for(unsigned int iMppc = 0 ; iMppc < MPPCfolders.size() ; iMppc++)
      {
        // std::cout << MPPCfolders[iMppc] << std::endl;
        gDirectory->cd(MPPCfolders[iMppc].c_str());
        TList *listMppc = gDirectory->GetListOfKeys();
        int nKeysMppc = listMppc->GetEntries();
        std::vector<std::string> keysMppcName;
        // fill a vector with the leaves names
        std::string crystal_prefix("Crystal");
        // std::string det_prefix("digitizerChannel");
        // std::string saturation_prefix("saturation");
        for(int i = 0 ; i < nKeysMppc ; i++){
          keysMppcName.push_back(listMppc->At(i)->GetName());
        }

        std::vector<std::string> CrystalFolders;

        for(unsigned int i = 0 ; i < keysMppcName.size() ; i++)
        {
          if (!keysMppcName[i].compare(0, crystal_prefix.size(), crystal_prefix))
          {
            CrystalFolders.push_back(keysMppcName[i]);
          }
          // if (!keysMppcName[i].compare(0, det_prefix.size(), det_prefix))
          // {
          //   std::stringstream snameCh;
          //   snameCh << ((TNamed*) gDirectory->Get(keysMppcName[i].c_str()))->GetTitle();
          //   tempDetector.digitizerChannel = atoi(snameCh.str().c_str());
          // }
          // if (!keysMppcName[i].compare(0, saturation_prefix.size(), saturation_prefix))
          // {
          //   std::stringstream snameCh;
          //   snameCh << ((TNamed*) gDirectory->Get(keysMppcName[i].c_str()))->GetTitle();
          //   tempDetector.saturation = atof(snameCh.str().c_str());
          // }
        }
        // detectorSaturation.push_back(tempDetector);

        for(unsigned int iCry = 0 ; iCry < CrystalFolders.size() ; iCry++)
        {
          //  std::cout << CrystalFolders[iCry] << std::endl;
          gDirectory->cd(CrystalFolders[iCry].c_str());

          Crystal_t temp_crystal;
          temp_crystal.number = -1;
          temp_crystal.CrystalCut = NULL;
          temp_crystal.CrystalCutWithoutCutG = NULL;
          temp_crystal.PhotopeakEnergyCut = NULL;
          temp_crystal.calibrationGraph = NULL;
          temp_crystal.simpleCTR = NULL;
          temp_crystal.centralCTR = NULL;
          temp_crystal.allCTR = NULL;
          temp_crystal.poliCorrCTR = NULL;
          temp_crystal.simpleCTR_norm = NULL;
          temp_crystal.centralCTR_norm = NULL;
          temp_crystal.allCTR_norm = NULL;
          temp_crystal.poliCorrCTR_norm = NULL;
          temp_crystal.wz = NULL;
          temp_crystal.accepted = true;
          temp_crystal.tw_correction = NULL;
          temp_crystal.polishedCorrection = false;



          //get crystal number
          temp_crystal.number = atoi((CrystalFolders[iCry].substr(crystal_prefix.size()+1,CrystalFolders[iCry].size()-crystal_prefix.size()-1)).c_str());

          TList *listCry = gDirectory->GetListOfKeys();
          int nKeysCry = listCry->GetEntries();
          std::vector<std::string> keysCryName;
          if(nKeysCry) //if directory not empty
          {

            for(int i = 0 ; i < nKeysCry ; i++){
              keysCryName.push_back(listCry->At(i)->GetName());
            }
            std::string CalibName;
            std::string CutName;
            std::vector<std::string> cutgNames;
            std::string cutG_prefix("cutg");
            std::string calibration_prefix("Calibration");
            std::string crystalCut_prefix("CrystalCut");
            std::string crystalCutWithoutCutG_prefix("CrystalCutWithoutCutG");
            std::string photopeakEnergyCut_prefix("PhotopeakEnergyCut");
            std::string channel_prefix("digitizerChannel");
            std::string w_channels_prefix("channelsNumRelevantForW");
            std::string timing_channel_prefix("timingChannel");
            std::string wz_prefix("w(z)");
            std::string t_channels_for_poli_prefix("tChannelsForPolishedCorrection");
            std::string mean_for_poli_prefix("meanForPolishedCorrection");
            std::string fwhm_for_poli_prefix("fwhmForPolishedCorrection");
            //  std::string correction_prefix("Correction Graph");
            //  std::string correction_rms_prefix("RMS Correction Graphs");
            for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
            {
              if(!keysCryName[i].compare(0,calibration_prefix.size(),calibration_prefix)) //find calibration graph
              {
                //  std::cout << keysCryName[i] << std::endl;
                TCanvas* C_graph = NULL;
                TGraph *calibGraph = NULL;
                C_graph = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
                if(C_graph)
                calibGraph = (TGraph*) C_graph->GetPrimitive(keysCryName[i].c_str());
                if(calibGraph)
                temp_crystal.calibrationGraph = calibGraph;
              }

              if(!keysCryName[i].compare(0,wz_prefix.size(),wz_prefix)) //find calibration graph
              {
                //  std::cout << keysCryName[i] << std::endl;
                TCanvas* C_graph = NULL;
                TGraph *calibGraph = NULL;
                C_graph = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
                if(C_graph)
                calibGraph = (TGraph*) C_graph->GetPrimitive(keysCryName[i].c_str());
                if(calibGraph)
                temp_crystal.wz = calibGraph;
              }

              if(!keysCryName[i].compare(0,crystalCut_prefix.size(),crystalCut_prefix)) // find tcut
              {
                //  std::cout << keysCryName[i] << std::endl;
                TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
                if(cut)
                temp_crystal.CrystalCut = cut;
              }
              if(!keysCryName[i].compare(0,crystalCutWithoutCutG_prefix.size(),crystalCutWithoutCutG_prefix)) // find tcut
              {
                //  std::cout << keysCryName[i] << std::endl;
                TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
                if(cut)
                temp_crystal.CrystalCutWithoutCutG = cut;
              }
              if(!keysCryName[i].compare(0,photopeakEnergyCut_prefix.size(),photopeakEnergyCut_prefix)) // find tcut
              {
                //  std::cout << keysCryName[i] << std::endl;
                TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
                if(cut)
                temp_crystal.PhotopeakEnergyCut = cut;
              }

              if(!keysCryName[i].compare(0,cutG_prefix.size(),cutG_prefix)) // find tcutgs
              {
                //  std::cout << keysCryName[i] << std::endl;
                TCutG* cutg = (TCutG*) gDirectory->Get(keysCryName[i].c_str());

                temp_crystal.cutg.push_back(cutg);
              }

              if(!keysCryName[i].compare(0,channel_prefix.size(),channel_prefix)) // find detector channel
              {
                //  std::cout << keysCryName[i] << std::endl;
                std::stringstream snameCh;
                snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
                //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
                //  istringstream()
                temp_crystal.detectorChannel = atoi(snameCh.str().c_str());
                //  std::cout <<temp_crystal.detectorChannel << std::endl;
                //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                //  << temp_crystal.detectorChannel << std::endl;
              }






              if(!keysCryName[i].compare(0,timing_channel_prefix.size(),timing_channel_prefix)) // find timing channel
              {
                //  std::cout << keysCryName[i] << std::endl;
                std::stringstream snameCh;
                snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
                //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
                //  istringstream()
                temp_crystal.timingChannel = atoi(snameCh.str().c_str());
                //  std::cout <<temp_crystal.detectorChannel << std::endl;
                //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                //  << temp_crystal.detectorChannel << std::endl;
              }

              if(!keysCryName[i].compare(0,w_channels_prefix.size(),w_channels_prefix)) // find detector channel
              {
                //  std::cout << keysCryName[i] << std::endl;
                std::stringstream snameCh;
                std::vector<int> *v;
                gDirectory->GetObject("channelsNumRelevantForW",v);
                // snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
                //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
                //  istringstream()
                temp_crystal.relevantForW = v[0];
                //  std::cout <<temp_crystal.detectorChannel << std::endl;
                //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                //  << temp_crystal.detectorChannel << std::endl;
              }

              if(!keysCryName[i].compare(0,t_channels_for_poli_prefix.size(),t_channels_for_poli_prefix))
              {
                std::vector<int> *v;
                gDirectory->GetObject("tChannelsForPolishedCorrection",v);
                temp_crystal.polishedCorrection = true;
                temp_crystal.tChannelsForPolishedCorrection = v[0];
              }

              if(!keysCryName[i].compare(0,mean_for_poli_prefix.size(),mean_for_poli_prefix))
              {
                std::vector<double> *v;
                gDirectory->GetObject("meanForPolishedCorrection",v);
                temp_crystal.meanForPolishedCorrection = v[0];
              }

              if(!keysCryName[i].compare(0,fwhm_for_poli_prefix.size(),fwhm_for_poli_prefix))
              {
                std::vector<double> *v;
                gDirectory->GetObject("fwhmForPolishedCorrection",v);
                temp_crystal.fwhmForPolishedCorrection = v[0];
              }



            }

            bool dirExists;
            std::stringstream sname;
            dirExists = gDirectory->cd("TimeCorrection");
            if(dirExists)
            {
              sname.str("");
              sname << "Central correction - Crystal " << temp_crystal.number;
              temp_crystal.centralCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
              sname.str("");
              sname << "Full correction - Crystal " << temp_crystal.number;
              temp_crystal.allCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
              sname.str("");

              temp_crystal.path = gDirectory->GetPath();
              TList *listTcorr = gDirectory->GetListOfKeys();
              int nKeysTcorr = listTcorr->GetEntries();
              std::vector<std::string> keysTcorrName;
              if(nKeysTcorr) //if directory not empty
              {
                for(int i = 0 ; i < nKeysTcorr ; i++){
                  keysTcorrName.push_back(listTcorr->At(i)->GetName());
                }

                std::string deltaWGraph_prefix = "DeltaW Graph";
                std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";
                std::string graph_delay_prefix = "Graph Delay ch_";
                std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
                std::string delay_timing_ch_prefix = "delayTimingChannels";
                for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
                {
                  if(!keysTcorrName[i].compare(0,deltaWGraph_prefix.size(),deltaWGraph_prefix))
                  {
                    TGraphErrors *calibGraph = NULL;
                    calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                    if(calibGraph)
                    {
                      temp_crystal.tw_correction = calibGraph;
                      //fit with straight line
                      TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                      calibGraph->Fit(line,"Q");
                      temp_crystal.tw_correction_line = line;

                    }

                  }

                  if(!keysTcorrName[i].compare(0,rms_deltaWGraph_prefix.size(),rms_deltaWGraph_prefix))
                  {
                    TGraphErrors *calibGraph = NULL;
                    calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                    if(calibGraph)
                    {
                      temp_crystal.rms_tw_correction = calibGraph;
                      //fit with straight line
                      TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                      calibGraph->Fit(line,"Q");
                      temp_crystal.rms_tw_correction_line = line;

                    }

                  }

                  if(!keysTcorrName[i].compare(0,graph_delay_prefix.size(),graph_delay_prefix))
                  {
                    TGraphErrors *calibGraph = NULL;
                    calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                    if(calibGraph)
                    {
                      // -- check if the channel is not excluded by the user
                      // extract next 2 characters
                      //
                      std::string str2 = keysTcorrName[i].substr (graph_delay_prefix.size(),6);     // take a string with next 6 characters after the prefix
                      std::size_t found = str2.find_first_of("ch_");                                  // find next "_"
                      std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                      // std::cout << keysTcorrName[i] << " " << str2 << " " << str3 << std::endl;     // output
                      int current_ch = atoi(str3.c_str());                                          // transform in int
                      // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output

                      bool acceptCh = true;
                      for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                      {
                        if(current_ch == forbidden_channels[iForb])
                        {
                          acceptCh = false;
                        }
                      }

                      if(acceptCh)               // add graph if the ch is accepted
                      {
                        temp_crystal.delay.push_back(calibGraph);
                        //fit with straight line
                        TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                        calibGraph->Fit(line,"Q");
                        temp_crystal.delay_line.push_back(line);
                      }



                    }

                  }

                  if(!keysTcorrName[i].compare(0,rms_graph_delay_prefix.size(),rms_graph_delay_prefix))
                  {
                    TGraphErrors *calibGraph = NULL;
                    calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());


                    if(calibGraph)
                    {
                      // -- check if the channel is not excluded by the user
                      // extract next 2 characters
                      //
                      std::string str2 = keysTcorrName[i].substr (rms_graph_delay_prefix.size(),6);     // take a string with next 6 characters after the prefix
                      std::size_t found = str2.find_first_of("_");                                  // find next "_"
                      std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                      int current_ch = atoi(str3.c_str());                                          // transform in int
                      // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output


                      bool acceptCh = true;
                      for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                      {
                        if(current_ch == forbidden_channels[iForb])
                        {
                          acceptCh = false;
                        }
                      }

                      if(acceptCh)               // add graph if the ch is accepted
                      {
                        temp_crystal.rms_delay.push_back(calibGraph);
                        //fit with straight line
                        TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                        calibGraph->Fit(line,"Q");
                        temp_crystal.rms_delay_line.push_back(line);
                      }
                    }

                  }

                  if(!keysTcorrName[i].compare(0,delay_timing_ch_prefix.size(),delay_timing_ch_prefix))
                  {

                    // // -- check if the channel is not excluded by the user
                    // // extract next 2 characters
                    // //
                    // std::string str2 = keysTcorrName[i].substr (delay_timing_ch_prefix.size(),6);     // take a string with next 6 characters after the prefix
                    // std::size_t found = str2.find_first_of("_");                                  // find next "_"
                    // std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                    // // std::cout << keysTcorrName[i] << " " << str2 << " " << str3 << std::endl;     // output
                    // int current_ch = atoi(str3.c_str());                                          // transform in int
                    // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output
                    //
                    // bool acceptCh = true;
                    // for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                    // {
                    //   if(current_ch == forbidden_channels[iForb])
                    //   {
                    //     acceptCh = false;
                    //   }
                    // }
                    //
                    // if(acceptCh)               // add graph if the ch is accepted
                    // {
                    //  std::cout << keysCryName[i] << std::endl;
                    std::stringstream snameCh;
                    std::vector<int> *v;
                    gDirectory->GetObject("delayTimingChannels",v);
                    // snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
                    //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
                    //  istringstream()
                    temp_crystal.delayTimingChannels = v[0];
                    //  std::cout <<temp_crystal.detectorChannel << std::endl;
                    //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                    //  << temp_crystal.detectorChannel << std::endl;
                    // }

                  }
                }
              }
              gDirectory->cd("..");
            }

            bool likeExists;
            // std::stringstream sname;
            likeExists = gDirectory->cd("LikelihoodCorrection");
            if(likeExists)
            {

              sname << "T vs w after likelihood correction - crystal " << temp_crystal.number;
              temp_crystal.scatter_likeT_w = new TH2F(sname.str().c_str(),sname.str().c_str(),250,0,1,histoBins,histoMin,histoMax);
              sname.str("");

              TList *listTcorr = gDirectory->GetListOfKeys();
              int nKeysTcorr = listTcorr->GetEntries();
              std::vector<std::string> keysTcorrName;
              if(nKeysTcorr) //if directory not empty
              {
                for(int i = 0 ; i < nKeysTcorr ; i++){
                  keysTcorrName.push_back(listTcorr->At(i)->GetName());
                }

                std::string Mu_W_prefix = "histo_d_W_ch_";

                for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
                {

                  if(!keysTcorrName[i].compare(0,Mu_W_prefix.size(),Mu_W_prefix))
                  {
                    TH2F *calibGraph = NULL;
                    calibGraph = (TH2F*) gDirectory->Get(keysTcorrName[i].c_str());
                    if(calibGraph)
                    {
                      // -- check if the channel is not excluded by the user
                      // extract next 2 characters
                      //
                      std::string str2 = keysTcorrName[i].substr (Mu_W_prefix.size(), keysTcorrName[i].size()-Mu_W_prefix.size());     // take a string with next 6 characters after the prefix
                      std::size_t found = str2.find_first_of("_");                                  // find next "_"
                      std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                      // std::cout << keysTcorrName[i] << " " << str2 << " " << str3 << std::endl;     // output
                      int current_ch = atoi(str3.c_str());                                          // transform in int
                      // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output

                      bool acceptCh = true;
                      for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                      {
                        if(current_ch == forbidden_channels[iForb])
                        {
                          acceptCh = false;
                        }
                      }

                      if(acceptCh)               // add graph if the ch is accepted
                      {
                        // find timing channels number
                        std::size_t found2 = str2.find_last_of("_");

                        std::string str4 = str2.substr(found2+1,str2.size()-found2);
                        int current_t =  atoi(str4.c_str());
                        // std::cout << str2 << "\t" << current_t << std::endl;
                        LikelihoodElement_t temp_like;
                        temp_like.t_channel = current_t;
                        temp_like.mu_w = calibGraph;
                        temp_crystal.likelihood_mu_w.push_back(temp_like);
                        //fit with straight line
                        // TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                        // calibGraph->Fit(line,"Q");
                        // temp_crystal.delay_line.push_back(line);
                      }

                    }
                  }
                }
              }
              gDirectory->cd("..");
            }



            sname << "No correction - Crystal " << temp_crystal.number;
            temp_crystal.simpleCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
            sname.str("");
            sname << "Polished correction - Crystal " << temp_crystal.number;
            temp_crystal.poliCorrCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
            sname.str("");
            TCut globalCut; // the cut for the formula
            globalCut += temp_crystal.CrystalCutWithoutCutG->GetTitle();     // this is BasicCut (XYZ and taggingPhotopeak) + CutTrigger (TriggerChannel and broadcut)
            globalCut += temp_crystal.PhotopeakEnergyCut->GetTitle();        // this is the cut on photopeak energy of the corrected spectrum for this crystal
            for(unsigned int iCutg = 0; iCutg < temp_crystal.cutg.size(); iCutg++)
            {
              globalCut += temp_crystal.cutg[iCutg]->GetName();              // these are the two cutg for this crystal
            }
            sname.str("");

            sname << "Formula" << temp_crystal.number;
            TTreeFormula* Formula = new TTreeFormula(sname.str().c_str(),globalCut,tree);
            formulas.Add(Formula);
            temp_crystal.Formula = Formula;
            sname.str("");

            // if(temp_crystal.calibrationGraph && temp_crystal.CrystalCutWithoutCutG && temp_crystal.PhotopeakEnergyCut && (temp_crystal.cutg.size() == 2))
            if(temp_crystal.calibrationGraph && temp_crystal.CrystalCutWithoutCutG && temp_crystal.PhotopeakEnergyCut)
            {
              crystal.push_back(temp_crystal);
            }

          }
          gDirectory->cd("..");
        }
        calibrationFile->cd("Module 0.0");
      }





      std::cout << "Calibration data found for crystals: " << std::endl;
      for(unsigned int i = 0 ;  i < crystal.size() ; i++)
      {
        if(crystal[i].accepted)
        {
          std::cout << crystal[i].number << std::endl;

          // and prepare correlation matrix scatter plots
          crystal[i].correlation_matrix_element = new TH2F**[crystal[i].likelihood_mu_w.size()];
          for (unsigned int iLike = 0 ; iLike < crystal[i].likelihood_mu_w.size(); iLike++)
          {
            crystal[i].correlation_matrix_element[iLike] = new TH2F*[crystal[i].likelihood_mu_w.size()];
            for (unsigned int jLike = 0 ; jLike < crystal[i].likelihood_mu_w.size(); jLike++)
            {
              std::stringstream tname;
              tname << "Cry_"<< crystal[i].number
                    << "_histo_c_" << iLike+1 << "_" << jLike+1
                    << "_t_" << crystal[i].likelihood_mu_w[iLike].t_channel
                    << "_t_" << crystal[i].likelihood_mu_w[jLike].t_channel;
              crystal[i].correlation_matrix_element[iLike][jLike] = new TH2F(tname.str().c_str(),tname.str().c_str(),250,0,1,corrBins,corrMin,corrMax);
            }
          }














          std::stringstream tname;


          tname << "Cry_"<< crystal[i].number << "_likelihood_ctr";




              // crystal[i].correlation_matrix_element[iLike][jLike] = new
          //and also final likelihood ctr histos
          crystal[i].likelihood_CTR = new TH1F(tname.str().c_str(),tname.str().c_str(),histoBins,histoMin,histoMax);

          // now prepare the w limits, based on slicing of z coordinate
          double doiLimitDownFraction = 0.01;
          double doiLimitUpFraction = 0.99;
          // int sections = 10;

          double crystaLenght = crystal[i].calibrationGraph->Eval(0);
          double min_z = doiLimitDownFraction*crystaLenght;
          double max_z = doiLimitUpFraction*crystaLenght;







          // double step = (max_z -min_z)/sections;
          // //remember that w and y go in opposite directions
          // // crystal[i].mu_histo = new TH1F**[sections];
          // for(int iSteps = sections ; iSteps > 0 ; iSteps--)
          // {
          //   //slicing in equal z spaces
          //   double yPos_low = min_z + (iSteps-1)*step;
          //   double yPos_up  = min_z + (iSteps)*step;
          //   double w_max  = crystal[i].wz->Eval(yPos_low);
          //   double w_min  = crystal[i].wz->Eval(yPos_up);
          //   crystal[i].w_limit_low.push_back(w_min);
          //   crystal[i].w_limit_up.push_back(w_max);
          //
          //   // slicing in equal w spaces
          //   // crystal[i].w_limit_low.push_back(equal_w_min + );
          //   // crystal[i].w_limit_up.push_back(w_max);
          //
          //
          //
          //
          //
          //   // // debug
          //   // std::cout << yPos_up
          //   //           << "\t"
          //   //           << w_min
          //   //           << yPos_low
          //   //           << "\t"
          //   //           << w_max
          //   //           << "\t"
          //   //           << std::endl;
          // }


          // slicing in equal w spaces
          double equal_w_min = crystal[i].wz->Eval(max_z);;
          double equal_w_max = crystal[i].wz->Eval(min_z);;
          double step_equal_w = (equal_w_max -equal_w_min)/sections;
          std::cout << step_equal_w << " " << equal_w_min << " " << equal_w_max << std::endl;

          for(int iSteps = 0 ; iSteps < sections ; iSteps++)
          {
            crystal[i].w_limit_low.push_back(equal_w_min + (step_equal_w*iSteps));
            crystal[i].w_limit_up.push_back(equal_w_min + (step_equal_w*(iSteps+1)));
          }



          //debug
          for(unsigned int iW = 0 ;iW < crystal[i].w_limit_low.size() ; iW++)
          {
            std::cout << crystal[i].w_limit_low[iW] << "\t";
          }
          std::cout << std::endl;
          for(unsigned int iW = 0 ;iW < crystal[i].w_limit_up.size() ; iW++)
          {
            std::cout << crystal[i].w_limit_up[iW] << "\t";
          }
          std::cout << std::endl;


          // h_correlation_matrix_element_w
          crystal[i].h_correlation_matrix_element_w = new TH1F***[crystal[i].w_limit_low.size()];
          for(unsigned int iW = 0 ;iW < crystal[i].w_limit_low.size() ; iW++)
          {
            crystal[i].h_correlation_matrix_element_w[iW] = new TH1F**[crystal[i].likelihood_mu_w.size()];
            for (unsigned int iLike = 0 ; iLike < crystal[i].likelihood_mu_w.size(); iLike++)
            {
              crystal[i].h_correlation_matrix_element_w[iW][iLike] = new TH1F*[crystal[i].likelihood_mu_w.size()];
              for (unsigned int jLike = 0 ; jLike < crystal[i].likelihood_mu_w.size(); jLike++)
              {
                std::stringstream tname;
                tname << "Cry_"<< crystal[i].number
                      << "_histo_c_" << iLike+1 << "_" << jLike+1
                      << "_w_" << iW
                      << "_t_" << crystal[i].likelihood_mu_w[iLike].t_channel
                      << "_t_" << crystal[i].likelihood_mu_w[jLike].t_channel;
                crystal[i].h_correlation_matrix_element_w[iW][iLike][jLike] = new TH1F(tname.str().c_str(),tname.str().c_str(),corrBins,corrMin,corrMax);
                // std::cout << tname.str() << std::endl;
              }
            }
          }


          // and prepare an matrix mu_i_w
          crystal[i].mu_histo =  new TH1F**[crystal[i].likelihood_mu_w.size()];
          crystal[i].mu_i_w = new double*[crystal[i].likelihood_mu_w.size()];
          crystal[i].mu_i_w_counter = new long int*[crystal[i].likelihood_mu_w.size()];
          for (unsigned int iLike = 0 ; iLike < crystal[i].likelihood_mu_w.size(); iLike++)
          {
            crystal[i].mu_histo[iLike] =  new TH1F*[sections];
            crystal[i].mu_i_w[iLike] = new double[sections];
            crystal[i].mu_i_w_counter[iLike] = new long int[sections];
            for(int iW = 0 ; iW < sections ; iW++)
            {
              std::stringstream sname;

              sname << "Mu histo - det " << iLike << " w " << iW;
              crystal[i].mu_i_w[iLike][iW] = 0;
              crystal[i].mu_i_w_counter[iLike][iW] = 0;
              crystal[i].mu_histo[iLike][iW] = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
              sname.str("");
            }
          }

          //and prepare the correlation matrix array for this crystal
          crystal[i].corr_mat = new gsl_matrix*[sections];
          crystal[i].inverse = new gsl_matrix*[sections];
          crystal[i].norm_corr_mat = new gsl_matrix*[sections];
          crystal[i].corr_mat_counter = new long int**[sections];
          crystal[i].normalization_factor = new double[sections];


          for (int iW = 0 ; iW < sections ; iW++)
          {
            crystal[i].corr_mat[iW] = gsl_matrix_alloc(crystal[i].likelihood_mu_w.size(), crystal[i].likelihood_mu_w.size());
            crystal[i].norm_corr_mat[iW] = gsl_matrix_alloc(crystal[i].likelihood_mu_w.size(), crystal[i].likelihood_mu_w.size());
            crystal[i].normalization_factor[iW] = 0;
            // crystal[i].inverse[iW] = gsl_matrix_alloc(crystal[i].likelihood_mu_w.size(), crystal[i].likelihood_mu_w.size());
            // and set to 0
            crystal[i].corr_mat_counter[iW] = new long int*[crystal[i].likelihood_mu_w.size()];
            for (unsigned int iLike = 0 ; iLike < crystal[i].likelihood_mu_w.size(); iLike++)
            {
              crystal[i].corr_mat_counter[iW][iLike] = new long int[crystal[i].likelihood_mu_w.size()];
              for (unsigned int jLike = 0 ; jLike < crystal[i].likelihood_mu_w.size(); jLike++)
              {
                crystal[i].corr_mat_counter[iW][iLike][jLike] = 0;
                size_t iValue = (size_t) iLike;
                size_t jValue = (size_t) jLike;
                gsl_matrix_set(crystal[i].corr_mat[iW], iValue, jValue, 0);
                gsl_matrix_set(crystal[i].norm_corr_mat[iW], iValue, jValue, 0);
              }
            }


          }
          // end of debug
        }
      }

      //notify TTreeFormula(s) to TChain
      tree->SetNotify(&formulas);

      std::cout << "taggingCrystalTimingChannel = " << taggingCrystalTimingChannel << std::endl;

      // LOOP on all calibration events
      // to calculate the mu_i for each detector
      long long int nevent = tree->GetEntries();
      std::cout << "Total number of events = " << nevent << std::endl;
      long int goodEvents = 0;
      long int counter = 0;
      // for (long long int i=0;i<1000000;i++)
      for (long long int i=0;i<nevent;i++)
      {
        // std::cout << "Event " << i << std::endl;
        tree->GetEvent(i);              //read complete accepted event in memory
        for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
        {
          if(crystal[iCry].accepted)
          {
            if(FormulaTag->EvalInstance())
            {
              if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
              {

                goodEvents++;

                //calculate FloodZ...
                Float_t FloodZ;
                float centralChargeOriginal;
                float centralSaturation;
                float centralPedestal;
                Float_t division = 0.0;

                centralChargeOriginal = charge[crystal[iCry].detectorChannel];
                for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                {
                  if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
                  {
                    centralSaturation = detectorSaturation[iSat].saturation;
                    centralPedestal = detectorSaturation[iSat].pedestal;
                  }
                }
                float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );

                for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
                {
                  float originalCh = charge[crystal[iCry].relevantForW[iW]];
                  float saturationCh;
                  float pedestalCorr;
                  for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                  {
                    if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
                    {
                      saturationCh = detectorSaturation[iSat].saturation;
                      pedestalCorr = detectorSaturation[iSat].pedestal;
                    }
                  }
                  division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
                }

                FloodZ = centralChargeCorr / division;
                // std::cout << FloodZ << std::endl;

                for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
                {
                  if( (FloodZ > crystal[iCry].w_limit_low[iW]) && (FloodZ <= crystal[iCry].w_limit_up[iW])) // look for w slice for this w
                  {
                    for(unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
                    {
                      int iTiming =  crystal[iCry].likelihood_mu_w[iLike].t_channel;
                      if((timeStamp[iTiming] != 0) &&
                         (timeStamp[taggingCrystalTimingChannel] != 0) &&
                         ((timeStamp[iTiming] - timeStamp[taggingCrystalTimingChannel]) > histoMin) &&
                         ((timeStamp[iTiming] - timeStamp[taggingCrystalTimingChannel]) < histoMax)
                        )
                      {
                        crystal[iCry].mu_i_w_counter[iLike][iW]++; // increase count of w in this detector and w slice
                        crystal[iCry].mu_i_w[iLike][iW] += timeStamp[iTiming] - timeStamp[taggingCrystalTimingChannel];      // increase the partial sum
                        crystal[iCry].mu_histo[iLike][iW]->Fill(timeStamp[iTiming] - timeStamp[taggingCrystalTimingChannel]);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        counter++;

        int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
        if( (perc % 10) == 0 )
        {
          std::cout << "\r";
          std::cout << perc << "% done... ";
          //std::cout << counter << std::endl;
        }
      }
      std::cout << "Good events = " << goodEvents << std::endl;



      for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
      {
        if(crystal[iCry].accepted)
        {
          std::cout << "Crystal " << crystal[iCry].number << std::endl;

          // // debug - output the counters
          // for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
          // {
          //   // std::cout << "------------------------------------" << std::endl;
          //   // std::cout << "w slice " << iW << std::endl;
          //   for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
          //   {
          //     std::cout << crystal[iCry].mu_i_w_counter[iLike][iW] << " ";
          //     // crystal[iCry].mu_i_w[iLike][iW] = crystal[iCry].mu_i_w[iLike][iW] / crystal[iCry].mu_i_w_counter[iLike][iW];
          //   }
          //   std::cout << std::endl;
          // }

          // divide partial sum by counters
          for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
          {
            for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
            {
              crystal[iCry].mu_i_w[iLike][iW] = crystal[iCry].mu_i_w[iLike][iW] / crystal[iCry].mu_i_w_counter[iLike][iW];
            }
          }

          // debug - output final matrix
          for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
          {
            // std::cout << "------------------------------------" << std::endl;
            // std::cout << "w slice " << iW << std::endl;
            std::cout << std::setw(8) << (crystal[iCry].w_limit_up[iW] + crystal[iCry].w_limit_low[iW]) / 2.0 << " ";
            for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
            {
              std::cout << std::setw(8) <<crystal[iCry].mu_i_w[iLike][iW] << " ";
              // crystal[iCry].mu_i_w[iLike][iW] = crystal[iCry].mu_i_w[iLike][iW] / crystal[iCry].mu_i_w_counter[iLike][iW];
            }
            std::cout << std::endl;
          }
        }
      }


      //debug
      // output the mu_i_w




      // AGAIN LOOP
      // now loop and calculate the correlation matrices, one per w, each matrix is 9x9 elements
      counter = 0;
      goodEvents = 0;
      // for (long long int i=0;i<1000000;i++)
      for (long long int i=0;i<nevent;i++)
      {
        // std::cout << "Event " << i << std::endl;
        tree->GetEvent(i);              //read complete accepted event in memory
        for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
        {
          if(crystal[iCry].accepted)
          {
            if(FormulaTag->EvalInstance() ) // if in photopeak of tagging crystal
            {
              if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
              {
                goodEvents++;


                //calculate FloodZ...
                Float_t FloodZ;
                float centralChargeOriginal;
                float centralSaturation;
                float centralPedestal;
                Float_t division = 0.0;

                centralChargeOriginal = charge[crystal[iCry].detectorChannel];
                for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                {
                  if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
                  {
                    centralSaturation = detectorSaturation[iSat].saturation;
                    centralPedestal = detectorSaturation[iSat].pedestal;
                  }
                }
                float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );

                for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
                {
                  float originalCh = charge[crystal[iCry].relevantForW[iW]];
                  float saturationCh;
                  float pedestalCorr;
                  for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                  {
                    if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
                    {
                      saturationCh = detectorSaturation[iSat].saturation;
                      pedestalCorr = detectorSaturation[iSat].pedestal;
                    }
                  }
                  division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
                }

                FloodZ = centralChargeCorr / division;
                // std::cout << FloodZ << std::endl;


                for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
                {
                  if( (FloodZ > crystal[iCry].w_limit_low[iW]) && (FloodZ <= crystal[iCry].w_limit_up[iW])) // look for w slice for this w
                  {
                    for(unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
                    {
                      int iTiming =  crystal[iCry].likelihood_mu_w[iLike].t_channel; // extract t det num
                      double mu_i = crystal[iCry].mu_i_w[iLike][iW];

                      for(unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++) // run on detectors
                      {
                        int jTiming =  crystal[iCry].likelihood_mu_w[jLike].t_channel; // extract t det num
                        double mu_j = crystal[iCry].mu_i_w[jLike][iW];
                        size_t iValue = (size_t) iLike;
                        size_t jValue = (size_t) jLike;

                        // if...
                        if((timeStamp[iTiming] != 0) &&
                           (timeStamp[taggingCrystalTimingChannel] != 0) &&
                           ((timeStamp[iTiming] - timeStamp[taggingCrystalTimingChannel]) > histoMin) &&
                           ((timeStamp[iTiming] - timeStamp[taggingCrystalTimingChannel]) < histoMax) &&
                           (timeStamp[jTiming] != 0) &&
                           // (timeStamp[taggingCrystalTimingChannel] != 0) &&
                           ((timeStamp[jTiming] - timeStamp[taggingCrystalTimingChannel]) > histoMin) &&
                           ((timeStamp[jTiming] - timeStamp[taggingCrystalTimingChannel]) < histoMax)
                          )
                        {
                          // increase counter
                          crystal[iCry].corr_mat_counter[iW][iLike][jLike]++;
                          // extract current value
                          double prev = gsl_matrix_get(crystal[iCry].corr_mat[iW], iValue, jValue);
                          // calc increment
                          double c_ij_partial = ((timeStamp[iTiming] - timeStamp[taggingCrystalTimingChannel]) - mu_i)*((timeStamp[jTiming] - timeStamp[taggingCrystalTimingChannel]) - mu_j);
                          //update value
                          gsl_matrix_set(crystal[iCry].corr_mat[iW], iValue, jValue, prev + c_ij_partial);
                          crystal[iCry].h_correlation_matrix_element_w[iW][iLike][jLike]->Fill(c_ij_partial);
                        }
                      }
                    }
                  }
                }




                // //build correlation matrix for this event...
                // gsl_matrix *corr_mat = gsl_matrix_alloc(crystal[iCry].likelihood_mu_w.size(),
                //                                    crystal[iCry].likelihood_mu_w.size());
                // // Double_t **c_ij_mat = new Double_t*[crystal[iCry].likelihood_mu_w.size()];
                // // Float_t** correlation_matrix = new Float_t*[crystal[iCry].likelihood_mu_w.size()];
                // for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++)
                // {
                //
                //   int iTiming =  crystal[iCry].likelihood_mu_w[iLike].t_channel;
                //
                //
                //   for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++)
                //   {
                //     // extract j timing channel
                //     int jTiming =  crystal[iCry].likelihood_mu_w[jLike].t_channel;
                //
                //     // calc c_ij
                //     Int_t w_bin = crystal[iCry].correlation_matrix_element[iLike][jLike]->GetXaxis()->FindBin(FloodZ);
                //
                //     Int_t w_bin_min = w_bin - 0;
                //     Int_t w_bin_max = w_bin + 0;
                //     TH1D* projh2Y_slice = crystal[iCry].correlation_matrix_element[iLike][jLike]->ProjectionY("slice",w_bin_min,w_bin_max);
                //     Double_t c_ij = (Double_t) projh2Y_slice->GetMean();
                //
                //     // std::cout << iLike << "\t" << jLike << "\t" << c_ij << std::endl;
                //     size_t iValue = (size_t) iLike;
                //     size_t jValue = (size_t) jLike;
                //
                //     // c_ij_mat[iLike][jLike] = c_ij;
                //
                //
                //     // correlation_matrix[iLike][jLike] = c_ij;
                //     gsl_matrix_set(corr_mat, iValue, jValue, (double) c_ij);
                //
                //
                //
                //
                //   }
                // }
                //
                // // print correlation matrix
                //
                // // std::cout << "-------------------"<< std::endl;
                // // for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++)
                // // {
                // //   std::cout << crystal[iCry].likelihood_mu_w[iLike].t_channel << "\t";
                // // }
                // // std::cout <<  std::endl;
                // // print_mat_contents(corr_mat,crystal[iCry].likelihood_mu_w.size());
                // // std::cout << std::endl;
                //
                // //check if matrix is singular, in which case skip the rest
                // // double det = determinantOfMatrix(c_ij_mat, crystal[iCry].likelihood_mu_w.size());
                // // std::cout << det << std::endl;
                // // double det;
                // int signum;
                // gsl_permutation *p = gsl_permutation_alloc(crystal[iCry].likelihood_mu_w.size());
                // gsl_matrix *tmpA = gsl_matrix_alloc(crystal[iCry].likelihood_mu_w.size(),crystal[iCry].likelihood_mu_w.size());
                // gsl_matrix_memcpy(tmpA , corr_mat);
                //
                // gsl_linalg_LU_decomp(tmpA , p , &signum);
                // double det = gsl_linalg_LU_det(tmpA , signum);
                // gsl_permutation_free(p);
                // // std::cout << det << std::endl;
                //
                //
                // // invert correlation matrix
                // if(det != 0)
                // {
                //   gsl_matrix *inverse = invert_a_matrix(corr_mat,crystal[iCry].likelihood_mu_w.size());
                //
                //   // print inverse matrix
                //   // std::cout << std::endl;
                //   // print_mat_contents(inverse,crystal[iCry].likelihood_mu_w.size());
                //   // std::cout << std::endl;
                //
                //   //calc normalization factor
                //   // just sum of all elements in the inverse matrix...
                //   double normalization_factor = 0.0;
                //   for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++)
                //   {
                //     for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++)
                //     {
                //       size_t iValue = (size_t) iLike;
                //       size_t jValue = (size_t) jLike;
                //       normalization_factor += gsl_matrix_get(inverse, iValue, jValue);
                //     }
                //   }
                //
                //   // calc best time estimator
                //   double best_estimator = 0.0;
                //   for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++)
                //   {
                //     // extract i timing channel
                //     int jTiming =  crystal[iCry].likelihood_mu_w[jLike].t_channel;
                //     Int_t w_bin_j = crystal[iCry].likelihood_mu_w[jLike].mu_w->GetXaxis()->FindBin(FloodZ);
                //
                //     //calc mu_i
                //     // set some limits for acceptance of w. let's start with +/- 2 bins
                //     Int_t w_bin_j_min = w_bin_j - 0;
                //     Int_t w_bin_j_max = w_bin_j + 0;
                //
                //     size_t jValue = (size_t) jLike;
                //
                //     TH1D* projh2Y_j_slice = crystal[iCry].likelihood_mu_w[jLike].mu_w->ProjectionY("slicej",w_bin_j_min,w_bin_j_max);
                //     Float_t mu_j = projh2Y_j_slice->GetMean();
                //
                //     double weight = 0.0;
                //     for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++)
                //     {
                //       size_t iValue = (size_t) iLike;
                //       weight += gsl_matrix_get(inverse, iValue, jValue);
                //     }
                //     best_estimator += weight * ((timeStamp[jTiming] - timeStamp[taggingCrystalTimingChannel]) - mu_j);
                //   }
                //   best_estimator = best_estimator / normalization_factor;
                //   // std::cout << best_estimator << std::endl;
                //   crystal[iCry].likelihood_CTR->Fill(best_estimator);
                // }

              }
            }
          }
        }
        counter++;

        int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
        if( (perc % 10) == 0 )
        {
          std::cout << "\r";
          std::cout << perc << "% done... ";
          //std::cout << counter << std::endl;
        }
      }
      std::cout << "Good events = " << goodEvents << std::endl;


      // calc final corr matrix, dividing by counter
      for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
      {
        if(crystal[iCry].accepted)
        {
          std::cout << "Crystal " << crystal[iCry].number << std::endl;

          // // debug - output the counters
          // for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
          // {
          //   // std::cout << "------------------------------------" << std::endl;
          //   // std::cout << "w slice " << iW << std::endl;
          //   for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
          //   {
          //     std::cout << crystal[iCry].mu_i_w_counter[iLike][iW] << " ";
          //     // crystal[iCry].mu_i_w[iLike][iW] = crystal[iCry].mu_i_w[iLike][iW] / crystal[iCry].mu_i_w_counter[iLike][iW];
          //   }
          //   std::cout << std::endl;
          // }

          // divide partial sum by counters
          for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
          {
            for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
            {
              for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++) // run on detectors
              {
                size_t iValue = (size_t) iLike;
                size_t jValue = (size_t) jLike;
                // extract current value
                double prev = gsl_matrix_get(crystal[iCry].corr_mat[iW], iValue, jValue);
                // calc increment
                double final_v = prev / crystal[iCry].corr_mat_counter[iW][iValue][jValue];
                //update value
                // method 1, using the average
                // gsl_matrix_set(crystal[iCry].corr_mat[iW], iValue, jValue, final_v);
                // method 2, taking mean of Histograms
                double mean_from_histo = crystal[iCry].h_correlation_matrix_element_w[iW][iLike][jLike]->GetMean();
                gsl_matrix_set(crystal[iCry].corr_mat[iW], iValue, jValue, mean_from_histo);
                // crystal[iCry].h_correlation_matrix_element_w[iW][iLike][jLike]->Fill(c_ij_partial);

              }
            }
          }

          // debug - output final matrix

          std::cout << "|----------------------------------|" << std::endl;
          std::cout << "|          Covariance matrix       |" << std::endl;
          std::cout << "|----------------------------------|" << std::endl;
          for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
          {
            std::cout << "------------------------------------" << std::endl;
            std::cout << "w slice " << iW << std::endl;
            for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
            {
              for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++) // run on detectors
              {
                size_t iValue = (size_t) iLike;
                size_t jValue = (size_t) jLike;
                std::cout << gsl_matrix_get(crystal[iCry].corr_mat[iW], iValue, jValue) << "\t";
              }
              std::cout << std::endl;
              // std::cout << crystal[iCry].mu_i_w[iLike][iW] << " ";
              // crystal[iCry].mu_i_w[iLike][iW] = crystal[iCry].mu_i_w[iLike][iW] / crystal[iCry].mu_i_w_counter[iLike][iW];
            }
            std::cout << std::endl;
          }


          std::cout << "|----------------------------------|" << std::endl;
          std::cout << "|          norm covariance         |" << std::endl;
          std::cout << "|----------------------------------|" << std::endl;
          for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
          {
            // std::cout << "------------------------------------" << std::endl;
            // std::cout << "w slice " << iW << std::endl;
            for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
            {
              for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++) // run on detectors
              {
                size_t iValue = (size_t) iLike;
                size_t jValue = (size_t) jLike;
                double cov = gsl_matrix_get(crystal[iCry].corr_mat[iW], iValue, jValue);
                double sigma_i =  gsl_matrix_get(crystal[iCry].corr_mat[iW], iValue, iValue);
                double sigma_j =  gsl_matrix_get(crystal[iCry].corr_mat[iW], jValue, jValue);
                double div_factor = sqrt(sigma_i*sigma_j);
                gsl_matrix_set(crystal[iCry].norm_corr_mat[iW], iValue, jValue, cov/div_factor);
                // std::cout << gsl_matrix_get(crystal[iCry].corr_mat[iW], iValue, jValue) << "\t";
              }
              // std::cout << std::endl;
              // std::cout << crystal[iCry].mu_i_w[iLike][iW] << " ";
              // crystal[iCry].mu_i_w[iLike][iW] = crystal[iCry].mu_i_w[iLike][iW] / crystal[iCry].mu_i_w_counter[iLike][iW];
            }
            // std::cout << std::endl;
          }


          for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
          {
            //now in principle the determinants should never be 0...
            // crystal[iCry].inverse[iW] = invert_a_matrix(crystal[iCry].corr_mat[iW],crystal[iCry].likelihood_mu_w.size());
            // and print it

            std::cout << "------------------------------------" << std::endl;
            std::cout << "w slice " << iW << std::endl;
            for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
            {
              for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++) // run on detectors
              {
                size_t iValue = (size_t) iLike;
                size_t jValue = (size_t) jLike;
                std::cout << std::setw(8) << gsl_matrix_get(crystal[iCry].norm_corr_mat[iW], iValue, jValue) << "\t";
              }
              std::cout << std::endl;
              // std::cout << crystal[iCry].mu_i_w[iLike][iW] << " ";
              // crystal[iCry].mu_i_w[iLike][iW] = crystal[iCry].mu_i_w[iLike][iW] / crystal[iCry].mu_i_w_counter[iLike][iW];
            }
            std::cout << std::endl;
          }





          // invert correlation matrices
          for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
          {
            //now in principle the determinants should never be 0...
            crystal[iCry].inverse[iW] = invert_a_matrix(crystal[iCry].corr_mat[iW],crystal[iCry].likelihood_mu_w.size());
            // and calc the sum of all elements (normalization_factor)

            // std::cout << "------------------------------------" << std::endl;
            // std::cout << "w slice " << iW << std::endl;
            for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
            {
              for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++) // run on detectors
              {
                size_t iValue = (size_t) iLike;
                size_t jValue = (size_t) jLike;
                crystal[iCry].normalization_factor[iW] += gsl_matrix_get(crystal[iCry].inverse[iW], iValue, jValue);
            //     std::cout << gsl_matrix_get(crystal[iCry].corr_mat[iW], iValue, jValue) << "\t";
              }
            //   std::cout << std::endl;
            //   // std::cout << crystal[iCry].mu_i_w[iLike][iW] << " ";
            //   // crystal[iCry].mu_i_w[iLike][iW] = crystal[iCry].mu_i_w[iLike][iW] / crystal[iCry].mu_i_w_counter[iLike][iW];
            }
            // std::cout << std::endl;

          }


          // debug - output final matrix
          std::cout << "|----------------------------------|" << std::endl;
          std::cout << "|          Inverse matrix          |" << std::endl;
          std::cout << "|----------------------------------|" << std::endl;
          for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
          {
            //now in principle the determinants should never be 0...
            // crystal[iCry].inverse[iW] = invert_a_matrix(crystal[iCry].corr_mat[iW],crystal[iCry].likelihood_mu_w.size());
            // and print it

            std::cout << "------------------------------------" << std::endl;
            std::cout << "w slice " << iW << std::endl;
            for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++) // run on detectors
            {
              for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++) // run on detectors
              {
                size_t iValue = (size_t) iLike;
                size_t jValue = (size_t) jLike;
                std::cout << gsl_matrix_get(crystal[iCry].inverse[iW], iValue, jValue) << "\t";
              }
              std::cout << std::endl;
              // std::cout << crystal[iCry].mu_i_w[iLike][iW] << " ";
              // crystal[iCry].mu_i_w[iLike][iW] = crystal[iCry].mu_i_w[iLike][iW] / crystal[iCry].mu_i_w_counter[iLike][iW];
            }
            std::cout << std::endl;
          }
        }
      }


      // final LOOP, calculate the bet time estimators
      counter = 0;
      goodEvents = 0;
      // for (long long int i=0;i<1000000;i++)
      for (long long int i=0;i<nevent;i++)
      {
        // std::cout << "Event " << i << std::endl;
        tree->GetEvent(i);              //read complete accepted event in memory
        for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
        {
          if(crystal[iCry].accepted)
          {
            if(FormulaTag->EvalInstance() ) // if in photopeak of tagging crystal
            {
              if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
              {
                // if( (timeStamp[jTiming] != 0) &&
                //     (timeStamp[taggingCrystalTimingChannel] != 0) &&
                //     ( (timeStamp[jTiming] - timeStamp[taggingCrystalTimingChannel]) > histoMin ) &&
                //     ( (timeStamp[jTiming] - timeStamp[taggingCrystalTimingChannel] ) < histoMax )
                //   )
                //   {
                //
                //   }
                bool acceptEvent = true;

                for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++)
                {
                  int iTiming =  crystal[iCry].likelihood_mu_w[iLike].t_channel;
                  if((timeStamp[iTiming] == 0) ||
                     ((timeStamp[iTiming] - timeStamp[taggingCrystalTimingChannel]) < histoMin)  ||
                     ((timeStamp[iTiming] - timeStamp[taggingCrystalTimingChannel]) > histoMax ))
                  {
                    acceptEvent = false;
                  }
                }
                if(timeStamp[taggingCrystalTimingChannel] == 0)
                {
                  acceptEvent = false;
                }
                if(acceptEvent)
                {
                  goodEvents++;

                  double simpleCTR = timeStamp[crystal[iCry].timingChannel] - timeStamp[taggingCrystalTimingChannel];
                  crystal[iCry].simpleCTR->Fill(simpleCTR);

                  double avgCTR = 0.0;
                  double sumWeight = 0.0;
                  for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++)
                  {
                    int jTiming =  crystal[iCry].likelihood_mu_w[jLike].t_channel;
                    // size_t jValue = (size_t) jLike;
                    avgCTR += timeStamp[jTiming] - timeStamp[taggingCrystalTimingChannel];
                  }
                  avgCTR = avgCTR/crystal[iCry].likelihood_mu_w.size();
                  crystal[iCry].poliCorrCTR->Fill(avgCTR);


                  //calculate FloodZ...
                  Float_t FloodZ;
                  float centralChargeOriginal;
                  float centralSaturation;
                  float centralPedestal;
                  Float_t division = 0.0;

                  centralChargeOriginal = charge[crystal[iCry].detectorChannel];
                  for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                  {
                    if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
                    {
                      centralSaturation = detectorSaturation[iSat].saturation;
                      centralPedestal = detectorSaturation[iSat].pedestal;
                    }
                  }
                  float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );

                  for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
                  {
                    float originalCh = charge[crystal[iCry].relevantForW[iW]];
                    float saturationCh;
                    float pedestalCorr;
                    for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                    {
                      if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
                      {
                        saturationCh = detectorSaturation[iSat].saturation;
                        pedestalCorr = detectorSaturation[iSat].pedestal;
                      }
                    }
                    division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
                  }

                  FloodZ = centralChargeCorr / division;
                  // std::cout << FloodZ << std::endl;


                  for (unsigned int iW = 0 ; iW < crystal[iCry].w_limit_low.size() ; iW++) // run on w slices
                  {
                    if( (FloodZ > crystal[iCry].w_limit_low[iW]) && (FloodZ <= crystal[iCry].w_limit_up[iW])) // look for w slice for this w
                    {

                        // calc best time estimator using the appropriate inverse[iW]
                        double best_estimator = 0.0;
                        for (unsigned int jLike = 0 ; jLike < crystal[iCry].likelihood_mu_w.size(); jLike++)
                        {
                          // extract i timing channel
                          int jTiming =  crystal[iCry].likelihood_mu_w[jLike].t_channel;
                          size_t jValue = (size_t) jLike;


                          double mu_j = crystal[iCry].mu_i_w[jLike][iW];
                          double d_measure_j = timeStamp[jTiming] - timeStamp[taggingCrystalTimingChannel];

                          double weight = 0.0;
                          for (unsigned int iLike = 0 ; iLike < crystal[iCry].likelihood_mu_w.size(); iLike++)
                          {
                            size_t iValue = (size_t) iLike;
                            weight += gsl_matrix_get(crystal[iCry].inverse[iW], iValue, jValue);
                          }
                          best_estimator += weight * (d_measure_j - mu_j);
                        }
                        best_estimator = best_estimator / crystal[iCry].normalization_factor[iW];

                        crystal[iCry].likelihood_CTR->Fill(best_estimator);
                        crystal[iCry].scatter_likeT_w->Fill(FloodZ,best_estimator);

                        // std::cout << best_estimator << std::endl;

                    }
                  }
                }

              }
            }
          }
        }
        counter++;

        int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
        if( (perc % 10) == 0 )
        {
          std::cout << "\r";
          std::cout << perc << "% done... ";
          //std::cout << counter << std::endl;
        }
      }
      std::cout << "Good events = " << goodEvents << std::endl;




      std::sort(crystal.begin(), crystal.end(), compareByNumber);

      TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
      outputFile->cd();
      for(unsigned int i = 0 ;  i < crystal.size() ; i++)
      {
        crystal[i].simpleCTR->Write();
        crystal[i].poliCorrCTR->Write();
        crystal[i].likelihood_CTR->Write();
        crystal[i].scatter_likeT_w->Write();

        for (unsigned int iW = 0 ; iW < crystal[i].w_limit_low.size() ; iW++) // run on w slices
        {
          for(unsigned int iLike = 0 ; iLike < crystal[i].likelihood_mu_w.size(); iLike++) // run on detectors
          {
            crystal[i].mu_histo[iLike][iW]->Write();
          }
        }

        for(unsigned int iW = 0 ;iW < crystal[i].w_limit_low.size() ; iW++)
        {
          // crystal[i].h_correlation_matrix_element_w = new TH1F**[crystal[i].likelihood_mu_w.size()];
          for (unsigned int iLike = 0 ; iLike < crystal[i].likelihood_mu_w.size(); iLike++)
          {
            // crystal[i].h_correlation_matrix_element_w[iLike] = new TH1F*[crystal[i].likelihood_mu_w.size()];
            for (unsigned int jLike = 0 ; jLike < crystal[i].likelihood_mu_w.size(); jLike++)
            {
              // std::stringstream tname;
              // tname << "Cry_"<< crystal[i].number
                    // << "_histo_c_" << iLike+1 << "_" << jLike+1
                    // << "_w_" << iW
                    // << "_t_" << crystal[i].likelihood_mu_w[iLike].t_channel
                    // << "_t_" << crystal[i].likelihood_mu_w[jLike].t_channel;
              crystal[i].h_correlation_matrix_element_w[iW][iLike][jLike]->Write();
            }
          }
        }


      }
      outputFile->Close();

      calibrationFile->Close();
      // outputFile->Close();
      // textfile.close();
      // std::cout << std::endl;
      // std::cout << "Histograms saved in   " << outputFileName << std::endl;
      // std::cout << "Text summary saved in " << textFileName << std::endl;
      return 0;
    }
