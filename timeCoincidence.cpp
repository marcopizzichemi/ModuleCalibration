// compile with
// g++ -o ../build/timeCoincidence timeCoincidence.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

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
#include "TMatrixD.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <getopt.h>
#include <algorithm>    // std::sort
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "Crystal.h"
// #include "./include/ConfigFile.h"


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


//------------------------------------------//
//      FUNCTIONS INPUT PARSING             //
//------------------------------------------//
struct split_t
{
  enum empties_t { empties_ok, no_empties };
};

template <typename Container>
Container& split(Container& result,const typename Container::value_type& s,const typename Container::value_type& delimiters,split_t::empties_t empties = split_t::empties_ok)
{
  // splits the strings into individual fields
  // useful if you want to pass many parameters in the same string
  result.clear();
  size_t current;
  size_t next = -1;
  do
  {
    if (empties == split_t::no_empties)
    {
      next = s.find_first_not_of( delimiters, next + 1 );
      if (next == Container::value_type::npos) break;
      next -= 1;
    }
    current = next + 1;
    next = s.find_first_of( delimiters, current );
    result.push_back( s.substr( current, next - current ) );
  }
  while (next != Container::value_type::npos);
  return result;
}

void trim( std::string& s )
{
  // Remove leading and trailing whitespace
  static const char whitespace[] = " \n\t\v\r\f";
  s.erase( 0, s.find_first_not_of(whitespace) );
  s.erase( s.find_last_not_of(whitespace) + 1U );
}
//------------------------------------------//
//      end of FUNCTIONS INPUT PARSING      //
//------------------------------------------//




//------------------------------------------//
//      FUNCTIONS FOR FWHM CALCULATION      //
//------------------------------------------//

//FIXME both to be modified, no need for tagFwhm here!

void extractWithEMG(TH1F* histo,double fitPercMin,double fitPercMax,int divs, double tagFwhm,double* res,double* fitRes)
{
  // std::cout << "aaaaa" << std::endl;
  // preliminary gauss fit
  TCanvas *cTemp  = new TCanvas("temp","temp");
  TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  // resctrict the fitting range of gauss function

  gaussDummy->SetLineColor(kRed);
  double fitGaussMin = histo->GetMean()-2.0*histo->GetRMS();
  double fitGaussMax = histo->GetMean()+2.0*histo->GetRMS();
  double f1min = histo->GetXaxis()->GetXmin();
  double f1max = histo->GetXaxis()->GetXmax();
  if(fitGaussMin < f1min)
  {
    fitGaussMin = f1min;
  }
  if(fitGaussMax > f1max)
  {
    fitGaussMax = f1max;
  }
  TFitResultPtr rGauss = histo->Fit(gaussDummy,"QNS","",fitGaussMin,fitGaussMax);
  Int_t fitStatusGauss= rGauss;

  //NB fit results converted to int gives the fit status:
  // fitStatusGauss == 0 -> fit OK
  // fitStatusGauss != 0 -> fit FAILED

  double fitMin;
  double fitMax;
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    // use the histogram values
    fitMin = fitGaussMin;
    fitMax = fitGaussMax;
  }
  else
  {
    // use fit values
    fitMin = gaussDummy->GetParameter(1) - fitPercMin*(gaussDummy->GetParameter(2));
    fitMax = gaussDummy->GetParameter(1) + fitPercMax*(gaussDummy->GetParameter(2));
  }

  // chech that they are not outside the limits defined by user
  if(fitMin < f1min)
  {
    fitMin = f1min;
  }
  if(fitMax > f1max)
  {
    fitMax = f1max;
  }

  //define the standard gauss plus exp and the same but with inverted sign for x and mu
  //the standard definition should be equivalent to the classical definition of the Exponentially Modified Gaussian function
  TF1* gexp      = new TF1("gexp","[0]/sqrt(2)*exp([2]^2/2/[3]^2-(x-[1])/[3])*(1-TMath::Erf(([1]-x+[2]^2/[3])/(sqrt(2*[2]^2))))",f1min,f1max);
  gexp->SetLineColor(kGreen);
  gexp->SetParName(0,"N");
  gexp->SetParName(1,"mu");
  gexp->SetParName(2,"Sigma");
  gexp->SetParName(3,"tau");
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    gexp->SetParameters(histo->GetEntries(),histo->GetMean(),histo->GetRMS(),histo->GetRMS());
  }
  else
  {
    // use fit values
    gexp->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),gaussDummy->GetParameter(2));
  }
  TFitResultPtr r_gexp = histo->Fit(gexp,"QS","",fitMin,fitMax);

  TF1* gexp_inv  = new TF1("gexp_inv","[0]/sqrt(2)*exp([2]^2/2/[3]^2-((-x)-(-[1]))/[3])*(1-TMath::Erf(((-[1])-(-x)+[2]^2/[3])/(sqrt(2*[2]^2))))",f1min,f1max);
  gexp_inv->SetLineColor(kBlue);
  gexp_inv->SetParName(0,"N");
  gexp_inv->SetParName(1,"mu");
  gexp_inv->SetParName(2,"Sigma");
  gexp_inv->SetParName(3,"tau");
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    gexp_inv->SetParameters(histo->GetEntries(),histo->GetMean(),histo->GetRMS(),histo->GetRMS());
  }
  else
  {
    // use fit values
    gexp_inv->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),gaussDummy->GetParameter(2));
  }
  TFitResultPtr r_gexp_inv = histo->Fit(gexp_inv,"QS","",fitMin,fitMax);

  //try both, see what fits better

  Int_t fitStatusGexp = r_gexp;
  Int_t fitStatusGexp_inv = r_gexp_inv;

  double chi2gexp;
  double chi2gexp_inv;

  //NB fit results converted to int gives the fit status:
  // fitStatus == 0 -> fit OK
  // fitStatus != 0 -> fit FAILED

  if(fitStatusGexp == 0) // if gexp worked
  {
    chi2gexp = r_gexp->Chi2();
  }
  if(fitStatusGexp_inv == 0)// if gexp_inv worked
  {
    chi2gexp_inv   = r_gexp_inv->Chi2();
  }

  // now remember:
  // mean = mu + tau
  // sig = sqrt(Sigma^2 + tau^2)
  // and errors from error propagation
  // (meanErr)^2 = errMu^2 + errTau^2

  float mean = 0.0;
  float sigma = 0.0;
  float meanErr = 0.0;
  float sigmaErr = 0.0;
  TF1 *fitFunction = NULL;

  if(fitStatusGexp == 0 && fitStatusGexp_inv != 0) //if gexp worked and gexp_inv didn't worked
  {
    // std::cout << "fitStatusGexp     = " << fitStatusGexp << std::endl;
    // std::cout << "fitStatusGexp_inv = " << fitStatusGexp_inv << std::endl;
    // std::cout << "gexp worked and gexp_inv didn't worked" << std::endl;
    // fit again just to draw it
    // histo->Fit(gexp,"Q","",fitMin,fitMax);
    // delete the other function
    // delete gexp_inv;
    fitFunction = gexp;
  }
  else
  {
    if(fitStatusGexp != 0 && fitStatusGexp_inv == 0) //if gexp didn't work and gexp_inv worked
    {
      // delete gexp;
      // std::cout << "fitStatusGexp     = " << fitStatusGexp << std::endl;
      // std::cout << "fitStatusGexp_inv = " << fitStatusGexp_inv << std::endl;
      // std::cout << "gexp didn't work and gexp_inv worked" << std::endl;
      fitFunction = gexp_inv;
    }
    else // both worked or nothing worked
    {
      if(fitStatusGexp == 0 && fitStatusGexp_inv == 0)
      {
        // std::cout << "fitStatusGexp     = " << fitStatusGexp << std::endl;
        // std::cout << "fitStatusGexp_inv = " << fitStatusGexp_inv << std::endl;
        // std::cout << "both worked" << std::endl;
        // std::cout << "chi2gexp = " << chi2gexp << std::endl;
        // std::cout << "chi2gexp_inv = " << chi2gexp_inv << std::endl;
        if(chi2gexp > chi2gexp_inv) // gexp_inv better than gexp
        {
          // delete gexp;
          // std::cout << "chi2gexp > chi2gexp_inv" << std::endl;
          // std::cout << "using gexp_inv" << std::endl;
          fitFunction = gexp_inv;
        }
        else // gexp better than gexp_inv
        {
          // std::cout << "chi2gexp < chi2gexp_inv" << std::endl;
          // std::cout << "using gexp" << std::endl;
          // delete gexp_inv;
          fitFunction = gexp;
        }
      }
      else  // nothing worked
      {
        // leave values untouched
        // delete all func
        // std::cout << "fitStatusGexp     = " << fitStatusGexp << std::endl;
        // std::cout << "fitStatusGexp_inv = " << fitStatusGexp_inv << std::endl;
        // std::cout << "nothing worked" << std::endl;
        // try with gauss...
        TF1 *gaussCTR = new TF1("gaussCTR","gaus");
        gaussCTR->SetParameters(histo->GetEntries(),histo->GetMean(),histo->GetRMS());
        TFitResultPtr gCTR = histo->Fit(gaussCTR,"Q","",fitMin,fitMax); // re-fit just to store only the good one
        Int_t gRes = gCTR;
        if(gRes == 0)
        {
          fitFunction = gaussCTR;
          // std::cout << "Gauss fit worked" << std::endl;
        }
        else
        {
          fitFunction = NULL;
          // std::cout << "Not even gauss fit worked" << std::endl;
        }
        delete gexp_inv;
        delete gexp;
      }
    }
  }


  if(fitFunction == NULL)
  {
    res[0] = 0;
    res[1] = 0;
    res[2] = 0;
    res[3] = 0;
  }
  else
  {
    // SetParName(0,"N");
    // SetParName(1,"mu");
    // SetParName(2,"Sigma");
    // SetParName(3,"tau");
    //
    //
    histo->Fit(fitFunction,"Q","",fitMin,fitMax); // re-fit just to store only the good one
    // // write variables or it's gonna be a mess
    // float mu = fitFunction->GetParameter(1);
    // float e_mu = fitFunction->GetParError(1);
    // float s = fitFunction->GetParameter(2);
    // float e_s = fitFunction->GetParError(2);
    // float tau = fitFunction->GetParameter(3);
    // float e_tau = fitFunction->GetParError(3);
    //
    // mean = mu + tau;
    // sigma = TMath::Sqrt(TMath::Power(s,2) + TMath::Power(tau,2));
    // meanErr = TMath::Sqrt(TMath::Power(e_mu,2) + TMath::Power(e_tau,2));
    // sigmaErr = TMath::Sqrt(TMath::Power(s*e_s/sigma,2) + TMath::Power(tau*e_tau/sigma,2));
    //
    // res[0] = mean;
    // res[1] = sigma;
    // res[2] = meanErr;
    // res[3] = sigmaErr;


    double min,max,min10,max10;
    // int divs = 3000;
    double step = (f1max-f1min)/divs;
    double funcMax = fitFunction->GetMaximum(fitMin,fitMax);
    for(int i = 0 ; i < divs ; i++)
    {
      if( (fitFunction->Eval(f1min + i*step) < funcMax/2.0) && (fitFunction->Eval(f1min + (i+1)*step) > funcMax/2.0) )
      {
        min = f1min + (i+0.5)*step;
      }
      if( (fitFunction->Eval(f1min + i*step) > funcMax/2.0) && (fitFunction->Eval(f1min + (i+1)*step) < funcMax/2.0) )
      {
        max = f1min + (i+0.5)*step;
      }
      if( (fitFunction->Eval(f1min + i*step) < funcMax/10.0) && (fitFunction->Eval(f1min + (i+1)*step) > funcMax/10.0) )
      {
        min10 = f1min + (i+0.5)*step;
      }
      if( (fitFunction->Eval(f1min + i*step) > funcMax/10.0) && (fitFunction->Eval(f1min + (i+1)*step) < funcMax/10.0) )
      {
        max10 = f1min + (i+0.5)*step;
      }
    }
    res[0] = sqrt(2)*sqrt(pow((max-min),2)-pow(tagFwhm,2));
    res[1] = sqrt(2)*sqrt(pow((max10-min10),2)-pow((tagFwhm/2.355)*4.29,2));

    fitRes[0] = fitFunction->GetChisquare();
    fitRes[1] = fitFunction->GetNDF();
    fitRes[2] = fitFunction->GetProb();
  }
  delete cTemp;
}

void extractCTR(TH1F* histo,double fitPercMin,double fitPercMax, int divs, double tagFwhm, double* res, double* fitRes)
{

  // preliminary gauss fit
  TCanvas *cTemp  = new TCanvas("temp","temp");
  TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  // resctrict the fitting range of gauss function

  gaussDummy->SetLineColor(kRed);
  double fitGaussMin = histo->GetMean()-2.0*histo->GetRMS();
  double fitGaussMax = histo->GetMean()+2.0*histo->GetRMS();
  double f1min = histo->GetXaxis()->GetXmin();
  double f1max = histo->GetXaxis()->GetXmax();
  if(fitGaussMin < f1min)
  {
    fitGaussMin = f1min;
  }
  if(fitGaussMax > f1max)
  {
    fitGaussMax = f1max;
  }
  TFitResultPtr rGauss = histo->Fit(gaussDummy,"QNS","",fitGaussMin,fitGaussMax);
  Int_t fitStatusGauss= rGauss;

  //NB fit results converted to int gives the fit status:
  // fitStatusGauss == 0 -> fit OK
  // fitStatusGauss != 0 -> fit FAILED

  double fitMin;
  double fitMax;
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    // use the histogram values
    fitMin = fitGaussMin;
    fitMax = fitGaussMax;
  }
  else
  {
    // use fit values
    fitMin = gaussDummy->GetParameter(1) - fitPercMin*(gaussDummy->GetParameter(2));
    fitMax = gaussDummy->GetParameter(1) + fitPercMax*(gaussDummy->GetParameter(2));
  }

  // chech that they are not outside the limits defined by user
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
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    // use the histogram values
    cb->SetParameters(histo->GetEntries(),histo->GetMean(),histo->GetRMS(),1,3);
  }
  else
  {
    // use fit values
    cb->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  }
  TFitResultPtr rCb = histo->Fit(cb,"QNS","",fitMin,fitMax);

  //fit with gauss + exp
  TF1* gexp  = new TF1("gexp","[0]/sqrt(2)*exp([2]^2/2/[3]^2-(x-[1])/[3])*(1-TMath::Erf(([1]-x+[2]^2/[3])/(sqrt(2*[2]^2))))",f1min,f1max);
  gexp->SetLineColor(kGreen);
  gexp->SetParName(0,"N");
  gexp->SetParName(1,"Mean");
  gexp->SetParName(2,"Sigma");
  gexp->SetParName(3,"tau");
  // f1->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    // use the histogram values
    gexp->SetParameter(0,histo->GetEntries());
    gexp->SetParameter(1,histo->GetMean());
    gexp->SetParameter(2,histo->GetRMS());
    gexp->SetParameter(3,histo->GetRMS()); // ROOT really needs all parameters initialized, and a "good" guess for tau is the sigma of the previous fit...
  }
  else
  {
    // use fit values
    gexp->SetParameter(0,gaussDummy->GetParameter(0));
    gexp->SetParameter(1,gaussDummy->GetParameter(1));
    gexp->SetParameter(2,gaussDummy->GetParameter(2));
    gexp->SetParameter(3,gaussDummy->GetParameter(2)); // ROOT really needs all parameters initialized, and a "good" guess for tau is the sigma of the previous fit...
  }
  TFitResultPtr rGexp = histo->Fit(gexp,"QNS","",fitMin,fitMax);

  Int_t fitStatusCb = rCb;
  Int_t fitStatusGexp = rGexp;

  double chi2gexp;
  double chi2cb;

  if(fitStatusGexp == 0) // if Gexp worked
  {
    chi2gexp = rGexp->Chi2();
  }
  if(fitStatusCb == 0)// if cb worked
  {
    chi2cb   = rCb->Chi2();
  }


  //set function to measure ctr etc...
  TF1 *f1;
  if((fitStatusGexp  != 0) && (fitStatusCb != 0) && (fitStatusGauss != 0)) // all fit didn't work, just set everything to 0
  {
    res[0] = 0;
    res[1] = 0;

    fitRes[0] = 0;
    fitRes[1] = 0;
    fitRes[2] = 0;
    // res[2] = 0;
    // res[3] = 0;
  }
  else
  {
    if((fitStatusGexp  != 0) && (fitStatusCb != 0) && (fitStatusGauss == 0)) // only gauss worked
    {
      f1 = gaussDummy;
      f1->SetLineColor(kRed);
      histo->Fit(f1,"Q","",fitGaussMin,fitGaussMax);
      res[0] = sqrt(2)*sqrt(pow((2.355*f1->GetParameter(2)),2)-pow(tagFwhm,2));
      res[1] = sqrt(2)*sqrt(pow((4.29*f1->GetParameter(2)),2)-pow((tagFwhm/2.355)*4.29,2));

      fitRes[0] = f1->GetChisquare();
      fitRes[1] = f1->GetNDF();
      fitRes[2] = f1->GetProb();

      delete gexp;
      delete cb;
    }
    else
    {
      if((fitStatusGexp  != 0) && (fitStatusCb == 0)) // only cb worked
      {
        f1 = cb;
        f1->SetLineColor(kRed);
        histo->Fit(f1,"Q","",fitMin,fitMax);
        delete gexp;

      }
      else if((fitStatusGexp  == 0) && (fitStatusCb != 0)) // only gexp worked
      {
        f1 = gexp;
        f1->SetLineColor(kRed);
        histo->Fit(f1,"Q","",fitMin,fitMax);
        delete cb;
      }
      else // both worked
      {
        if(chi2gexp > chi2cb)
        {
          f1 = cb;
          f1->SetLineColor(kRed);
          histo->Fit(f1,"Q","",fitMin,fitMax);
          delete gexp;
        }
        else
        {
          f1 = gexp;
          f1->SetLineColor(kRed);
          histo->Fit(f1,"Q","",fitMin,fitMax);
          delete cb;
        }
      }

      double min,max,min10,max10;
      // int divs = 3000;
      double step = (f1max-f1min)/divs;
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
      res[0] = sqrt(2)*sqrt(pow((max-min),2)-pow(tagFwhm,2));
      res[1] = sqrt(2)*sqrt(pow((max10-min10),2)-pow((tagFwhm/2.355)*4.29,2));

      fitRes[0] = f1->GetChisquare();
      fitRes[1] = f1->GetNDF();
      fitRes[2] = f1->GetProb();
      // std::cout << f1->GetChisquare()/f1->GetNDF() << std::endl;
      delete cTemp;
    }
  }
}

//------------------------------------------//
//  end of FUNCTIONS FOR FWHM CALCULATION   //
//------------------------------------------//


bool compareByNumber(const Crystal_t &a,const Crystal_t  &b)
{
  return a.number < b.number;
}

// function to check if file exists
inline bool fileExists(const std::string& name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}


void usage()
{
  std::cout << "\t\t" << "[--folder] <path> [-i|--prefix] <file_prefix> [--calibration] <file1.root,file2.root,...> [OPTIONS]" << std::endl
            << "\t\t\t" << "<path>                                             - path to input files folder - default = \"./\" " << std::endl
            << "\t\t\t" << "<file_prefix>                                      - prefix of TTree files to analyze"   << std::endl
            << "\t\t\t" << "<file1.root,file2.root,...>                        - csv list of calibration files"   << std::endl
            << "\t\t" << "<output.root>                                      - output file name"   << std::endl
            << "\t\t" << "--histoMin <value>                                 - lower limit of CTR spectra, in sec - default = -15e-9"  << std::endl
            << "\t\t" << "--histoMax <value>                                 - upper limit of CTR spectra, in sec - default = 15e-9"  << std::endl
            << "\t\t" << "--histoBins <value>                                - n of bins for CTR spectra - default = 500"  << std::endl
            << "\t\t" << "--length <value>                                   - crystal length in mm, default = 15.0"  << std::endl
            << "\t\t" << "--doiFraction <value>                              - fraction of DOI length towards which the time stamps are corrected (from 0 to 1)"  << std::endl
            << "\t\t" << "                                                   - 0 = front of the crystal (DOI close to detector) "  << std::endl
            << "\t\t" << "                                                   - 1 = back of the crystal (DOI far from detector) "  << std::endl
            << "\t\t" << "--start-time                                       - acq time from which events are accepted [h]  - default = 0"  << std::endl
            << "\t\t" << "--sliced                                           - if given, it's a slice acq                   - default = not given"  << std::endl
            << std::endl;
}



void read_calibration(std::vector<Crystal_t> *crystal, std::vector<detector_t> *detectorSaturation, TFile* calibrationFile, TChain* tree,TList* formulas,Float_t histoMin , Float_t histoMax ,Int_t histoBins)
{
  // STANDARD CALIBRATION FILE

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
  float marginWZgraph = 0.1;
  float WrangeBinsForTiming = 1;
  std::string marginWZgraph_prefix("marginWZgraph");
  std::string WrangeBinsForTiming_prefix("WrangeBinsForTiming");
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
    detectorSaturation->push_back(tempDetector);
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
    if(!keysModName[i].compare(0,marginWZgraph_prefix.size(),marginWZgraph_prefix)) // find tcut
    {
      std::stringstream snameCh;
      snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
      marginWZgraph = atof(snameCh.str().c_str());
      std::cout << "margin cut set to " << marginWZgraph << std::endl;
    }
    if(!keysModName[i].compare(0,WrangeBinsForTiming_prefix.size(),WrangeBinsForTiming_prefix)) // find tcut
    {
      std::stringstream snameCh;
      snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
      // binningForWCut =
      // applyBinRestriction = true;
      // std::cout << "Applying bin restriction, binningForWCut set to " << binningForWCut << std::endl;
      WrangeBinsForTiming = atof(snameCh.str().c_str()); // hardcode it to this
    }

  }

  // std::stringstream sformulaname;
  // sformulaname << "FormulaTag";
  TCut taggingPhotopeakCutName;
  taggingPhotopeakCutName = taggingPhotopeakCut->GetTitle();
  // std::cout << "FormulaTag ------------- \n" << taggingPhotopeakCutName << std::endl;


  // TTreeFormula* FormulaTag = new TTreeFormula("FormulaTag",taggingPhotopeakCutName,tree);
  // formulas.Add(FormulaTag);

  // TTreeFormula* FormulaTagAnalysis = new TTreeFormula("FormulaTagAnalysis",taggingPhotopeakCutName,tree);
  // formulasAnalysis.Add(FormulaTagAnalysis);




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
       temp_crystal.minAcceptedW = 0;
       temp_crystal.maxAcceptedW = 1;
       temp_crystal.wMinSlicing = 0;
       temp_crystal.wMaxSlicing = 1;
       temp_crystal.wStepSlicing = 1;
       temp_crystal.marginWZgraph = marginWZgraph;
       temp_crystal.WrangeBinsForTiming = WrangeBinsForTiming;
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
       temp_crystal.likeCTR = NULL;
       temp_crystal.likeCTR_norm = NULL;
       temp_crystal.hybridCTR = NULL;
       temp_crystal.hybridCTR_norm = NULL;
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
         std::string t_channels_for_poli_mean_prefix("tChannelsForPolishedCorrectionMean");
         std::string t_channels_for_poli_FWHM_prefix("tChannelsForPolishedCorrectionFWHM");
         std::string mean_for_poli_prefix("meanForPolishedCorrection");
         std::string fwhm_for_poli_prefix("fwhmForPolishedCorrection");

         std::string lightCentral_prefix("Light collected in trigger crystal");
         std::string lightAll_prefix("Sum spectrum highlighted");
         std::string basicCTR_prefix("Basic CTR histogram");
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

           if(!keysCryName[i].compare(0,lightCentral_prefix.size(),lightCentral_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TH1F* aHisto = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
             temp_crystal.lightCentralHisto = aHisto;
           }

           if(!keysCryName[i].compare(0,lightAll_prefix.size(),lightAll_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TH1F* aHisto = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
             temp_crystal.lightAllHisto = aHisto;
           }

           if(!keysCryName[i].compare(0,basicCTR_prefix.size(),basicCTR_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TH1F* aHisto = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
             temp_crystal.basicCTRhisto = aHisto;
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

           if(!keysCryName[i].compare(0,t_channels_for_poli_mean_prefix.size(),t_channels_for_poli_mean_prefix))
           {
             std::vector<int> *v;
             gDirectory->GetObject("tChannelsForPolishedCorrectionMean",v);
             temp_crystal.polishedCorrection = true;
             temp_crystal.tChannelsForPolishedCorrectionMean = v[0];
           }

           if(!keysCryName[i].compare(0,t_channels_for_poli_FWHM_prefix.size(),t_channels_for_poli_FWHM_prefix))
           {
             std::vector<int> *v;
             gDirectory->GetObject("tChannelsForPolishedCorrectionFWHM",v);
             temp_crystal.polishedCorrection = true;
             temp_crystal.tChannelsForPolishedCorrectionFWHM = v[0];
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

           // if(likelihood)
           // {
           //   sname.str("");
           //   sname << "Likelihood correction - Crystal " << temp_crystal.number;
           //   temp_crystal.likeCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
           //   sname.str("");
           // }
           // if(hybridCorrection)
           // {
           //   sname.str("");
           //   sname << "Hybrid correction - Crystal " << temp_crystal.number;
           //   temp_crystal.hybridCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
           //   sname.str("");
           // }

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
                   std::size_t found = str2.find_first_of("_");                                  // find next "_"
                   std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                   // std::cout << keysTcorrName[i] << " " << str2 << " " << str3 << std::endl;     // output
                   int current_ch = atoi(str3.c_str());                                          // transform in int
                   // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output

                   bool acceptCh = true;
                   // for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                   // {
                   //   if(current_ch == forbidden_channels[iForb])
                   //   {
                   //     acceptCh = false;
                   //   }
                   // }

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
                   // for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                   // {
                   //   if(current_ch == forbidden_channels[iForb])
                   //   {
                   //     acceptCh = false;
                   //   }
                   // }

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

         bool dirDelayExists;
         sname.str("");
         dirDelayExists = gDirectory->cd("DelayDir");
         if(dirDelayExists)
         {
           TList *listTcorr = gDirectory->GetListOfKeys();
           int nKeysTcorr = listTcorr->GetEntries();
           std::vector<std::string> keysTcorrName;
           if(nKeysTcorr) //if directory not empty
           {
             for(int i = 0 ; i < nKeysTcorr ; i++){
               keysTcorrName.push_back(listTcorr->At(i)->GetName());
             }

             // std::string deltaWGraph_prefix = "DeltaW Graph";
             std::string graph_prefix = "Delay graph t";
             // std::string graph_delay_prefix = "Graph Delay ch_";
             // std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
             // std::string delay_timing_ch_prefix = "delayTimingChannels";
             for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
             {
               if(!keysTcorrName[i].compare(0,graph_prefix.size(),graph_prefix))
               {
                 TGraphErrors *calibGraph = NULL;
                 calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                 {
                   graphs_t temp_graph;
                   std::string rmsName = calibGraph->GetName();

                   std::size_t found = rmsName.find_last_of(graph_prefix);
                   std::string tChannelStringFromRms   = rmsName.substr(found  +1);

                   int ch   = atoi( tChannelStringFromRms.c_str() );

                   temp_graph.timingChannel = ch;
                   temp_graph.graph = calibGraph;

                   temp_crystal.delay_graphs.push_back(temp_graph);
                 }
               }
             }
           }
           gDirectory->cd("..");
         }



         bool dirRMSExists;
         sname.str("");
         dirRMSExists = gDirectory->cd("RMSdir");
         if(dirRMSExists)
         {
           TList *listTcorr = gDirectory->GetListOfKeys();
           int nKeysTcorr = listTcorr->GetEntries();
           std::vector<std::string> keysTcorrName;
           if(nKeysTcorr) //if directory not empty
           {
             for(int i = 0 ; i < nKeysTcorr ; i++){
               keysTcorrName.push_back(listTcorr->At(i)->GetName());
             }

             // std::string deltaWGraph_prefix = "DeltaW Graph";
             std::string graph_prefix = "RMS graph t";
             // std::string graph_delay_prefix = "Graph Delay ch_";
             // std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
             // std::string delay_timing_ch_prefix = "delayTimingChannels";
             for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
             {
               if(!keysTcorrName[i].compare(0,graph_prefix.size(),graph_prefix))
               {
                 TGraphErrors *calibGraph = NULL;
                 calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                 {
                   graphs_t temp_graph;
                   std::string rmsName = calibGraph->GetName();

                   std::size_t found = rmsName.find_last_of(graph_prefix);
                   std::string tChannelStringFromRms   = rmsName.substr(found  +1);

                   int rmsCh   = atoi( tChannelStringFromRms.c_str() );

                   temp_graph.timingChannel = rmsCh;
                   temp_graph.graph = calibGraph;

                   temp_crystal.rms_graphs.push_back(temp_graph);

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
         // globalCut += temp_crystal.CrystalCutWithoutCutG->GetTitle();     // this is BasicCut (XYZ and taggingPhotopeak) + CutTrigger (TriggerChannel and broadcut)
         globalCut += temp_crystal.PhotopeakEnergyCut->GetTitle();        // this is the cut on photopeak energy of the corrected spectrum for this crystal
         for(unsigned int iCutg = 0; iCutg < temp_crystal.cutg.size(); iCutg++)
         {
           globalCut += temp_crystal.cutg[iCutg]->GetName();              // these are the two cutg for this crystal
         }
         sname.str("");

         sname << "FormulaAnalysis" << temp_crystal.number;
         TTreeFormula* FormulaAnalysis = new TTreeFormula(sname.str().c_str(),globalCut,tree);
         formulas->Add(FormulaAnalysis);
         temp_crystal.FormulaAnalysis = FormulaAnalysis;
         sname.str("");

         if(temp_crystal.calibrationGraph && temp_crystal.FormulaAnalysis)
         {
           crystal->push_back(temp_crystal);
         }

       }
       gDirectory->cd("..");
    }
    calibrationFile->cd("Module 0.0");
  }
  // calibrationFile->Close();
}





//----------------//
//  MAIN PROGRAM  //
//----------------//
int main (int argc, char** argv)
{

  // check if there are args, otherwise print the usage info
  if(argc < 2)
  {
    std::cout << argv[0];
    usage();
    return 1;
  }

  // save the entire command line
  std::stringstream streamCommand;
  for(int i=0 ; i < argc; i++)
  {
    streamCommand << argv[i] << " ";
  }

  // set stat and fit information level in root files
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);


  //---------------------------------------//
  // PARSE COMMAND LINE INPUTs             //
  //---------------------------------------//

  std::string inputFileName = "";
  std::string analysis_folder = "./";
  std::string outputFileName = "outputFile.root";
  std::string calibrationFileNames = "";
  // std::string calibration_folder = "./";

  // std::string calibration_files = "";
  //
  //
  // // std::string coincidenceCalibrationFileName = "";
  //
  // std::string exclude_channels = "";
  // bool exclude = false;
  // bool simulation = false;
  Float_t length = 15.0; //mm
  Float_t doiFraction = 0.5;
  // Float_t tagFwhm = 70.0e-12; //s
  // Float_t rmsLow = 1.75;
  // Float_t rmsHigh = 1.75;
  Float_t histoMin = -15e-9;//s
  Float_t histoMax = 15e-9;//s
  // Float_t fitPercMin = 5;
  // Float_t fitPercMax = 6;
  // Float_t likeMin = -15e-9;//s
  // Float_t likeMax = 15e-9;//s
  // int likeBins = 500;
  // int divs       = 10000;
  int histoBins = 500;
  // int smooth = 0; //
  // int bins = 40;
  // double minCTR = 100;
  // double maxCTR = 500;
  // int func = 0;
  // bool unbinned = false;
  // bool fitCorrection = false;
  // bool basicLikelihood = false;
  // bool hybridCorrection = false;
  double start_time = 0;
  bool sliced = false;
  // bool likelihood = false;
  // bool likelihoodLine = false;
  // int WrangeBinsForTiming = 10;
  // float marginWZgraph = 0.1; // and then we read from modulecalib file
  // int binningForWCut = -1;
  // bool applyBinRestriction = false;
  // float marginWZgraph = 0.1;

  // parse arguments
  static struct option longOptions[] =
  {
			{ "prefix", required_argument, 0, 0 },
      { "folder", required_argument, 0, 0 },
      { "calibration", required_argument, 0, 0 },
      { "histoMin", required_argument, 0, 0 },
      { "histoMax", required_argument, 0, 0 },
      { "histoBins", required_argument, 0, 0 },
      { "length", required_argument, 0, 0 },
      { "doiFraction", required_argument, 0, 0 },
      { "start-time", required_argument, 0, 0 },
      { "sliced", no_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { NULL, 0, 0, 0 }

      // { "calibration", required_argument, 0, 0 },
      // { "simulation", no_argument, 0, 0 },

      // { "tagFwhm", required_argument, 0, 0 },
      // { "rmsLow", required_argument, 0, 0 },
      // { "rmsHigh", required_argument, 0, 0 },

      // { "smooth", required_argument, 0, 0 },
      // { "fitPercMin", required_argument, 0, 0 },
      // { "fitPercMax", required_argument, 0, 0 },
      // { "divs", required_argument, 0, 0 },
      // { "bins", required_argument, 0, 0 },
      // { "func", required_argument, 0, 0 },
      // { "unbinned", no_argument, 0, 0 },
      // { "fitCorrection", no_argument, 0, 0 },
      // { "exclude-channels", required_argument, 0, 0 },
      // { "start-time", required_argument, 0, 0 },
      // { "sliced", no_argument, 0, 0 },
      // { "likelihood", no_argument, 0, 0 },
      // { "likeMin", required_argument, 0, 0 },
      // { "likeMax", required_argument, 0, 0 },
      // { "likeBins", required_argument, 0, 0 },
      // { "basicLikelihood", no_argument, 0, 0 },
      // { "likelihoodLine", no_argument, 0, 0 },
      // { "hybridCorrection", no_argument, 0, 0 },
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			inputFileName = (char *)optarg;
    }
		else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      analysis_folder = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      calibrationFileNames = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 3){
      histoMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      histoMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 5){
      histoBins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6){
      length = atof((char *)optarg);;
    }
    else if (c == 0 && optionIndex == 7){
      doiFraction = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      start_time = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      sliced = true;
    }
    else if (c == 0 && optionIndex == 10){
      outputFileName = (char *)optarg;
    }
    else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
  }

  // else if (c == 0 && optionIndex == 2){
  //   calibrationFileName = (char *)optarg;
  // }
  // else if (c == 0 && optionIndex == 3){
  //   std::cout << "Dataset from simulation " << std::endl;
  //   simulation = true;
  // }

  // // else if (c == 0 && optionIndex == 6){
  // //   coincidenceCalibrationFileName = (char *)optarg;
  // // }
  // else if (c == 0 && optionIndex == 6){
  //   tagFwhm = atof((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 7){
  //   rmsLow = atof((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 8){
  //   rmsHigh = atof((char *)optarg);
  // }

  // else if (c == 0 && optionIndex == 12){
  //   smooth = atoi((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 13){
  //   fitPercMin = atof((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 14){
  //   fitPercMax = atof((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 15){
  //   divs = atoi((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 16){
  //   bins = atoi((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 17){
  //   func = atoi((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 18){
  //   unbinned = true;
  // }
  // else if (c == 0 && optionIndex == 19){
  //   fitCorrection = true;
  // }
  // else if (c == 0 && optionIndex == 20){
  //   exclude = true;
  //   exclude_channels = (char *)optarg;
  // }

  // else if (c == 0 && optionIndex == 23){
  //   likelihood = true;
  // }
  // else if (c == 0 && optionIndex == 24){
  //   likeMin = atof((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 25){
  //   likeMax = atof((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 26){
  //   likeBins = atoi((char *)optarg);
  // }
  // else if (c == 0 && optionIndex == 27){
  //   basicLikelihood = true;
  // }
  // else if (c == 0 && optionIndex == 28){
  //   likelihoodLine = true;
  // }
  // else if (c == 0 && optionIndex == 29){
  //   hybridCorrection = true;
  // }

  // check if required are given and files actually exists
  // first, input given and not empty
  if(inputFileName == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide the prefix of input files!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  if(calibrationFileNames == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide calibration files!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }

  // get input files list
  std::vector<std::string> v;
  read_directory(analysis_folder, v);
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;

  for(unsigned int i = 0 ; i < v.size() ; i++)
  {
    if(!v[i].compare(0,inputFileName.size(),inputFileName))
    {
      listInputFiles.push_back(analysis_folder + "/" + v[i]);
    }
  }
  // check if it's empty
  if(listInputFiles.size() == 0)
  {
    std::cout << std::endl;
    std::cout << "ERROR! Some input files do not exists! Aborting." << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }


  //calibration files
  std::vector<std::string> listCalibrationFiles;
  std::vector<std::string> calibrationFiles;
  split( listCalibrationFiles, calibrationFileNames, "," );  // split the entry
  bool calibrationFilesExist = true;
  for(unsigned int i = 0 ; i < listCalibrationFiles.size() ; i++)
  {
    if(!fileExists(listCalibrationFiles[i]))
    {
      calibrationFilesExist = false;
      std::cout << "ERROR! File " << listCalibrationFiles[i] << " does NOT exist!!!" << std::endl;
    }
  }

  if(calibrationFilesExist == false)
  {
    std::cout << std::endl;
    std::cout << "ERROR! Some input files do not exists! Aborting." << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }

  //---------------------------------------//
  // end of PARSE COMMAND LINE INPUTs      //
  //---------------------------------------//



  //---------------------------------------//
  // FEEDBACK PARAMETERS                   //
  //---------------------------------------//
  std::cout << std::endl;
  std::cout << "//-------------------------------------//"  << std::endl;
  std::cout << "// INPUT PARAMETERS                    //"  << std::endl;
  std::cout << "//-------------------------------------//"  << std::endl;
  std::cout << "Analysis folder      = " << analysis_folder << std::endl;
  std::cout << "File prefix          = " << inputFileName   << std::endl;
  std::cout << "Calibration files    = " ;
  for(unsigned int i = 0 ; i < listCalibrationFiles.size() ; i++)
  {
    std::cout << listCalibrationFiles[i];
    if(i < (listCalibrationFiles.size() -1))
    {
      std::cout << ",";
    }
  }
  std::cout << "histoMin             = " << histoMin      << std::endl;
  std::cout << "histoMax             = " << histoMax      << std::endl;
  std::cout << "histoBins            = " << histoBins     << std::endl;
  std::cout << "length               = " << length        << std::endl;
  std::cout << "DOI fraction         = " << doiFraction   << std::endl;
  std::cout << "Output file          = " << outputFileName << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  //---------------------------------------//
  // end of FEEDBACK PARAMETERS            //
  //---------------------------------------//


  // INPUT
  // READ ANALYSIS DATASET
  std::cout << std::endl;
  std::cout << "//-------------------------------------//" << std::endl;
  std::cout << "//         ANALYSIS FILES              //" << std::endl;
  std::cout << "//-------------------------------------//" << std::endl;

  std::cout << std::endl;

  //----------------------------------------------------------//
  //  Get TChain of analysis TTree files                   //
  //----------------------------------------------------------//

  TChain* tree = new TChain("adc");  // create the input tchain and the analysis ttree
  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    std::cout << "Adding file " << listInputFiles[i] << std::endl;
    tree->Add(listInputFiles[i].c_str());
  }
  std::cout << "|----------------------------------------|" << std::endl;
  std::cout << std::endl;

  // Add several TTreeFormula to the list;
  // TList formulasAnalysis;

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



  // INPUT
  // READ CALIBRATION DATASET
  std::cout << std::endl;
  std::cout << "//-------------------------------------//" << std::endl;
  std::cout << "//         CALIBRATION FILES           //" << std::endl;
  std::cout << "//-------------------------------------//" << std::endl;
  std::cout << std::endl;

  std::vector<Crystal_t> *pCrystal = new std::vector<Crystal_t>();
  std::vector<detector_t> *pDetectorSaturation = new std::vector<detector_t>();
  TList *formulas = new TList();
  std::vector<TFile*> calibrationFile;
  for(unsigned int i = 0 ; i < listCalibrationFiles.size() ; i++)
  {
    TFile* pCalibrationFile = new TFile(listCalibrationFiles[i].c_str());
    calibrationFile.push_back(pCalibrationFile);
  }

  for(unsigned int i = 0 ; i < listCalibrationFiles.size() ; i++)
  {
    TFile *calibrationFile = new TFile(listCalibrationFiles[i].c_str());
    read_calibration(pCrystal,pDetectorSaturation,calibrationFile,tree,formulas,histoMin,histoMax,histoBins);
  }

  std::vector<Crystal_t> crystal = *pCrystal;
  std::vector<detector_t> detectorSaturation = *pDetectorSaturation;
  std::cout << std::endl;
  std::cout << "|--------------------------------------------|" << std::endl;
  std::cout << "|                 w and z cuts               |" << std::endl;
  std::cout << "|--------------------------------------------|" << std::endl;
  for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  {
    if(crystal[iCry].accepted)
    {
      float beginW = crystal[iCry].wz->Eval(length - crystal[iCry].marginWZgraph);
      float endW = crystal[iCry].wz->Eval(crystal[iCry].marginWZgraph);


      // first wmean
      Float_t wmin,wmin2;
      Float_t wmax,wmax2;
      Float_t wmean,wmean2;

      int iBin = 0;

      wmin = beginW + ((iBin*(endW - beginW))/crystal[iCry].WrangeBinsForTiming);
      wmax = beginW + (((iBin+1)*(endW - beginW))/crystal[iCry].WrangeBinsForTiming);
      wmean = (wmax + wmin) / 2.0;
      crystal[iCry].wStepSlicing = (wmax - wmin);
      crystal[iCry].minAcceptedW = wmean;
      crystal[iCry].wMinSlicing = wmin;

      iBin = crystal[iCry].WrangeBinsForTiming -1;
      wmin2 = beginW + ((iBin*(endW - beginW))/crystal[iCry].WrangeBinsForTiming);
      wmax2 = beginW + (((iBin+1)*(endW - beginW))/crystal[iCry].WrangeBinsForTiming);
      wmean2 = (wmax2 + wmin2) / 2.0;
      crystal[iCry].maxAcceptedW = wmean2;
      crystal[iCry].wMaxSlicing = wmax2;


      //DEBUG
      std::cout << "Crystal = " << crystal[iCry].number << std::endl;
      std::cout << "Step lenght in w, for slicing = " << crystal[iCry].wStepSlicing << std::endl;
      std::cout << "Steps in w, for slicing = " << crystal[iCry].WrangeBinsForTiming << std::endl;
      std::cout << "min accepted w in slicing = " << crystal[iCry].wMinSlicing
                << " corresponding to " << crystal[iCry].calibrationGraph->Eval(crystal[iCry].wMinSlicing)
                << std::endl;
      std::cout << "min accepted w in Eval = " << crystal[iCry].minAcceptedW
                << " corresponding to " << crystal[iCry].calibrationGraph->Eval(crystal[iCry].minAcceptedW)
                << std::endl;
      std::cout << "max accepted w in slicing = " << crystal[iCry].wMaxSlicing
                << " corresponding to " << crystal[iCry].calibrationGraph->Eval(crystal[iCry].wMaxSlicing)
                << std::endl;
      std::cout << "max accepted w in Eval = " << crystal[iCry].maxAcceptedW
                << " corresponding to " << crystal[iCry].calibrationGraph->Eval(crystal[iCry].maxAcceptedW)
                << std::endl;

      //
      // fix the delay and rms graphs
      // rms are always one more than delay (because there is no delay for central channel, but there is a weight, the rms of the basic ctr), so run on rms vector
      for(unsigned int iRms = 0; iRms < crystal[iCry].rms_graphs.size(); iRms++ )
      {
        correction_graphs_t temp_corr;
        if(crystal[iCry].rms_graphs[iRms].timingChannel == crystal[iCry].timingChannel)
        {
          temp_corr.isMainChannel = true;
          temp_corr.delay = NULL;
        }
        else
        {
          temp_corr.isMainChannel = false;
          for(unsigned int iDel = 0; iDel < crystal[iCry].delay_graphs.size(); iDel++)
          {
            if(crystal[iCry].rms_graphs[iRms].timingChannel == crystal[iCry].delay_graphs[iDel].timingChannel)
            {
              temp_corr.delay = crystal[iCry].delay_graphs[iDel].graph;
            }

          }

        }
        temp_corr.timingChannel = crystal[iCry].rms_graphs[iRms].timingChannel;
        temp_corr.rms = crystal[iCry].rms_graphs[iRms].graph;
        crystal[iCry].correction_graphs.push_back(temp_corr);
      }
      // std::cout << crystal[iCry].delay_graphs.size() << std::endl;
      // std::cout << crystal[iCry].rms_graphs.size() << std::end
      std::cout << "Time correction graphs " << crystal[iCry].correction_graphs.size() << std::endl;

      // set the polished correction vectors
      // look in the mean vectors, see if there is correspondance in the rms one, if yes add to the vector
      for(unsigned int iPoli = 0; iPoli < crystal[iCry].tChannelsForPolishedCorrectionMean.size(); iPoli++)
      {
        for(unsigned int jPoli = 0; jPoli < crystal[iCry].tChannelsForPolishedCorrectionFWHM.size(); jPoli++)
        {
          if(crystal[iCry].tChannelsForPolishedCorrectionMean[iPoli] == crystal[iCry].tChannelsForPolishedCorrectionFWHM[jPoli])
          {
            polished_correction_t tempPoli;
            tempPoli.timingChannel = crystal[iCry].tChannelsForPolishedCorrectionMean[iPoli];
            tempPoli.mean          = crystal[iCry].meanForPolishedCorrection[iPoli];
            tempPoli.rms           = crystal[iCry].fwhmForPolishedCorrection[jPoli];
            crystal[iCry].polished_correction.push_back(tempPoli);
          }
        }
      }


    }
  }
  std::cout << "|--------------------------------------------|" << std::endl;
  std::cout << std::endl;

  for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  {
    if(crystal[iCry].accepted)
    {
      std::cout << "crystal " << crystal[iCry].number << std::endl;
      for(unsigned int iCor = 0; iCor < crystal[iCry].correction_graphs.size(); iCor++)
      {
        std::cout << crystal[iCry].correction_graphs[iCor].isMainChannel << " "
                  << crystal[iCry].correction_graphs[iCor].timingChannel << " ";
        if(crystal[iCry].correction_graphs[iCor].isMainChannel != true)
        {
          std::cout << crystal[iCry].correction_graphs[iCor].delay->GetName() << " ";
        }

        std::cout << crystal[iCry].correction_graphs[iCor].rms->GetName() << " "
                  << std::endl;
      }
      for(unsigned int iPoli = 0; iPoli < crystal[iCry].polished_correction.size(); iPoli++)
      {
        std::cout << crystal[iCry].polished_correction[iPoli].timingChannel << " "
                  << crystal[iCry].polished_correction[iPoli].mean << " "
                  << crystal[iCry].polished_correction[iPoli].rms << std::endl;
      }
    }
  }

  std::cout << "detectorSaturation = " <<  detectorSaturation.size() << std::endl;
  for(unsigned int iDet = 0 ;  iDet < detectorSaturation.size() ; iDet++)
  {
    std::cout<< detectorSaturation[iDet].digitizerChannel << std::endl;
  }


  // FIXME FOR NOW, restrict to hitogram of 2 crystals per run (no time to be more sofisticated...)
  std::stringstream sname;
  sname << "Basic Time Histogram";
  TH1F* basicCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
  sname.str("");
  sname << "Corrected Time Histogram";
  TH1F* correctedCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
  sname.str("");

  long long int counter = 0;

  // EVENTS LOOP

  tree->SetNotify(formulas);
  long long int nevent = tree->GetEntries();
  std::cout << "Total number of events in analysis file = " << nevent << std::endl;
  long int goodEvents = 0;

  double tStartAnalysis  = (tree->GetMinimum("ExtendedTimeTag"))/(1e9*3600); // t start in h
  double tEndAnalysis    = (tree->GetMaximum("ExtendedTimeTag") - tree->GetMinimum("ExtendedTimeTag"))/(1e9*3600); // t length in h

  double tStart2Analysis = (tree->GetMinimum("DeltaTimeTag"))/(1e9*3600);
  double tEnd2Analysis   = (tree->GetMaximum("DeltaTimeTag") - tree->GetMinimum("DeltaTimeTag"))/(1e9*3600);
  // for (long long int i=0;i<1000000;i++)
  for (long long int i=0;i<nevent;i++)
  {
    // std::cout << "Event " << i << std::endl;
    tree->GetEvent(i);              //read complete accepted event in memory
    //skip data if user say so
    bool keepEvent = true;
    if(sliced)
    {
      if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStartAnalysis) < start_time)
      {
        keepEvent = false;
      }
    }
    else
    {
      if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2Analysis) < start_time)
      {
        keepEvent = false;
      }
    }


    if(keepEvent)
    {
      for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
      {
        if(crystal[iCry].accepted)
        {
          if(crystal[iCry].FormulaAnalysis->EvalInstance())  //if in cutg and photopeak cuts for one crystal
          {

            //find w of this event
            //calculate FloodZ aka w
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
              // std::cout << crystal[iCry].relevantForW[iW] << std::endl;
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
              // std::cout << originalCh << " "
              //           << saturationCh << " "
              //           << pedestalCorr << " "
              //           << std::endl;
              division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
            }

            FloodZ = centralChargeCorr / division;



            // if one is accepted, run on the others, but then remember to stop!!!
            for(unsigned int jCry = 0 ;  jCry < crystal.size() ; jCry++)
            {
              if(iCry != jCry) // don't try to get a coincidence with the same crystal..
              {
                if(crystal[jCry].accepted)
                {
                  if(crystal[jCry].FormulaAnalysis->EvalInstance())  //if in cutg and photopeak cuts for another crystal
                  {


                    if((timeStamp[crystal[iCry].timingChannel] != 0) && (timeStamp[crystal[jCry].timingChannel] != 0)) // no zeroes
                    {
                      // basic CTR
                      double simpleCTR = timeStamp[crystal[iCry].timingChannel] - timeStamp[crystal[jCry].timingChannel];
                      basicCTR->Fill(simpleCTR);

                      // FIXME
                      // corrected CTR if both have correction graphs





                    }
                    goodEvents++;
                  }
                }
              }
            }
            // don't look for other coincidences (i.e. let's assume random coincidences are = 0 ...)
            break;
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
  std::cout << std::endl;
  std::cout << "Coincidence events found = " << goodEvents << std::endl;

  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  basicCTR->Write();
  correctedCTR->Write();
  outputFile->Close();
  std::cout << "Histrograms saved to file " << outputFileName << std::endl;




  return 0;
}
