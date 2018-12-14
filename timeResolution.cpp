// compile with
// g++ -o ../build/timeResolution timeResolution.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

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
#include <dirent.h>

#include "Crystal.h"

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


gsl_matrix *
invert_a_matrix(gsl_matrix *matrix,size_t size)
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
  std::cout << "\t\t" << "[--input-folder] <path> [-i|--input] <file_prefix>  [-o|--output] <output.root> [-c|--calibration] calibration.root [OPTIONS]" << std::endl
            << "\t\t" << "<path>                                            - path to the folder where  TTrees.root for the analysis are  - default = \"./\" " << std::endl
            << "\t\t" << "<file_prefix>                                      - prefix of TTree files to analyze"   << std::endl
            << "\t\t" << "<output.root>                                      - output file name"   << std::endl
            << "\t\t" << "<calibration.root>                                 - calibration file " << std::endl
            << "\t\t" << "--simulation                                       - the datast is from a simulation (therefore the tagging photopeak is ignored)" << std::endl
            << "\t\t" << "--length <value>                                   - crystal length in mm, default = 15.0"  << std::endl
            << "\t\t" << "--doiFraction <value>                              - fraction of DOI length towards which the time stamps are corrected (from 0 to 1)"  << std::endl
            << "\t\t" << "                                                   - 0 = front of the crystal (DOI close to detector) "  << std::endl
            << "\t\t" << "                                                   - 1 = back of the crystal (DOI far from detector) "  << std::endl
            << "\t\t" << "--tagFwhm <value>                                  - FWHM timing resolution of reference board, in sec - default = 70e-12"  << std::endl
            << "\t\t" << "--rmsLow <value>                                   - lower bound of CTR fit -> mean - rmsLow*mean - default = 1.75"  << std::endl
            << "\t\t" << "--rmsHigh <value>                                  - upper bound of CTR fit -> mean + rmsHigh*mean - default = 1.75"  << std::endl
            << "\t\t" << "--histoMin <value>                                 - lower limit of CTR spectra, in sec - default = -15e-9"  << std::endl
            << "\t\t" << "--histoMax <value>                                 - upper limit of CTR spectra, in sec - default = 15e-9"  << std::endl
            << "\t\t" << "--histoBins <value>                                - n of bins for CTR spectra - default = 500"  << std::endl
            << "\t\t" << "--smooth <value>                                   - n of iteration in CTR histograms smoothing - default = 0 (no smoothing)"  << std::endl
            << "\t\t" << "--fitPercMin <value>                               - time fit min is set to ((gauss fit mean) - fitPercMin*(gauss fit sigma))  - default = 5"  << std::endl
            << "\t\t" << "--fitPercMax <value>                               - time fit max is set to ((gauus fit mean) - fitPercMax*(gauss fit sigma))  - default = 6" << std::endl
            << "\t\t" << "--divs <value>                                     - n of divisions when looking for FWHM - default = 10000"  << std::endl
            << "\t\t" << "--bins <value>                                     - n of bins in summary CTR histograms - deafult 40"  << std::endl
            << "\t\t" << "--func <value>                                     - function for fitting (default = 0)"  << std::endl
            << "\t\t" << "                                                   - 0 = crystalball "  << std::endl
            << "\t\t" << "                                                   - 1 = gauss+exp "  << std::endl
            << "\t\t" << "--unbinned                                         - use also the unbinned method to calculate CTR - default = 0 (false)"  << std::endl
            << "\t\t" << "--fitCorrection                                    - use line fit to perform correction   - default = not given (false)"  << std::endl
            << "\t\t" << "--exclude-channels                                 - channels to exclude from time correction, comma separated - default = "" "  << std::endl
            << "\t\t" << "--start-time                                       - acq time from which events are accepted [h]  - default = 0"  << std::endl
            << "\t\t" << "--sliced                                           - if given, it's a slice acq                   - default = not given"  << std::endl
            << "\t\t" << "--likelihood                                       - if given, perform likelihood correction                   - default = not given"  << std::endl
            << "\t\t" << "--likeMin <value>                                  - lower limit of likelihood spectra, in sec - default = -5e-9"  << std::endl
            << "\t\t" << "--likeMax <value>                                  - upper limit of likelihood spectra, in sec - default = 5e-9"  << std::endl
            << "\t\t" << "--likeBins <value>                                 - n of bins for likelihood spectra - default = 500"  << std::endl
            << "\t\t" << "--basicLikelihood                                  - likelihood without line fits   - default = not given (false)"  << std::endl
            << "\t\t" << "--likelihoodLine                                   - using line fit of inverse covariance matrix, instead of tgraph eval. valid only if basicLikelihood is false - default = not given (false)"  << std::endl
            << "\t\t" << "--hybridCorrection                                 - performing hybrid correction - default = not given (false)"  << std::endl
            << "\t\t" << std::endl;
}

//----------------//
//  MAIN PROGRAM  //
//----------------//
int main (int argc, char** argv)
{
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }

  std::stringstream streamCommand;
  for(int i=0 ; i < argc; i++)
  {
    streamCommand << argv[i] << " ";
  }

  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);

  std::string inputFileName = "";
  std::string outputFileName = "";
  std::string calibrationFileName = "";
  std::string calibration_folder = "./";
  std::string analysis_folder = "./";
  std::string calibration_files = "";


  // std::string coincidenceCalibrationFileName = "";

  std::string exclude_channels = "";
  bool exclude = false;
  bool simulation = false;
  Float_t length = 15.0; //mm
  Float_t doiFraction = 0.5;
  Float_t tagFwhm = 70.0e-12; //s
  Float_t rmsLow = 1.75;
  Float_t rmsHigh = 1.75;
  Float_t histoMin = -15e-9;//s
  Float_t histoMax = 15e-9;//s
  Float_t fitPercMin = 5;
  Float_t fitPercMax = 6;
  Float_t likeMin = -15e-9;//s
  Float_t likeMax = 15e-9;//s
  int likeBins = 500;
  int divs       = 10000;
  int histoBins = 500;
  int smooth = 0; //
  int bins = 40;
  double minCTR = 100;
  double maxCTR = 500;
  int func = 0;
  bool unbinned = false;
  bool fitCorrection = false;
  bool basicLikelihood = false;
  bool hybridCorrection = false;
  double start_time = 0;
  bool sliced = false;
  bool likelihood = false;
  bool likelihoodLine = false;
  // int WrangeBinsForTiming = 10;
  // float marginWZgraph = 0.1; // and then we read from modulecalib file
  // int binningForWCut = -1;
  // bool applyBinRestriction = false;
  // float marginWZgraph = 0.1;

  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "calibration", required_argument, 0, 0 },
      { "simulation", no_argument, 0, 0 },
      { "length", required_argument, 0, 0 },
      { "doiFraction", required_argument, 0, 0 },
      { "tagFwhm", required_argument, 0, 0 },
      { "rmsLow", required_argument, 0, 0 },
      { "rmsHigh", required_argument, 0, 0 },
      { "histoMin", required_argument, 0, 0 },
      { "histoMax", required_argument, 0, 0 },
      { "histoBins", required_argument, 0, 0 },
      { "smooth", required_argument, 0, 0 },
      { "fitPercMin", required_argument, 0, 0 },
      { "fitPercMax", required_argument, 0, 0 },
      { "divs", required_argument, 0, 0 },
      { "bins", required_argument, 0, 0 },
      { "func", required_argument, 0, 0 },
      { "unbinned", no_argument, 0, 0 },
      { "fitCorrection", no_argument, 0, 0 },
      { "exclude-channels", required_argument, 0, 0 },
      { "start-time", required_argument, 0, 0 },
      { "sliced", no_argument, 0, 0 },
      { "likelihood", no_argument, 0, 0 },
      { "likeMin", required_argument, 0, 0 },
      { "likeMax", required_argument, 0, 0 },
      { "likeBins", required_argument, 0, 0 },
      { "basicLikelihood", no_argument, 0, 0 },
      { "likelihoodLine", no_argument, 0, 0 },
      { "hybridCorrection", no_argument, 0, 0 },
      { "input-folder", required_argument, 0, 0 },
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
    else if (c == 0 && optionIndex == 3){
      std::cout << "Dataset from simulation " << std::endl;
      simulation = true;
    }
    else if (c == 0 && optionIndex == 4){
      length = atof((char *)optarg);;
    }
    else if (c == 0 && optionIndex == 5){
      doiFraction = atof((char *)optarg);;
    }
    // else if (c == 0 && optionIndex == 6){
    //   coincidenceCalibrationFileName = (char *)optarg;
    // }
    else if (c == 0 && optionIndex == 6){
      tagFwhm = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      rmsLow = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      rmsHigh = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      histoMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      histoMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      histoBins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 12){
      smooth = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 13){
      fitPercMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 14){
      fitPercMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 15){
      divs = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 16){
      bins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 17){
      func = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 18){
      unbinned = true;
    }
    else if (c == 0 && optionIndex == 19){
      fitCorrection = true;
    }
    else if (c == 0 && optionIndex == 20){
      exclude = true;
      exclude_channels = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 21){
      start_time = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 22){
      sliced = true;
    }
    else if (c == 0 && optionIndex == 23){
      likelihood = true;
    }
    else if (c == 0 && optionIndex == 24){
      likeMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 25){
      likeMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 26){
      likeBins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 27){
      basicLikelihood = true;
    }
    else if (c == 0 && optionIndex == 28){
      likelihoodLine = true;
    }
    else if (c == 0 && optionIndex == 29){
      hybridCorrection = true;
    }
    else if (c == 0 && optionIndex == 30){
      analysis_folder = (char *)optarg;
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}

  if(calibration_files == "")
  {
    calibration_files = inputFileName;
  }

  std::cout << "Analysis folder       = " << analysis_folder << std::endl;
  std::cout << "Analysis file prefix  = " << inputFileName << std::endl;
  std::cout << "Calibration file      = " << calibrationFileName << std::endl;

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

  std::cout << "Chosen (length * doiFraction) = " << length * doiFraction << std::endl;
  if(fitCorrection)
  {
    std::cout << "Using linear fits to perform time correction" << std::endl;
  }
  if(likelihood)
  {
    std::cout << "Performing likelihood correction " << std::endl;
  }
  if(basicLikelihood)
  {
    std::cout << "Likelihood correction uses sliced arrays" << std::endl;
  }
  else
  {
    if(likelihoodLine)
    {
      std::cout << "Likelihood correction interpolates arrays with lines" << std::endl;
    }
    else
    {
      std::cout << "Likelihood correction interpolates arrays with TGraphs" << std::endl;
    }
  }
  if(hybridCorrection)
  {
    std::cout << "Performing hybrid correction" << std::endl;
  }

  //prepare output text file
  std::string textFileName = outputFileName.substr(0,outputFileName.size()-5);
  textFileName += ".txt";
  // std::cout << textFileName << std::endl;

  std::ofstream textfile;
  textfile.open (textFileName.c_str(),std::ofstream::out);




  // INPUT

  // ANALYSIS DATASET
  // read file in dir
  std::cout << std::endl;
  std::cout << "|----------------------------------------|" << std::endl;
  std::cout << "|         ANALYSIS FILES                 |" << std::endl;
  std::cout << "|----------------------------------------|" << std::endl;
  std::cout << std::endl;
  std::vector<std::string> v;
  read_directory(analysis_folder, v);
  // std::copy(v.begin(), v.end(),std::ostream_iterator<std::string>(std::cout, "\n"));
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;

  for(unsigned int i = 0 ; i < v.size() ; i++)
  {
    if(!v[i].compare(0,inputFileName.size(),inputFileName))
    {
      listInputFiles.push_back(analysis_folder + "/" + v[i]);
    }
  }

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
  TList formulasAnalysis;

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


  // supposedly the calibration and analysis ttrees have the same format, or you are in big trouble...

  // TObjArray *leavescopy = tree->GetListOfLeaves();
  // int nLeaves = leavescopy->GetEntries();
  // std::vector<std::string> leavesName;
  // // fill a vector with the leaves names
  // for(int i = 0 ; i < nLeaves ; i++)
  // {
  //   leavesName.push_back(leavescopy->At(i)->GetName());
  // }
  // // count the entries that start with "ch"
  // int numOfCh = 0;
  // // int numOfCry = 0;
  // std::string ch_prefix("ch");
  // std::string t_prefix("t");
  //
  // // std::string cry_prefix("cry");
  // for(int i = 0 ; i < nLeaves ; i++)
  // {
  //   //     leavesName.push_back(leavescopy->At(i)->GetName());
  //   if (!leavesName[i].compare(0, ch_prefix.size(), ch_prefix))
  //   {
  //     numOfCh++;
  //     detector_channels.push_back(atoi( (leavesName[i].substr(ch_prefix.size(),leavesName[i].size()-ch_prefix.size())).c_str() )) ;
  //   }
  //   // if (!leavesName[i].compare(0, cry_prefix.size(), cry_prefix))
  //   // numOfCry++;
  // }
  // //the string "cry" appears 4 times per crystal..
  // // numOfCry = numOfCry / 4;
  // std::cout << "Detector Channels \t= " << numOfCh << std::endl;
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
  // ULong64_t     ChainExtendedTimeTag_analysis;                                // extended time tag
  // ULong64_t     ChainDeltaTimeTag_analysis;                                   // delta tag from previous
  // UShort_t      *charge_analysis;
  // Float_t      *timeStamp_analysis;
  // TBranch      *bChainExtendedTimeTag_analysis;                               // branches for above data
  // TBranch      *bChainDeltaTimeTag_analysis;                                  // branches for above data
  // TBranch      **bCharge_analysis;
  // TBranch      **btimeStamp_analysis;
  //
  // charge_analysis = new UShort_t[numOfCh];
  // timeStamp_analysis = new Float_t[numOfCh];
  // bCharge_analysis = new TBranch*[numOfCh];
  // btimeStamp_analysis = new TBranch*[numOfCh];
  //
  // // set branches for reading the input files
  // tree->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag_analysis, &bChainExtendedTimeTag_analysis);
  // tree->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag_analysis, &bChainDeltaTimeTag_analysis);
  // // if(usingRealSimData)
  // // {
  // //   tree->SetBranchAddress("RealX", &RealX, &bRealX);
  // //   tree->SetBranchAddress("RealY", &RealY, &bRealY);
  // //   tree->SetBranchAddress("RealZ", &RealZ, &bRealZ);
  // //   // fchain->SetBranchAddress("Tagging", &simTaggingCharge, &bsimTaggingCharge);
  // //   // fchain->SetBranchAddress("TaggingTimeStamp", &simTaggingTime, &bsimTaggingTime);
  // //   tree->SetBranchAddress("CrystalsHit",&CrystalsHit, &bCrystalsHit);
  // //   tree->SetBranchAddress("NumbOfInteractions",&NumbOfInteractions, &bNumbOfInteractions);
  // //   // fchain->SetBranchAddress("TotalCryEnergy",&TotalCryEnergy, &bTotalCryEnergy);
  // // }
  // for (int i = 0 ; i < detector_channels.size() ; i++)
  // {
  //   //empty the stringstreams
  //   std::stringstream sname;
  //   sname << "ch" << detector_channels[i];
  //   // stype << "ch" << i << "/F";
  //   tree->SetBranchAddress(sname.str().c_str(),&charge_analysis[detector_channels[i]],&bCharge_analysis[detector_channels[i]]);
  //   sname.str("");
  //   // stype.str("");
  //
  //   sname << "t" << detector_channels[i];
  //   // stype << "t" << i << "/F";
  //   tree->SetBranchAddress(sname.str().c_str(),&timeStamp_analysis[detector_channels[i]],&btimeStamp_analysis[detector_channels[i]]);
  //   sname.str("");
  //   // stype.str("");
  // }

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
  float marginWZgraph = 0.1;
  float WrangeBinsForTiming = 0.1;
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
      // binningForWCut = atof(snameCh.str().c_str());
      // applyBinRestriction = true;
      // std::cout << "Applying bin restriction, binningForWCut set to " << binningForWCut << std::endl;
      WrangeBinsForTiming = atof(snameCh.str().c_str());
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

  // std::stringstream sformulaname;
  // sformulaname << "FormulaTag";
  TCut taggingPhotopeakCutName;
  taggingPhotopeakCutName = taggingPhotopeakCut->GetTitle();
  // std::cout << "FormulaTag ------------- \n" << taggingPhotopeakCutName << std::endl;


  // TTreeFormula* FormulaTag = new TTreeFormula("FormulaTag",taggingPhotopeakCutName,tree);
  // formulas.Add(FormulaTag);

  TTreeFormula* FormulaTagAnalysis = new TTreeFormula("FormulaTagAnalysis",taggingPhotopeakCutName,tree);
  formulasAnalysis.Add(FormulaTagAnalysis);




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

           if(likelihood)
           {
             sname.str("");
             sname << "Likelihood correction - Crystal " << temp_crystal.number;
             temp_crystal.likeCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
             sname.str("");
           }
           if(hybridCorrection)
           {
             sname.str("");
             sname << "Hybrid correction - Crystal " << temp_crystal.number;
             temp_crystal.hybridCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
             sname.str("");
           }

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
         globalCut += temp_crystal.CrystalCutWithoutCutG->GetTitle();     // this is BasicCut (XYZ and taggingPhotopeak) + CutTrigger (TriggerChannel and broadcut)
         globalCut += temp_crystal.PhotopeakEnergyCut->GetTitle();        // this is the cut on photopeak energy of the corrected spectrum for this crystal
         for(unsigned int iCutg = 0; iCutg < temp_crystal.cutg.size(); iCutg++)
         {
           globalCut += temp_crystal.cutg[iCutg]->GetName();              // these are the two cutg for this crystal
         }
         sname.str("");

         sname << "FormulaAnalysis" << temp_crystal.number;
         TTreeFormula* FormulaAnalysis = new TTreeFormula(sname.str().c_str(),globalCut,tree);
         formulasAnalysis.Add(FormulaAnalysis);
         temp_crystal.FormulaAnalysis = FormulaAnalysis;
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

  // //DEBUG
  // for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  // {
  //   std::cout << detectorSaturation[iSat].digitizerChannel << " " << detectorSaturation[iSat].saturation << " " << detectorSaturation[iSat].pedestal << std::endl;
  // }

  // COINCIDENCE CALIBRATION FILE
  //
  // TFile *coincidenceCalibrationFile = new TFile(coincidenceCalibrationFileName.c_str());
  // coincidenceCalibrationFile->cd("Module 0.0");
  // TList *listModuleCoinc = gDirectory->GetListOfKeys();
  // int nKeysModCoinc = listModuleCoinc->GetEntries();
  // std::vector<std::string> keysModNameCoinc;
  // // fill a vector with the leaves names
  // // std::string mppc_prefixCoinc("MPPC");
  // for(int i = 0 ; i < nKeysModCoinc ; i++){
  //   keysModNameCoinc.push_back(listModuleCoinc->At(i)->GetName());
  // }

//   TCut *taggingPhotopeakCut;
//   std::string taggingPhotopeakCut_prefix("taggingPhotopeakCut");
  // std::vector<std::string> MPPCfoldersCoinc;
//   for(unsigned int i = 0 ; i < keysModNameCoinc.size() ; i++)
//   {
    // if (!keysModNameCoinc[i].compare(0, mppc_prefixCoinc.size(), mppc_prefixCoinc))
    // {
    //   MPPCfoldersCoinc.push_back(keysModNameCoinc[i]);
    // }
//     if(!keysModNameCoinc[i].compare(0,taggingPhotopeakCut_prefix.size(),taggingPhotopeakCut_prefix)) // find tcut
//     {
     //  std::cout << keysCryName[i] << std::endl;
//       taggingPhotopeakCut = (TCut*) gDirectory->Get( keysModName[i].c_str());
      // if(cut)
      //   temp_crystal.CrystalCut = cut;
//     }
//   }

  // std::stringstream sformulaname;
  // sformulaname << "FormulaTag";
  // TCut taggingPhotopeakCutName;
  // taggingPhotopeakCutName = taggingPhotopeakCut->GetTitle();
  // TTreeFormula* FormulaTag = new TTreeFormula("FormulaTag",taggingPhotopeakCutName,tree);

  //look only into crystals found in the other calibration file
//   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
//   {
//     bool dirFound = gDirectory->cd(crystal[iCry].path);
//     if(dirFound)
//     {
//       crystal[iCry].accepted = true;
//       TList *listTcorr = gDirectory->GetListOfKeys();
//       int nKeysTcorr = listTcorr->GetEntries();
//       std::vector<std::string> keysTcorrName;
//       if(nKeysTcorr) //if directory not empty
//       {
//         for(int i = 0 ; i < nKeysTcorr ; i++){
//           keysTcorrName.push_back(listTcorr->At(i)->GetName());
//         }
//
//         std::string deltaWGraph_prefix = "DeltaW Graph";
//         std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";
//         // std::string graph_delay_prefix = "Graph Delay ch_";
//         // std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
//
//         for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
//         {
//           if(!keysTcorrName[i].compare(0,deltaWGraph_prefix.size(),deltaWGraph_prefix))
//           {
//             TGraph *calibGraph = NULL;
//             calibGraph = (TGraph*) gDirectory->Get(keysTcorrName[i].c_str());
//             if(calibGraph)
//               crystal[iCry].tw_correction = calibGraph;
//           }
//
//           if(!keysTcorrName[i].compare(0,rms_deltaWGraph_prefix.size(),rms_deltaWGraph_prefix))
//           {
//            TGraph *calibGraph = NULL;
//            calibGraph = (TGraph*) gDirectory->Get(keysTcorrName[i].c_str());
//            if(calibGraph)
//               crystal[iCry].rms_tw_correction = calibGraph;
//           }
//         }
//       }
//     }
//   }



  // // //BEGIN of DEBUG
  // for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  // {
  //   std::cout << crystal[i].number << "\t"
  //             // << crystal[i].cut->GetTitle() << "\t"
  //             << crystal[i].calibrationGraph->GetName() << "\t";
  //   for(unsigned int j = 0 ; j < crystal[i].cutg.size(); j++)
  //   {
  //     std::cout << crystal[i].cutg[j]->GetName() << "\t";
  //   }
  //   std::cout << std::endl;
  // }
  // // //END of DEBUG


  //------------------------------------------------//
  //------------------------------------------------//
  //------------------------------------------------//

  // the w of the crystal is sliced in modulecalibration, by doing
  // zMin = marginWZgraph           --> calc beginW using w(z) tgraph
  // zMax = length - marginWZgraph  --> calc endW using w(z) tgraph
  // then the range  (endW - beginW) is divided by WrangeBinsForTiming
  // so each bin is long (endW - beginW)/WrangeBinsForTiming
  // starts from
  // wmin = beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming)
  // and ends in
  // wmax = beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming)
  // then on this bin we calc whatever we are calculating, and the values are
  // associated to the mean w of the win, simply
  // wmean = (wmax + wmin) / 2.0
  // as a consequence all the tgraphs can be used only from the first wmean to the las wmean,
  // whatever is outside this range has to be thrown away
  // this has to be written in the crystal as a limit on the acceptance of w

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
      std::cout << "Correction graphs " << crystal[iCry].correction_graphs.size() << std::endl;

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



  // for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  // {
  //   if(crystal[iCry].accepted)
  //   {
  //     std::cout << "Crystal " << crystal[iCry].number << std::endl;
  //
  //     //count the delay plots, add 1, to get the relevant detector channels
  //     int nDetectors = crystal[iCry].delay.size() + 1;
  //     for(int iDet = 0 ; iDet < nDetectors ; iDet++)
  //     {
  //       ctr_aligned_t temp_ctr_aligned;
  //       int timingChannel;
  //       if(iDet == 0)
  //       {
  //         temp_ctr_aligned.isMainChannel = true;
  //         timingChannel = crystal[iCry].timingChannel;
  //         temp_ctr_aligned.delay_graph = NULL;
  //         temp_ctr_aligned.delay_rms_graph = NULL;
  //
  //         std::stringstream sname;
  //         sname << "Original_scatter - t" << timingChannel << "-t" << taggingCrystalTimingChannel << "vs.W_cry" << crystal[iCry].number;
  //         temp_ctr_aligned.original_scatter = new TH2F(sname.str().c_str(),sname.str().c_str(),1000,0,1,200,-1e-9,1e-9); //FIXME
  //         sname.str("");
  //         sname << "Aligned_scatter - t" << timingChannel << "-t" << taggingCrystalTimingChannel << "vs.W_cry" << crystal[iCry].number;
  //         temp_ctr_aligned.aligned_scatter  = new TH2F(sname.str().c_str(),sname.str().c_str(),1000,0,1,200,-1e-9,1e-9); //FIXME
  //         sname.str("");
  //       }
  //       else
  //       {
  //
  //         temp_ctr_aligned.isMainChannel = false;
  //
  //         std::string graphName = crystal[iCry].delay[iDet-1]->GetName();
  //         std::string rmsName = crystal[iCry].rms_delay[iDet-1]->GetName();
  //         std::size_t foundGraph = graphName.find_last_of("_");
  //         std::size_t foundRms   = rmsName.find_last_of("_");
  //
  //         std::string tChannelStringFromGraph = graphName.substr(foundGraph+1);
  //         std::string tChannelStringFromRms   = rmsName.substr(foundRms  +1);
  //         int graphCh = atoi(tChannelStringFromGraph.c_str() );
  //         int rmsCh   = atoi( tChannelStringFromRms.c_str() );
  //
  //         if(graphCh != rmsCh)
  //         {
  //           std::cout << "ERROR! TGraphs of delay and rms are from different timing channels!!!!" << std::endl;
  //           break;
  //         }
  //
  //         timingChannel = graphCh;
  //
  //         temp_ctr_aligned.delay_graph = crystal[iCry].delay[iDet-1];
  //         temp_ctr_aligned.delay_rms_graph = crystal[iCry].rms_delay[iDet-1];
  //
  //         std::stringstream sname;
  //         sname << "Original_scatter - t" << timingChannel << "-t" << taggingCrystalTimingChannel << "vs.W_cry" << crystal[iCry].number;
  //         temp_ctr_aligned.original_scatter = new TH2F(sname.str().c_str(),sname.str().c_str(),1000,0,1,200,-1e-9,4e-9); //FIXME
  //         sname.str("");
  //         sname << "Aligned_scatter - t" << timingChannel << "-t" << taggingCrystalTimingChannel << "vs.W_cry" << crystal[iCry].number;
  //         temp_ctr_aligned.aligned_scatter  = new TH2F(sname.str().c_str(),sname.str().c_str(),1000,0,1,200,-1e-9,4e-9); //FIXME
  //         sname.str("");
  //       }
  //
  //       temp_ctr_aligned.timingChannel = timingChannel;
  //       temp_ctr_aligned.ctr_aligned_graph = NULL;
  //       temp_ctr_aligned.ctr_aligned_rms_graph = NULL;
  //
  //       crystal[iCry].ctr_aligned.push_back(temp_ctr_aligned);
  //     }
  //   }
  // }



  std::cout << "Calibration data found for crystals: " << std::endl;
  for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  {
    if(crystal[i].accepted)
    {
      std::cout << crystal[i].number << std::endl;
    }
  }









  // ------------------------------------------- //
  // START OF LOOPS FOR COMPLETING CALIBRATION
  // ------------------------------------------- //

  //notify TTreeFormula(s) to TChain
  long long int nevent = tree->GetEntries();
  // ULong64_t tStart = tree->GetMinimum("ExtendedTimeTag");

  // std::cout << "Total number of events in calibration files = " << nevent << std::endl;
  long int goodEvents = 0;
  long int correlationMatrixEvents = 0;
  long int counter = 0;

  double tStart  = (tree->GetMinimum("ExtendedTimeTag"))/(1e9*3600); // t start in h
  double tEnd    = (tree->GetMaximum("ExtendedTimeTag") - tree->GetMinimum("ExtendedTimeTag"))/(1e9*3600); // t length in h

  double tStart2 = (tree->GetMinimum("DeltaTimeTag"))/(1e9*3600);
  double tEnd2   = (tree->GetMaximum("DeltaTimeTag") - tree->GetMinimum("DeltaTimeTag"))/(1e9*3600);

  if(likelihood)
  {
    std::cout << "LIKELIHOOD CORRECTION" <<  std::endl;
    //create the correlation matrix and its inverse in each accepted crystal
    std::cout << "Preparing arrays... " <<  std::endl;
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        std::cout << "Crystal " << crystal[iCry].number << std::endl;

        //count the delay plots, add 1, to get the likelihood channels
        int nLike = crystal[iCry].delay.size() + 1;


        float beginW = crystal[iCry].wz->Eval(length - crystal[iCry].marginWZgraph);
        float endW = crystal[iCry].wz->Eval(crystal[iCry].marginWZgraph);

        std::cout << "Preparing slices..." << std::endl;

        // since WrangeBinsForTiming is taken from the calibration file, there is no need to check here for the
        // min and maxAcceptedW
        for(int iBin = 0; iBin < crystal[iCry].WrangeBinsForTiming; iBin++) //
        {
          slice_t temp_slice;

          Float_t wmin = beginW + ((iBin*(endW - beginW))/crystal[iCry].WrangeBinsForTiming);
          Float_t wmax = beginW + (((iBin+1)*(endW - beginW))/crystal[iCry].WrangeBinsForTiming);
          Float_t wmean = (wmax + wmin) / 2.0;
          Float_t werr = (wmax-wmin)/TMath::Sqrt(12.0);

          temp_slice.wmean = wmean;
          temp_slice.werr = werr;
          temp_slice.wmin = wmin;
          temp_slice.wmax = wmax;
          temp_slice.entries = 0;
          temp_slice.covariance = new Float_t*[nLike];
          temp_slice.inverse_covariance = new Float_t*[nLike];
          temp_slice.entries_covariance = new long long int*[nLike];
          temp_slice.normalized_covariance = new Float_t*[nLike];

          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            temp_slice.covariance[iCorr] = new Float_t[nLike];
            temp_slice.inverse_covariance[iCorr] = new Float_t[nLike];
            temp_slice.entries_covariance[iCorr] = new long long int[nLike];
            temp_slice.normalized_covariance[iCorr] = new Float_t[nLike];
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              temp_slice.covariance[iCorr][jCorr] = 0;
              temp_slice.inverse_covariance[iCorr][jCorr] = 0;
              temp_slice.entries_covariance[iCorr][jCorr] = 0;
              temp_slice.normalized_covariance[iCorr][jCorr] = 0;
            }
            if(iCorr == 0) // crystal detector
            {
              temp_slice.tChannel.push_back(crystal[iCry].timingChannel);
              temp_slice.averageDelay.push_back(0); // no delay from cry to cry, by definition

            }
            else
            {
              std::string graphName = crystal[iCry].delay[iCorr-1]->GetName();
              std::size_t foundGraph = graphName.find_last_of("_");
              std::string tChannelStringFromGraph = graphName.substr(foundGraph+1);
              int graphCh = atoi(tChannelStringFromGraph.c_str() );
              temp_slice.tChannel.push_back(graphCh);
              temp_slice.averageDelay.push_back(crystal[iCry].delay[iCorr-1]->Eval(wmean));
            }
            temp_slice.averageDeltaT.push_back(0);
            temp_slice.varianceDeltaT.push_back(0); // inizialize delta T to 0
            temp_slice.nVarianceDeltaT.push_back(0); // inizialize delta T to 0n
            temp_slice.nDeltaT.push_back(0); //inizialize deltaT to 0
          }
          crystal[iCry].slice.push_back(temp_slice);
        }
      }
    }

    std::cout << "Calculating average delta t for each slice..." << std::endl;
    for (long long int i=0;i<nevent;i++)
    {
      tree->GetEvent(i);              //read complete accepted event in memory

      //skip data if user say so
      bool keepEvent = true;
      if(sliced)
      {
        if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStart) < start_time)
        {
          keepEvent = false;
        }
      }
      else
      {
        if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2) < start_time)
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
            if(FormulaTagAnalysis->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
            {
              if(crystal[iCry].FormulaAnalysis->EvalInstance())  //if in global cut of crystal
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

                // std::cout << FloodZ << std::endl;


                // no eval here, so no need for min and maxAcceptedW
                for(int iSlice = 0 ; iSlice <  crystal[iCry].slice.size() ; iSlice++)
                {
                  if((FloodZ >= crystal[iCry].slice[iSlice].wmin )&&( FloodZ < crystal[iCry].slice[iSlice].wmax ) )
                  {
                    // this is the slice
                    // run on all detectors
                    int nLike = crystal[iCry].slice[iSlice].tChannel.size();
                    for(int iCorr = 0; iCorr < nLike ; iCorr++)
                    {
                      if((timeStamp[crystal[iCry].slice[iSlice].tChannel[iCorr]] != 0) &&
                         (timeStamp[taggingCrystalTimingChannel] != 0))
                      {
                        Float_t deltaT = (timeStamp[crystal[iCry].slice[iSlice].tChannel[iCorr]] - timeStamp[taggingCrystalTimingChannel]) - crystal[iCry].slice[iSlice].averageDelay[iCorr];
                        crystal[iCry].slice[iSlice].averageDeltaT[iCorr] += deltaT;
                        crystal[iCry].slice[iSlice].nDeltaT[iCorr] += 1;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    std::cout << "Averaging the delta Ts" << std::endl;
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        for(int iSlice = 0 ; iSlice < crystal[iCry].WrangeBinsForTiming ; iSlice++)
        {
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            crystal[iCry].slice[iSlice].averageDeltaT[iCorr] = crystal[iCry].slice[iSlice].averageDeltaT[iCorr] / crystal[iCry].slice[iSlice].nDeltaT[iCorr];
          }
        }
      }
    }

    // //DEBUG
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "|----------------------------------------|" << std::endl;
    std::cout << "|        After averaging deltaTs         |" << std::endl;
    std::cout << "|----------------------------------------|" << std::endl;
    std::cout << std::endl;
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < crystal[iCry].WrangeBinsForTiming ; iSlice++)
        {
          std::cout << "slice = " << iSlice << std::endl;
          std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();

          std::cout << "covariance" << std::endl;

          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "inverse_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "entries_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "normalized_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            std::cout << crystal[iCry].slice[iSlice].tChannel[iCorr] << "\t"
                      << crystal[iCry].slice[iSlice].averageDelay[iCorr] << "\t"
                      << crystal[iCry].slice[iSlice].averageDeltaT[iCorr] << "\t"
                      << crystal[iCry].slice[iSlice].nDeltaT[iCorr] << "\t"
                      << crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] << "\t"
                      << crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr] << "\t"
                      << std::endl;
          }
        }
      }
    }
    //END OF DEBUG

// for(unsigned int jCorr = 0; jCorr < nLike ; jCorr++)
          // {
          //   temp_covariance[jCorr] = new Float_t[k];
          // }
          // for(unsigned int iCorr = 0; iCorr < k ; iCorr++)
          // {
          //   for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
          //   {
          //     temp_covariance[iCorr][jCorr] = 0;
          //   }
          // }
    std::cout << "Calculating covariance matrix..." << std::endl;

    for (long long int i=0;i<nevent;i++)
    {
      tree->GetEvent(i);              //read complete accepted event in memory
      //skip data if user say so
      bool keepEvent = true;
      if(sliced)
      {
        if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStart) < start_time)
        {
          keepEvent = false;
        }
      }
      else
      {
        if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2) < start_time)
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
            if(FormulaTagAnalysis->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
            {
              if(crystal[iCry].FormulaAnalysis->EvalInstance())  //if in global cut of crystal
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

                // std::cout << FloodZ << std::endl;



                for(int iSlice = 0 ; iSlice <  crystal[iCry].slice.size() ; iSlice++)
                {
                  if((FloodZ >= crystal[iCry].slice[iSlice].wmin )&&( FloodZ < crystal[iCry].slice[iSlice].wmax ) )
                  {
                    // this is the slice
                    // run on all detectors
                    int nLike = crystal[iCry].slice[iSlice].tChannel.size();
                    for(int iCorr = 0; iCorr < nLike ; iCorr++)
                    {
                      if((timeStamp[crystal[iCry].slice[iSlice].tChannel[iCorr]] != 0) &&
                         (timeStamp[taggingCrystalTimingChannel] != 0))
                      {

                        Float_t deltaT_i = (timeStamp[crystal[iCry].slice[iSlice].tChannel[iCorr]] - timeStamp[taggingCrystalTimingChannel]) - crystal[iCry].slice[iSlice].averageDelay[iCorr];
                        Float_t element_i = deltaT_i - crystal[iCry].slice[iSlice].averageDeltaT[iCorr];
                        //update variance
                        crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] += pow(element_i,2);
                        crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr]++;

                        for(int jCorr = 0; jCorr < nLike ; jCorr++)
                        {
                          if((timeStamp[crystal[iCry].slice[iSlice].tChannel[jCorr]] != 0) &&
                             (timeStamp[taggingCrystalTimingChannel] != 0))
                          {

                            Float_t deltaT_j = (timeStamp[crystal[iCry].slice[iSlice].tChannel[jCorr]] - timeStamp[taggingCrystalTimingChannel]) - crystal[iCry].slice[iSlice].averageDelay[jCorr];
                            Float_t element_j = deltaT_j - crystal[iCry].slice[iSlice].averageDeltaT[jCorr];
                            crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] += element_i * element_j;
                            crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr]++;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }


    std::cout << "Finalizing covariace matrix..." << std::endl;

    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        // std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < crystal[iCry].WrangeBinsForTiming ; iSlice++)
        {
          // std::cout << "slice = " << iSlice << std::endl;
          // std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          // std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          // std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          // std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          // std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();
          // std::cout << "covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] = crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] / crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr];
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] = crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] / crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr];
            }
          }
        }
      }
    }


    // //DEBUG
    // for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    // {
    //   if(crystal[iCry].accepted)
    //   {
    //     std::cout << "crystal " << crystal[iCry].number << std::endl;
    //     for(int iSlice = 0 ; iSlice < crystal[iCry].WrangeBinsForTiming ; iSlice++)
    //     {
    //       std::cout << "slice = " << iSlice << std::endl;
    //       std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
    //       std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
    //       std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
    //       std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
    //       std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
    //       int nLike = crystal[iCry].slice[iSlice].tChannel.size();
    //
    //       std::cout << "covariance" << std::endl;
    //
    //       for(int iCorr = 0; iCorr < nLike ; iCorr++)
    //       {
    //         for(int jCorr = 0; jCorr < nLike ; jCorr++)
    //         {
    //           std::cout << crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] << "\t";
    //         }
    //         std::cout << std::endl;
    //       }
    //
    //       std::cout << "inverse_covariance" << std::endl;
    //       for(int iCorr = 0; iCorr < nLike ; iCorr++)
    //       {
    //         for(int jCorr = 0; jCorr < nLike ; jCorr++)
    //         {
    //           std::cout << crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr] << "\t";
    //         }
    //         std::cout << std::endl;
    //       }
    //
    //       std::cout << "entries_covariance" << std::endl;
    //       for(int iCorr = 0; iCorr < nLike ; iCorr++)
    //       {
    //         for(int jCorr = 0; jCorr < nLike ; jCorr++)
    //         {
    //           std::cout << crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr] << "\t";
    //         }
    //         std::cout << std::endl;
    //       }
    //
    //       std::cout << "normalized_covariance" << std::endl;
    //       for(int iCorr = 0; iCorr < nLike ; iCorr++)
    //       {
    //         for(int jCorr = 0; jCorr < nLike ; jCorr++)
    //         {
    //           std::cout << crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] << "\t";
    //         }
    //         std::cout << std::endl;
    //       }
    //       for(int iCorr = 0; iCorr < nLike ; iCorr++)
    //       {
    //         std::cout << crystal[iCry].slice[iSlice].tChannel[iCorr] << "\t"
    //                   << crystal[iCry].slice[iSlice].averageDelay[iCorr] << "\t"
    //                   << crystal[iCry].slice[iSlice].averageDeltaT[iCorr] << "\t"
    //                   << crystal[iCry].slice[iSlice].nDeltaT[iCorr] << "\t"
    //                   << crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] << "\t"
    //                   << crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr] << "\t"
    //                   << std::endl;
    //       }
    //     }
    //   }
    // }
    // //END OF DEBUG

    std::cout << "Calculating normalized covariance matrix..." << std::endl;

    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        // std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < crystal[iCry].WrangeBinsForTiming ; iSlice++)
        {
          // std::cout << "slice = " << iSlice << std::endl;
          // std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          // std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          // std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          // std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          // std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();
          // std::cout << "covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            // crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] = crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] / crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr];
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] = crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] / (TMath::Sqrt( crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] * crystal[iCry].slice[iSlice].varianceDeltaT[jCorr]));
            }
          }
        }
      }
    }


    //DEBUG
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < crystal[iCry].WrangeBinsForTiming ; iSlice++)
        {
          std::cout << "slice = " << iSlice << std::endl;
          // std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          // std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          // std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          // std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          // std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();

          std::cout << "covariance" << std::endl;

          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "inverse_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "entries_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "normalized_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          // for(int iCorr = 0; iCorr < nLike ; iCorr++)
          // {
          //   std::cout << crystal[iCry].slice[iSlice].tChannel[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].averageDelay[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].averageDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].nDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr] << "\t"
          //             << std::endl;
          // }
        }
      }
    }
    //END OF DEBUG


    std::cout << "Inverting covariance matrix..." << std::endl;

    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        // std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < crystal[iCry].WrangeBinsForTiming ; iSlice++)
        {
          // std::cout << "slice = " << iSlice << std::endl;
          // std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          // std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          // std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          // std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          // std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();
          // use TMatrixF/D?
          // TMatrixD matrix(nLike,nLike); // never again, TMatrix simply doesn't work
          // use gls
          gsl_matrix *matrix = gsl_matrix_alloc(nLike, nLike);



          //copy covariance matrix to tmatrix
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            // crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] = crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] / crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr];
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              gsl_matrix_set(matrix,iCorr,jCorr,crystal[iCry].slice[iSlice].covariance[iCorr][jCorr]);
              // matrix[iCorr][jCorr] = crystal[iCry].slice[iSlice].covariance[iCorr][jCorr];
              // crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] = crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] / (TMath::Sqrt( crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] * crystal[iCry].slice[iSlice].varianceDeltaT[jCorr]));
            }
          }

          size_t size = (size_t) nLike;

          gsl_matrix *inverse_matrix = invert_a_matrix(matrix,size);

          //invert tmatrix
          // double det;
          // matrix.Invert(&det);

          //copy matrix into inverted covariance
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            // crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] = crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] / crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr];
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
               crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr] = gsl_matrix_get(inverse_matrix,iCorr,jCorr);
              // crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] = crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] / (TMath::Sqrt( crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] * crystal[iCry].slice[iSlice].varianceDeltaT[jCorr]));
            }
          }
        }
      }
    }


    //DEBUG
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < crystal[iCry].WrangeBinsForTiming ; iSlice++)
        {
          std::cout << "slice = " << iSlice << std::endl;
          // std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          // std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          // std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          // std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          // std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();

          std::cout << "covariance" << std::endl;

          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "inverse_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "entries_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "normalized_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          // for(int iCorr = 0; iCorr < nLike ; iCorr++)
          // {
          //   std::cout << crystal[iCry].slice[iSlice].tChannel[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].averageDelay[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].averageDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].nDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr] << "\t"
          //             << std::endl;
          // }
        }
      }
    }
    //END OF DEBUG


    std::cout << "Producing inverse covariance TGraphs..." << std::endl;

    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        //create the 2d matrix of tgraphs
        int nLike = crystal[iCry].slice[0].tChannel.size();
        int nSlice = crystal[iCry].slice.size();
        crystal[iCry].inverse_covariance_element = new TGraph**[nLike];
        crystal[iCry].inverse_covariance_element_line = new TF1**[nLike];
        for(int iCorr = 0; iCorr < nLike ; iCorr++)
        {
          crystal[iCry].inverse_covariance_element[iCorr] = new TGraph*[nLike];
          crystal[iCry].inverse_covariance_element_line[iCorr] = new TF1*[nLike];
          for(int jCorr = 0; jCorr < nLike ; jCorr++)
          {
            std::vector<Float_t> w;
            std::vector<Float_t> inv_s_ij;
            for(int iSlice = 0 ; iSlice < nSlice ; iSlice++)
            {
              w.push_back(crystal[iCry].slice[iSlice].wmean);
              inv_s_ij.push_back(crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr]);
            }
            crystal[iCry].inverse_covariance_element[iCorr][jCorr] = new TGraph(w.size(),&w[0],&inv_s_ij[0]);
            std::stringstream sname;
            sname << "cry_" << crystal[iCry].number <<"_inv_s_" << iCorr << "_" << jCorr;

            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->SetName(sname.str().c_str());
            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->SetTitle(sname.str().c_str());
            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->GetXaxis()->SetTitle("w");
            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->GetYaxis()->SetTitle("inv_cov element");

            crystal[iCry].inverse_covariance_element_line[iCorr][jCorr] = new TF1("line","[0]*x + [1]",0,1);
            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->Fit(crystal[iCry].inverse_covariance_element_line[iCorr][jCorr],"Q");

          }
        }

      }
    }
  } // end likelihood
















  // create the ctr_aligned 2d plots, then slice them to have a rms(w) estimation
  // std::cout << "Calculating ctr_aligned 2d plots..." << std::endl;
  //
  // for (long long int i=0;i<nevent;i++)
  // {
  //   tree->GetEvent(i);              //read complete accepted event in memory
  //   //skip data if user say so
  //
  //   bool keepEvent = true;
  //   if(sliced)
  //   {
  //     if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStart) < start_time)
  //     {
  //       keepEvent = false;
  //     }
  //   }
  //   else
  //   {
  //     if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2) < start_time)
  //     {
  //       keepEvent = false;
  //     }
  //   }
  //
  //   if(keepEvent)
  //   {
  //
  //     for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  //     {
  //       if(crystal[iCry].accepted)
  //       {
  //
  //         if(FormulaTag->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
  //         {
  //
  //           if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
  //           {
  //             std::cout << "enter" << std::endl;
  //
  //             //find w of this event
  //             //calculate FloodZ aka w
  //             Float_t FloodZ;
  //             float centralChargeOriginal;
  //             float centralSaturation;
  //             float centralPedestal;
  //             Float_t division = 0.0;
  //
  //             centralChargeOriginal = charge[crystal[iCry].detectorChannel];
  //             for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  //             {
  //               if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
  //               {
  //                 centralSaturation = detectorSaturation[iSat].saturation;
  //                 centralPedestal = detectorSaturation[iSat].pedestal;
  //               }
  //             }
  //             float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );
  //
  //             for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
  //             {
  //               // std::cout << crystal[iCry].relevantForW[iW] << std::endl;
  //               float originalCh = charge[crystal[iCry].relevantForW[iW]];
  //
  //               float saturationCh;
  //               float pedestalCorr;
  //               for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  //               {
  //                 if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
  //                 {
  //                   saturationCh = detectorSaturation[iSat].saturation;
  //                   pedestalCorr = detectorSaturation[iSat].pedestal;
  //                 }
  //               }
  //               // std::cout << originalCh << " "
  //               //           << saturationCh << " "
  //               //           << pedestalCorr << " "
  //               //           << std::endl;
  //               division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
  //             }
  //
  //             FloodZ = centralChargeCorr / division;
  //
  //             //in these scatter plots, accept events only if they come from the accepted range of w
  //             if(FloodZ > crystal[iCry].minAcceptedW && FloodZ < crystal[iCry].maxAcceptedW)
  //             {
  //               for(unsigned int iDet = 0; iDet < crystal[iCry].ctr_aligned.size(); iDet++ )
  //               {
  //                 int timingChannel = crystal[iCry].ctr_aligned[iDet].timingChannel;
  //                 if((timeStamp[timingChannel] != 0) &&
  //                 (timeStamp[taggingCrystalTimingChannel] != 0)) //only non zeroes
  //                 {
  //                   float delay;
  //                   if(crystal[iCry].ctr_aligned[iDet].isMainChannel) // main channel has no delay with itself...
  //                   {
  //                     delay = 0;
  //                   }
  //                   else
  //                   {
  //                     delay = crystal[iCry].ctr_aligned[iDet].delay_graph->Eval(FloodZ);
  //                   }
  //
  //                   crystal[iCry].ctr_aligned[iDet].original_scatter->Fill(FloodZ,timeStamp[timingChannel]-timeStamp[taggingCrystalTimingChannel]);
  //                   crystal[iCry].ctr_aligned[iDet].aligned_scatter->Fill(FloodZ,timeStamp[timingChannel]-timeStamp[taggingCrystalTimingChannel] - delay);
  //                 }
  //               }
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  //slice the aligned scatter to extract a tgraph per scatter and find rms(w) to use in weights
  // for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  // {
  //   if(crystal[iCry].accepted)
  //   {
  //     std::cout << "Crystal " << crystal[iCry].number << std::endl;
  //
  //     //count the delay plots, add 1, to get the relevant detector channels
  //     int nDetectors = crystal[iCry].delay.size() + 1;
  //     for(int iDet = 0 ; iDet < nDetectors ; iDet++) // slice everything, even main channel, just for simplicity of coding
  //     {
  //
  //       // float beginW = crystal[iCry].wz->Eval(length - marginWZgraph);
  //       // float endW = crystal[iCry].wz->Eval(marginWZgraph);
  //
  //       // std::cout << "Preparing slices..." << std::endl;
  //
  //       // since WrangeBinsForTiming is taken from the calibration file, there is no need to check here for the
  //       // min and maxAcceptedW
  //       for(int iBin = 0; iBin < WrangeBinsForTiming; iBin++) //
  //       {
  //         //TString name = crystal[iCry].ctr_aligned[iDet].aligned_scatter->GetName();
  //         float wmin = crystal[iCry].wMinSlicing + (iBin * crystal[iCry].wStepSlicing);
  //         float wmax = crystal[iCry].wMinSlicing + ((iBin+1) * crystal[iCry].wStepSlicing);
  //         int firstBin = crystal[iCry].ctr_aligned[iDet].aligned_scatter->GetXaxis()->FindBin(wmin);
  //         int lastBin = crystal[iCry].ctr_aligned[iDet].aligned_scatter->GetXaxis()->FindBin(wmax);
  //         // std::cout << iBin << " " << wmin << " " << wmax << " " << firstBin << " " << lastBin << std::endl;
  //
  //         TString name = crystal[iCry].ctr_aligned[iDet].aligned_scatter->GetName();
  //         std::stringstream sname;
  //         sname << name << "_" <<  firstBin << "_" << lastBin;
  //
  //         double fitPercMin = 5.0;
  //         double fitPercMax = 6.0;
  //         // int divisions = 10000;
  //         double res[4];
  //         TH1D *histo = crystal[iCry].ctr_aligned[iDet].aligned_scatter->ProjectionY(sname.str().c_str(),firstBin,lastBin);
  //         // std::cout << "debug 1" << std::endl;
  //
  //         extractWithEMG(histo,fitPercMin,fitPercMax,res);
  //         if(res[0] == 0 && res[1] == 0 && res[2] == 0 && res[3] == 0) // all fits failed, don't accept the channel
  //         {
  //
  //         }
  //         else
  //         {
  //           crystal[iCry].ctr_aligned[iDet].wmean.push_back((wmax+wmin)/2.0);
  //           crystal[iCry].ctr_aligned[iDet].werr.push_back(0.0);
  //           crystal[iCry].ctr_aligned[iDet].ctr_center.push_back(res[0]);
  //           crystal[iCry].ctr_aligned[iDet].ctr_err.push_back(res[2]);
  //           crystal[iCry].ctr_aligned[iDet].rms.push_back(res[1]);
  //           crystal[iCry].ctr_aligned[iDet].rms_err.push_back(res[3]);
  //         }
  //         crystal[iCry].ctr_aligned[iDet].slice.push_back(histo);
  //       }
  //
  //
  //
  //     }
  //
  //
  //
  //
  //   }
  // }
  //
  //
  // //slice the aligned scatter to extract a tgraph per scatter and find rms(w) to use in weights
  // for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  // {
  //   if(crystal[iCry].accepted)
  //   {
  //     std::cout << "Crystal " << crystal[iCry].number << std::endl;
  //
  //     //count the delay plots, add 1, to get the relevant detector channels
  //     int nDetectors = crystal[iCry].delay.size() + 1;
  //     for(int iDet = 0 ; iDet < nDetectors ; iDet++) // slice everything, even main channel, just for simplicity of coding
  //     {
  //       crystal[iCry].ctr_aligned[iDet].ctr_aligned_graph     = new TGraphErrors(crystal[iCry].ctr_aligned[iDet].wmean.size(),
  //                                                                                &crystal[iCry].ctr_aligned[iDet].wmean[0],
  //                                                                                &crystal[iCry].ctr_aligned[iDet].ctr_center[0],
  //                                                                                &crystal[iCry].ctr_aligned[iDet].werr[0],
  //                                                                                &crystal[iCry].ctr_aligned[iDet].ctr_err[0]);
  //       //
  //       std::stringstream sname;
  //       sname << "Ctr aligned slices - Detector " << iDet << "-t" << crystal[iCry].ctr_aligned[iDet].timingChannel <<  " - Crystal " << crystal[iCry].number;
  //       crystal[iCry].ctr_aligned[iDet].ctr_aligned_graph->SetTitle(sname.str().c_str());
  //       crystal[iCry].ctr_aligned[iDet].ctr_aligned_graph->SetName(sname.str().c_str());
  //       crystal[iCry].ctr_aligned[iDet].ctr_aligned_graph->GetXaxis()->SetTitle("W");
  //       crystal[iCry].ctr_aligned[iDet].ctr_aligned_graph->GetYaxis()->SetTitle("Time [S]");
  //       sname.str("");
  //
  //       crystal[iCry].ctr_aligned[iDet].ctr_aligned_rms_graph = new TGraphErrors(crystal[iCry].ctr_aligned[iDet].wmean.size(),
  //                                                                                &crystal[iCry].ctr_aligned[iDet].wmean[0],
  //                                                                                &crystal[iCry].ctr_aligned[iDet].rms[0],
  //                                                                                &crystal[iCry].ctr_aligned[iDet].werr[0],
  //                                                                                &crystal[iCry].ctr_aligned[iDet].rms_err[0]);
  //       //
  //       sname.str("");
  //       sname << "Ctr aligned slices rms - Crystal " << iDet << "-t" << crystal[iCry].ctr_aligned[iDet].timingChannel <<  " - Crystal " << crystal[iCry].number;
  //       crystal[iCry].ctr_aligned[iDet].ctr_aligned_rms_graph->SetTitle(sname.str().c_str());
  //       crystal[iCry].ctr_aligned[iDet].ctr_aligned_rms_graph->SetName(sname.str().c_str());
  //       crystal[iCry].ctr_aligned[iDet].ctr_aligned_rms_graph->GetXaxis()->SetTitle("W");
  //       crystal[iCry].ctr_aligned[iDet].ctr_aligned_rms_graph->GetYaxis()->SetTitle("Time [S]");
  //       sname.str("");
  //
  //     }
  //   }
  // }













  //   std::cout << "Calculating delay(doi) plots... ";
  //   //first, likelihood part if given
  //   for (long long int i=0;i<nevent;i++)
  //   {
  //     // std::cout << "Event " << i << std::endl;
  //     tree->GetEvent(i);              //read complete accepted event in memory
  //     //skip data if user say so
  //     bool keepEvent = true;
  //     if(sliced)
  //     {
  //       if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStart) < start_time)
  //       {
  //         keepEvent = false;
  //       }
  //     }
  //     else
  //     {
  //       if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2) < start_time)
  //       {
  //         keepEvent = false;
  //       }
  //     }
  //
  //
  //     if(keepEvent)
  //     {
  //       for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  //       {
  //         if(crystal[iCry].accepted)
  //         {
  //           if(FormulaTag->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
  //           {
  //             if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
  //             {
  //
  //               correlationMatrixEvents++;
  //
  //               if(crystal[iCry].tw_correction)
  //               {
  //
  //                 //calculate FloodZ aka w
  //                 Float_t FloodZ;
  //                 float centralChargeOriginal;
  //                 float centralSaturation;
  //                 float centralPedestal;
  //                 Float_t division = 0.0;
  //
  //                 centralChargeOriginal = charge[crystal[iCry].detectorChannel];
  //                 for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  //                 {
  //                   if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
  //                   {
  //                     centralSaturation = detectorSaturation[iSat].saturation;
  //                     centralPedestal = detectorSaturation[iSat].pedestal;
  //                   }
  //                 }
  //                 float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );
  //
  //                 for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
  //                 {
  //                   // std::cout << crystal[iCry].relevantForW[iW] << std::endl;
  //                   float originalCh = charge[crystal[iCry].relevantForW[iW]];
  //
  //                   float saturationCh;
  //                   float pedestalCorr;
  //                   for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  //                   {
  //                     if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
  //                     {
  //                       saturationCh = detectorSaturation[iSat].saturation;
  //                       pedestalCorr = detectorSaturation[iSat].pedestal;
  //                     }
  //                   }
  //                   // std::cout << originalCh << " "
  //                   //           << saturationCh << " "
  //                   //           << pedestalCorr << " "
  //                   //           << std::endl;
  //                   division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
  //                 }
  //
  //                 FloodZ = centralChargeCorr / division;
  //
  //                 // std::cout << FloodZ << std::endl;
  //
  //                 unsigned int k = crystal[iCry].LikelihoodChannels.size();
  //                 for(unsigned int iCorr = 0 ; iCorr < k; iCorr++)
  //                 {
  //
  //                   if( (timeStamp[crystal[iCry].LikelihoodChannels[iCorr].tChannel] != 0) &&
  //                       (timeStamp[taggingCrystalTimingChannel] != 0))
  //                   {
  //
  //                     Float_t iTimeStamp = timeStamp[crystal[iCry].LikelihoodChannels[iCorr].tChannel] - timeStamp[taggingCrystalTimingChannel];
  //                     // std::cout << iTimeStamp << std::endl;
  //                     Float_t iDelay = 0;
  //                     if(iCorr == 0)
  //                     {
  //                       iDelay = 0;
  //                     }
  //                     else
  //                     {
  //                       iDelay = crystal[iCry].LikelihoodChannels[iCorr].delay_line->Eval(FloodZ);
  //                     }
  //                     crystal[iCry].LikelihoodChannels[iCorr].deltaTscatter->Fill(FloodZ,iTimeStamp-iDelay);
  //                   }
  //                 }
  //               }
  //             }
  //           }
  //         }
  //       }
  //     }
  //     counter++;
  //
  //
  //   }
  //   // std::cout << "\r";
  //   std::cout << "done!" << std::endl;
  //   // std::cout << "Events for covariance matrix = " << correlationMatrixEvents << std::endl;
  //
  //
  //
  //   //slice the th2f plots
  //   std::cout << "Creating slice histrograms..." << std::endl;
  //   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++) // for each crystal
  //   {
  //     unsigned int k = crystal[iCry].LikelihoodChannels.size(); // for each 2d plot
  //     for(unsigned int iCorr = 0 ; iCorr < k; iCorr++)
  //     {
  //
  //       float beginW = crystal[iCry].wz->Eval(length - marginWZgraph);
  //       float endW = crystal[iCry].wz->Eval(marginWZgraph);
  //       for(int iBin = 0; iBin < WrangeBinsForTiming; iBin++) //
  //       {
  //         Float_t wmin = beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming);
  //         Float_t wmax = beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming);
  //         Float_t wmean = (wmax + wmin) / 2.0;
  //         std::stringstream sname;
  //         sname << "T" << crystal[iCry].LikelihoodChannels[iCorr].tChannel << "_like_slide_cry" <<  crystal[iCry].number << "_w_" << wmin << "_" << wmax;
  //         TH1F* temp_histo = new TH1F(sname.str().c_str(),sname.str().c_str(),likeBins,likeMin,likeMax);
  //         crystal[iCry].LikelihoodChannels[iCorr].deltaTslice.push_back(temp_histo);
  //         crystal[iCry].LikelihoodChannels[iCorr].wmean.push_back(wmean);
  //         crystal[iCry].LikelihoodChannels[iCorr].werr.push_back((wmax-wmin)/TMath::Sqrt(12.0));
  //         // crystal[iCry].LikelihoodChannels[iCorr].dx.push_back(wmean);
  //         // crystal[iCry].LikelihoodChannels[iCorr].dex.push_back((wmax-wmin)/TMath::Sqrt(12.0));
  //       }
  //     }
  //
  //
  //
  //   }
  //
  //   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++) // for each crystal
  //   {
  //     for(int iBin = 0; iBin < WrangeBinsForTiming; iBin++) //
  //     {
  //       unsigned int k = crystal[iCry].LikelihoodChannels.size(); // for each 2d plot
  //       Float_t **temp_covariance;
  //       temp_covariance = new Float_t*[k];
  //       for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //       {
  //         temp_covariance[jCorr] = new Float_t[k];
  //       }
  //       for(unsigned int iCorr = 0; iCorr < k ; iCorr++)
  //       {
  //         for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //         {
  //           temp_covariance[iCorr][jCorr] = 0;
  //         }
  //       }
  //       crystal[iCry].covariance.push_back(temp_covariance);
  //
  //       long long int **temp_entries_covariance;
  //       temp_entries_covariance = new long long int*[k];
  //       for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //       {
  //         temp_entries_covariance[jCorr] = new long long int[k];
  //       }
  //       for(unsigned int iCorr = 0; iCorr < k ; iCorr++)
  //       {
  //         for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //         {
  //           temp_entries_covariance[iCorr][jCorr] = 0;
  //         }
  //       }
  //       crystal[iCry].entries_covariance.push_back(temp_entries_covariance);
  //
  //       Float_t **inverse_covariance;
  //       inverse_covariance = new Float_t*[k];
  //       for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //       {
  //         inverse_covariance[jCorr] = new Float_t[k];
  //       }
  //       for(unsigned int iCorr = 0; iCorr < k ; iCorr++)
  //       {
  //         for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //         {
  //           inverse_covariance[iCorr][jCorr] = 0;
  //         }
  //       }
  //       crystal[iCry].inverse_covariance.push_back(inverse_covariance);
  //     }
  //   }
  //
  //
  //
  //   correlationMatrixEvents = 0;
  //   counter = 0;
  //   std::cout << "Filling slice histrograms..." << std::endl;
  //   for (long long int i=0;i<nevent;i++)
  //   {
  //     // std::cout << "Event " << i << std::endl;
  //     tree->GetEvent(i);              //read complete accepted event in memory
  //     //skip data if user say so
  //     bool keepEvent = true;
  //     if(sliced)
  //     {
  //       if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStart) < start_time)
  //       {
  //         keepEvent = false;
  //       }
  //     }
  //     else
  //     {
  //       if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2) < start_time)
  //       {
  //         keepEvent = false;
  //       }
  //     }
  //
  //
  //     if(keepEvent)
  //     {
  //       for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  //       {
  //         if(crystal[iCry].accepted)
  //         {
  //           if(FormulaTag->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
  //           {
  //             if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
  //             {
  //
  //               correlationMatrixEvents++;
  //
  //               if(crystal[iCry].tw_correction)
  //               {
  //
  //                 //calculate FloodZ aka w
  //                 Float_t FloodZ;
  //                 float centralChargeOriginal;
  //                 float centralSaturation;
  //                 float centralPedestal;
  //                 Float_t division = 0.0;
  //
  //                 centralChargeOriginal = charge[crystal[iCry].detectorChannel];
  //                 for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  //                 {
  //                   if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
  //                   {
  //                     centralSaturation = detectorSaturation[iSat].saturation;
  //                     centralPedestal = detectorSaturation[iSat].pedestal;
  //                   }
  //                 }
  //                 float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );
  //
  //                 for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
  //                 {
  //                   // std::cout << crystal[iCry].relevantForW[iW] << std::endl;
  //                   float originalCh = charge[crystal[iCry].relevantForW[iW]];
  //
  //                   float saturationCh;
  //                   float pedestalCorr;
  //                   for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  //                   {
  //                     if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
  //                     {
  //                       saturationCh = detectorSaturation[iSat].saturation;
  //                       pedestalCorr = detectorSaturation[iSat].pedestal;
  //                     }
  //                   }
  //                   // std::cout << originalCh << " "
  //                   //           << saturationCh << " "
  //                   //           << pedestalCorr << " "
  //                   //           << std::endl;
  //                   division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
  //                 }
  //
  //                 FloodZ = centralChargeCorr / division;
  //
  //
  //                 float beginW = crystal[iCry].wz->Eval(length - marginWZgraph);
  //                 float endW = crystal[iCry].wz->Eval(marginWZgraph);
  //                 for(int iBin = 0; iBin < WrangeBinsForTiming; iBin++) //
  //                 {
  //                   Float_t wmin = beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming);
  //                   Float_t wmax = beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming);
  //
  //                   if((FloodZ >=wmin) && (FloodZ < wmax) )
  //                   {
  //                     unsigned int k = crystal[iCry].LikelihoodChannels.size(); // for each 2d plot
  //                     for(unsigned int iCorr = 0 ; iCorr < k; iCorr++)
  //                     {
  //                       if( (timeStamp[crystal[iCry].LikelihoodChannels[iCorr].tChannel] != 0) &&
  //                       (timeStamp[taggingCrystalTimingChannel] != 0))
  //                       {
  //                         Float_t iTimeStamp = timeStamp[crystal[iCry].LikelihoodChannels[iCorr].tChannel] - timeStamp[taggingCrystalTimingChannel];
  //
  //                         Float_t iDelay = 0;
  //                         if(iCorr == 0)
  //                         {
  //                           iDelay = 0;
  //                         }
  //                         else
  //                         {
  //                           iDelay = crystal[iCry].LikelihoodChannels[iCorr].delay_line->Eval(FloodZ);
  //                         }
  //                         crystal[iCry].LikelihoodChannels[iCorr].deltaTslice[iBin]->Fill(iTimeStamp-iDelay);
  //                       }
  //                     }
  //                   }
  //                 }
  //               }
  //             }
  //           }
  //         }
  //       }
  //     }
  //     counter++;
  //
  //     // int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
  //     // if( (perc % 10) == 0 )
  //     // {
  //     //   std::cout << "\r";
  //     //   std::cout << "Calculating covariance matrix... " << perc << "%";
  //     //   //std::cout << counter << std::endl;
  //     // }
  //   }
  //   // std::cout << "\r";
  //   std::cout << "done!" << std::endl;
  //   // std::cout << "Events for covariance matrix = " << correlationMatrixEvents << std::endl;
  //
  //
  //   std::cout << "fitting slice histograms..." << std::endl;
  //   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++) // for each crystal
  //   {
  //     unsigned int k = crystal[iCry].LikelihoodChannels.size(); // for each 2d plot
  //     for(unsigned int iCorr = 0 ; iCorr < k; iCorr++)
  //     {
  //       for(unsigned int iSlice = 0 ; iSlice < crystal[iCry].LikelihoodChannels[iCorr].deltaTslice.size(); iSlice++ )
  //       {
  //         // crystal[iCry].LikelihoodChannels[iCorr].deltaTslice[iSlice];
  //         double res[4];
  //         int divisions = 100; //useless
  //         extractFromHisto(crystal[iCry].LikelihoodChannels[iCorr].deltaTslice[iSlice],fitPercMin,fitPercMax,divisions,res);
  //
  //         if(res[0] == 0 && res[1] == 0 && res[2] == 0 && res[3] == 0) //ignore point if fit didn't work
  //         {
  //                                 // skip point
  //         }
  //         else
  //         {
  //
  //           crystal[iCry].LikelihoodChannels[iCorr].dx.push_back(crystal[iCry].LikelihoodChannels[iCorr].wmean[iSlice]);
  //           crystal[iCry].LikelihoodChannels[iCorr].dex.push_back(crystal[iCry].LikelihoodChannels[iCorr].werr[iSlice]);
  //           // crystal[iCry].LikelihoodChannels[iCorr].dex.push_back(res[2]);
  //           crystal[iCry].LikelihoodChannels[iCorr].dy.push_back(res[0]);
  //           crystal[iCry].LikelihoodChannels[iCorr].dey.push_back(res[2]);
  //         }
  //       }
  //     }
  //   }
  //
  //
  //
  //   std::cout << "generating and fitting tgraphs..." << std::endl;
  //   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++) // for each crystal
  //   {
  //     unsigned int k = crystal[iCry].LikelihoodChannels.size(); // for each 2d plot
  //     for(unsigned int iCorr = 0 ; iCorr < k; iCorr++)
  //     {
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph = new TGraphErrors(
  //         crystal[iCry].LikelihoodChannels[iCorr].dx.size(),
  //         &crystal[iCry].LikelihoodChannels[iCorr].dx[0],
  //         &crystal[iCry].LikelihoodChannels[iCorr].dy[0],
  //         &crystal[iCry].LikelihoodChannels[iCorr].dex[0],
  //         &crystal[iCry].LikelihoodChannels[iCorr].dey[0]
  //       );
  //       std::stringstream sname;
  //       sname << "deltaT graph cry " <<  crystal[iCry].number << "_T" <<  crystal[iCry].LikelihoodChannels[iCorr].tChannel;
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph->SetName(sname.str().c_str());
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph->SetTitle(sname.str().c_str());
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph->GetXaxis()->SetTitle("W");
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph->GetYaxis()->SetTitle("DeltaT");
  //
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTline = new TF1("deltaTline",  "[0]*x + [1]",0,1);
  //       float m = (crystal[iCry].LikelihoodChannels[iCorr].dy[0]-
  //                 crystal[iCry].LikelihoodChannels[iCorr].dy[crystal[iCry].LikelihoodChannels[iCorr].dx.size()-1])/
  //                 (crystal[iCry].LikelihoodChannels[iCorr].dx[0]-
  //                 crystal[iCry].LikelihoodChannels[iCorr].dx[crystal[iCry].LikelihoodChannels[iCorr].dx.size()-1]);
  //       float q = crystal[iCry].LikelihoodChannels[iCorr].dy[0] -
  //                 m*crystal[iCry].LikelihoodChannels[iCorr].dx[0];
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTline->SetParameter(0,m);
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTline->SetParameter(1,q);
  //
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph->Fit(crystal[iCry].LikelihoodChannels[iCorr].deltaTline,"");
  //     }
  //   }
  //
  //   std::cout << "generating covariance scatter plots..." << std::endl;
  //
  //
  //
  //
  //

  counter = 0;










  // FINAL LOOP

  tree->SetNotify(&formulasAnalysis);
  long long int neventAnalysis = tree->GetEntries();
  // ULong64_t tStart = tree->GetMinimum("ExtendedTimeTag");

  std::cout << "Total number of events in analysis file = " << neventAnalysis << std::endl;
  long int goodEventsAnalysis = 0;
  // long int correlationMatrixEvents = 0;
  long int counterAnalysis = 0;

  double tStartAnalysis  = (tree->GetMinimum("ExtendedTimeTag"))/(1e9*3600); // t start in h
  double tEndAnalysis    = (tree->GetMaximum("ExtendedTimeTag") - tree->GetMinimum("ExtendedTimeTag"))/(1e9*3600); // t length in h

  double tStart2Analysis = (tree->GetMinimum("DeltaTimeTag"))/(1e9*3600);
  double tEnd2Analysis   = (tree->GetMaximum("DeltaTimeTag") - tree->GetMinimum("DeltaTimeTag"))/(1e9*3600);
  // for (long long int i=0;i<1000000;i++)
  for (long long int i=0;i<neventAnalysis;i++)
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
          if(FormulaTagAnalysis->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
          {
            if(crystal[iCry].FormulaAnalysis->EvalInstance())  //if in global cut of crystal
            {

              goodEventsAnalysis++;

              //temp commented
              Float_t centralcorrection = 0.0;
              Float_t zeroCorrection    = 0.0;
              //no corr
              double simpleCTR = timeStamp[crystal[iCry].timingChannel] - timeStamp[taggingCrystalTimingChannel];
              if((timeStamp[crystal[iCry].timingChannel] != 0) && (timeStamp[taggingCrystalTimingChannel] != 0)) // no zeroes
              {
                crystal[iCry].simpleCTR->Fill(simpleCTR);
                crystal[iCry].vSimple.push_back(simpleCTR);
              }

              if(crystal[iCry].tw_correction)
              {

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


                // float beginW = crystal[iCry].wz->Eval(length - marginWZgraph);
                // float endW = crystal[iCry].wz->Eval(marginWZgraph);

                //skip event if is cut by min and maxAcceptedW
                if(FloodZ > crystal[iCry].minAcceptedW && FloodZ < crystal[iCry].maxAcceptedW)
                {
                  double centralCTR;

                  if(fitCorrection)
                  {
                    centralcorrection = crystal[iCry].tw_correction_line->Eval(crystal[iCry].wz->Eval(length*doiFraction)) - crystal[iCry].tw_correction_line->Eval(FloodZ);
                  }
                  else
                  {
                    centralcorrection = crystal[iCry].tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)) - crystal[iCry].tw_correction->Eval(FloodZ);
                  }

                  centralCTR = (timeStamp[crystal[iCry].timingChannel] + (centralcorrection)) - timeStamp[taggingCrystalTimingChannel];
                  crystal[iCry].centralCTR->Fill(centralCTR);
                  crystal[iCry].vCentral.push_back(centralCTR);


                  if(crystal[iCry].delay.size())
                  {
                    // begin of new way
                    float averageTimeStamp = 0.0;
                    float totalWeight = 0.0;

                    //first quickly check if there are zeroes
                    bool noZeroes = true;
                    if(timeStamp[taggingCrystalTimingChannel] == 0)
                    {
                      noZeroes = false;
                    }
                    for(unsigned int iDet = 0; iDet < crystal[iCry].correction_graphs.size(); iDet++)
                    {
                      int timingChannel = crystal[iCry].correction_graphs[iDet].timingChannel;
                      if(timeStamp[timingChannel] == 0)
                      {
                        noZeroes = false;
                      }
                    }
                    for(unsigned int iDet = 0; iDet < crystal[iCry].correction_graphs.size(); iDet++)
                    {
                      int timingChannel = crystal[iCry].correction_graphs[iDet].timingChannel;
                      float delay;
                      if(crystal[iCry].correction_graphs[iDet].isMainChannel)
                      {
                        delay = 0;
                      }
                      else
                      {
                        delay = crystal[iCry].correction_graphs[iDet].delay->Eval(FloodZ);
                      }
                      float correctedElement = timeStamp[timingChannel] - timeStamp[taggingCrystalTimingChannel] - delay;
                      if(correctedElement <=  histoMin || correctedElement >= histoMax )
                      {

                        noZeroes = false;
                      }
                    }

                    if(noZeroes)
                    {
                      for(unsigned int iDet = 0; iDet < crystal[iCry].correction_graphs.size(); iDet++)
                      {
                        //run on all the detectors, included the main one, but remember not to correct the main one for delay!
                        int timingChannel = crystal[iCry].correction_graphs[iDet].timingChannel;
                        float delay;
                        if(crystal[iCry].correction_graphs[iDet].isMainChannel)
                        {
                          delay = 0;
                        }
                        else
                        {
                          delay = crystal[iCry].correction_graphs[iDet].delay->Eval(FloodZ);
                        }
                        float delta = timeStamp[timingChannel] - timeStamp[taggingCrystalTimingChannel] - delay;
                        float weight = (1.0)/( TMath::Power(crystal[iCry].correction_graphs[iDet].rms->Eval(FloodZ),2) );

                        totalWeight += weight;
                        averageTimeStamp += delta*weight;
                      }
                      averageTimeStamp = averageTimeStamp/totalWeight;
                      double allCTR = averageTimeStamp + centralcorrection;

                      crystal[iCry].allCTR->Fill(allCTR);



                      crystal[iCry].vAll.push_back(allCTR);

                      // // DEBUG of left tail
                      // if(averageTimeStamp < 0.2e-9)
                      // {
                      //   std::cout << i << " " << FloodZ << " ";
                      //   std::cout << timeStamp[taggingCrystalTimingChannel]  << " ";
                      //   for(unsigned int iDet = 0; iDet < crystal[iCry].correction_graphs.size(); iDet++)
                      //   {
                      //     int timingChannel = crystal[iCry].correction_graphs[iDet].timingChannel;
                      //     std::cout << timeStamp[timingChannel] << " ";
                      //   }
                      //
                      //   std::cout << std::endl;
                      // }
                      // // end of DEBUG of left tail
                    }

                    //
                    // if(noZeroes)
                    // {
                    //
                    //   for(unsigned int iDet = 0; iDet < crystal[iCry].ctr_aligned.size(); iDet++)
                    //   {
                    //     //run on all the detectors, included the main one, but remember not to correct the main one for delay!
                    //     int timingChannel = crystal[iCry].ctr_aligned[iDet].timingChannel;
                    //     float delay;
                    //     if(crystal[iCry].ctr_aligned[iDet].isMainChannel)
                    //     {
                    //       delay = 0;
                    //     }
                    //     else
                    //     {
                    //       delay = crystal[iCry].ctr_aligned[iDet].delay_graph->Eval(FloodZ);
                    //     }
                    //     float delta = timeStamp_analysis[timingChannel] - timeStamp_analysis[taggingCrystalTimingChannel] - delay;
                    //     float weight = (1.0)/( TMath::Power(crystal[iCry].ctr_aligned[iDet].ctr_aligned_rms_graph->Eval(FloodZ),2) );
                    //
                    //     // std::cout << delta << " " << weight << std::endl;
                    //
                    //     totalWeight += weight;
                    //     averageTimeStamp += delta*weight;
                    //   }
                    //   averageTimeStamp = averageTimeStamp/totalWeight;
                    //   double allCTR = averageTimeStamp + centralcorrection;
                    //   crystal[iCry].allCTR->Fill(allCTR);
                    //   crystal[iCry].vAll.push_back(allCTR);
                    // }
                    // //--- end of new way









                    if(likelihood)
                    {
                      //combine the detector results using inverse covariance matrix elements
                      // calc the sum of all elements

                      //first check the no timestamp is 0
                      bool noZeroes = true;
                      if(timeStamp[taggingCrystalTimingChannel] == 0)
                      {
                        noZeroes = false;   // and you can already stop
                      }
                      else // look for all the other channels involved
                      {
                        for(unsigned int iCorr = 0; iCorr < crystal[iCry].slice[0].tChannel.size() ; iCorr++)
                        {
                          if(timeStamp[crystal[iCry].slice[0].tChannel[iCorr]] == 0)
                          {
                            noZeroes = false;
                          }
                        }
                      }

                      if(noZeroes)
                      {
                        if(basicLikelihood)
                        {
                          // zero approach: just find the w slice and use that inv cov matrix
                          int nSlice = crystal[iCry].slice.size();
                          for(int iSlice = 0 ; iSlice < nSlice ; iSlice++)
                          {
                            if((FloodZ >= crystal[iCry].slice[iSlice].wmin )&&( FloodZ < crystal[iCry].slice[iSlice].wmax ) ) // this is the w slice
                            {

                                Float_t sum_inv_cov = 0;
                                Float_t Dt_best = 0;
                                int nLike = crystal[iCry].slice[iSlice].tChannel.size();
                                for(int jCorr = 0; jCorr < nLike ; jCorr++)
                                {
                                  int tChannel = crystal[iCry].slice[iSlice].tChannel[jCorr];
                                  // find delay as function of doi
                                  Float_t delay;
                                  if(jCorr == 0)
                                  {
                                    delay = 0; // no delay for direct detector
                                  }
                                  else
                                  {
                                    delay = crystal[iCry].delay_line[jCorr-1]->Eval(FloodZ);
                                  }
                                  // find Dt of this detector (t - tR) - delay;
                                  Float_t Dt_j = (timeStamp[tChannel] - timeStamp[taggingCrystalTimingChannel]) - delay;

                                  // now calc the weight_i
                                  // keeping j fixed, run on i and sum the inv cov elements
                                  Float_t weight_i = 0;
                                  for(int iCorr = 0; iCorr < nLike ; iCorr++)
                                  {
                                    weight_i += crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr];
                                    // weight_i +=  crystal[iCry].inverse_covariance_element[iCorr][jCorr]
                                    // at the same time also sum all the elements of inv cov matrix to get the final division
                                    // sum_inv_cov += crystal[iCry].inverse_covariance_element[iCorr][jCorr]->Eval(FloodZ);
                                    sum_inv_cov += crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr];
                                  }
                                  Dt_best += weight_i * Dt_j;
                                }

                                Dt_best = Dt_best / sum_inv_cov;
                                double likeCTR = Dt_best + centralcorrection ;
                                crystal[iCry].likeCTR->Fill(likeCTR);
                                crystal[iCry].vLike.push_back(likeCTR);
                            }
                          }
                        }
                        else
                        {
                          if(likelihoodLine)
                          {
                            Float_t sum_inv_cov = 0;
                            Float_t Dt_best = 0;
                            int nLike = crystal[iCry].slice[0].tChannel.size();
                            for(int jCorr = 0; jCorr < nLike ; jCorr++)
                            {
                              int tChannel = crystal[iCry].slice[0].tChannel[jCorr];
                              // find delay as function of doi
                              Float_t delay;
                              if(jCorr == 0)
                              {
                                delay = 0; // no delay for direct detector
                              }
                              else
                              {
                                delay = crystal[iCry].delay_line[jCorr-1]->Eval(FloodZ);
                              }
                              // find Dt of this detector (t - tR) - delay;
                              Float_t Dt_j = (timeStamp[tChannel] - timeStamp[taggingCrystalTimingChannel]) - delay;

                              // now calc the weight_i
                              // keeping j fixed, run on i and sum the inv cov elements
                              Float_t weight_i = 0;
                              for(int iCorr = 0; iCorr < nLike ; iCorr++)
                              {
                                weight_i += crystal[iCry].inverse_covariance_element_line[iCorr][jCorr]->Eval(FloodZ);
                                sum_inv_cov += crystal[iCry].inverse_covariance_element_line[iCorr][jCorr]->Eval(FloodZ);
                              }
                              Dt_best += weight_i * Dt_j;
                            }
                            Dt_best = Dt_best / sum_inv_cov;
                            double likeCTR = Dt_best + centralcorrection ;
                            crystal[iCry].likeCTR->Fill(likeCTR);
                            crystal[iCry].vLike.push_back(likeCTR);
                          }
                          else
                          {
                            Float_t sum_inv_cov = 0;
                            Float_t Dt_best = 0;
                            int nLike = crystal[iCry].slice[0].tChannel.size();
                            for(int jCorr = 0; jCorr < nLike ; jCorr++)
                            {
                              int tChannel = crystal[iCry].slice[0].tChannel[jCorr];
                              // find delay as function of doi
                              Float_t delay;
                              if(jCorr == 0)
                              {
                                delay = 0; // no delay for direct detector
                              }
                              else
                              {
                                delay = crystal[iCry].delay_line[jCorr-1]->Eval(FloodZ);
                              }
                              // find Dt of this detector (t - tR) - delay;
                              Float_t Dt_j = (timeStamp[tChannel] - timeStamp[taggingCrystalTimingChannel]) - delay;

                              // now calc the weight_i
                              // keeping j fixed, run on i and sum the inv cov elements
                              Float_t weight_i = 0;
                              for(int iCorr = 0; iCorr < nLike ; iCorr++)
                              {
                                weight_i += crystal[iCry].inverse_covariance_element[iCorr][jCorr]->Eval(FloodZ);
                                sum_inv_cov += crystal[iCry].inverse_covariance_element[iCorr][jCorr]->Eval(FloodZ);
                              }
                              Dt_best += weight_i * Dt_j;
                            }
                            Dt_best = Dt_best / sum_inv_cov;
                            double likeCTR = Dt_best + centralcorrection ;
                            crystal[iCry].likeCTR->Fill(likeCTR);
                            crystal[iCry].vLike.push_back(likeCTR);
                          }
                        }
                      }
                    }
                    // crystal[iCry].allCTR->Fill(averageTimeStamp + zeroCorrection  - timeStamp[taggingCrystalTimingChannel]);
                  }
                  // hybrid correction
                  //
                  // if(hybridCorrection)
                  // {
                  //   //central time stamp
                  //   // std::cout << "----------------" << std::endl;
                  //   float weight = 0.0;
                  //   float meanTimeStamp = 0.0;
                  //   float sumWeight = 0.0;
                  //   // std::cout << "crystal[iCry].fwhmForPolishedCorrection[0] = " << crystal[iCry].fwhmForPolishedCorrection[0]<< std::endl;
                  //
                  //   weight = pow(crystal[iCry].fwhmForPolishedCorrection[0],-2);
                  //   float t_0 = timeStamp[crystal[iCry].tChannelsForPolishedCorrection[0]];
                  //
                  //   meanTimeStamp += weight * t_0;
                  //   sumWeight += weight;
                  //   // std::cout << weight << " " << t_0 << " " << meanTimeStamp << " " << sumWeight << std::endl;
                  //
                  //   // std::cout << t_0 << "\t";
                  //
                  //   for(unsigned int iPoli = 1; iPoli < crystal[iCry].tChannelsForPolishedCorrection.size(); iPoli++)
                  //   {
                  //     // std::cout << iPoli << std::endl;
                  //     // std::cout << "timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] = " << timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] << std::endl;
                  //     // std::cout << "crystal[iCry].meanForPolishedCorrection[iPoli] = " << crystal[iCry].meanForPolishedCorrection[iPoli] << std::endl;
                  //
                  //     float t_i = timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] - crystal[iCry].meanForPolishedCorrection[iPoli];
                  //     float weight_i = pow(crystal[iCry].fwhmForPolishedCorrection[iPoli],-2);
                  //     meanTimeStamp += weight_i * t_i;
                  //     sumWeight += weight_i;
                  //     // std::cout << weight_i << " " << t_i << " " << meanTimeStamp << " " << sumWeight << std::endl;
                  //     // std::cout << t_i << "\t";
                  //   }
                  //   // std::cout << std::endl;
                  //   meanTimeStamp = meanTimeStamp / sumWeight;
                  //   // std::cout << std::endl
                  //   //           << meanTimeStamp << "\t"
                  //   //           << std::endl;
                  //   double hybridCorrCTR = meanTimeStamp - timeStamp[taggingCrystalTimingChannel];
                  //   //correct by doi
                  //
                  //   hybridCorrCTR = hybridCorrCTR + centralcorrection; //FIXME temp
                  //
                  //   crystal[iCry].hybridCTR->Fill(hybridCorrCTR);
                  //   crystal[iCry].vhybrid.push_back(hybridCorrCTR);
                  //
                  //
                  // }
                }
              }

              if(crystal[iCry].polishedCorrection) //FIXME bugged now..
              {
                //central time stamp
                // std::cout << "----------------" << std::endl;

                float meanTimeStamp = 0.0;
                float sumWeight = 0.0;
                // std::cout << "crystal[iCry].fwhmForPolishedCorrection[0] = " << crystal[iCry].fwhmForPolishedCorrection[0]<< std::endl;

                bool noZeroes = true;
                if(timeStamp[taggingCrystalTimingChannel] == 0)
                {
                  noZeroes = false;
                }
                for(unsigned int iDet = 0; iDet < crystal[iCry].polished_correction.size(); iDet++)
                {
                  int timingChannel = crystal[iCry].polished_correction[iDet].timingChannel;
                  if(timeStamp[timingChannel] == 0)
                  {
                    noZeroes = false;
                  }
                }
                for(unsigned int iDet = 0; iDet < crystal[iCry].polished_correction.size(); iDet++)
                {
                  int timingChannel = crystal[iCry].polished_correction[iDet].timingChannel;
                  float delay = 0;
                  if(crystal[iCry].polished_correction[iDet].timingChannel == crystal[iCry].timingChannel)
                  {
                    delay = 0;
                  }
                  else
                  {
                    delay = crystal[iCry].polished_correction[iDet].mean;
                  }
                  float correctedElement = timeStamp[timingChannel] - timeStamp[taggingCrystalTimingChannel] - delay;
                  if(correctedElement <=  histoMin || correctedElement >= histoMax )
                  {
                    noZeroes = false;
                  }
                }

                if(noZeroes)
                {
                  // std::cout << "in " << crystal[iCry].polished_correction.size() <<  std::endl;
                  for(unsigned int iPoli = 0; iPoli < crystal[iCry].polished_correction.size(); iPoli++)
                  {
                    // std::cout << "ciao " << crystal[iCry].polished_correction.size() <<  std::endl;
                    // std::cout << iPoli <<" ";
                    float delay = 0;
                    float weight = 0.0;
                    if(crystal[iCry].polished_correction[iPoli].timingChannel == crystal[iCry].timingChannel)
                    {
                      delay = 0;
                    }
                    else
                    {
                      delay = crystal[iCry].polished_correction[iPoli].mean;
                    }
                    float rms = crystal[iCry].polished_correction[iPoli].rms;
                    weight = pow(rms,-2);
                    sumWeight += weight;
                    float correctedTimepstamp = timeStamp[crystal[iCry].polished_correction[iPoli].timingChannel] - timeStamp[taggingCrystalTimingChannel] - delay;
                    meanTimeStamp += weight * correctedTimepstamp;

                    // std::cout << rms << "\t"
                    //           << weight << "\t"
                    //           << sumWeight << "\t"
                    //           << correctedTimepstamp << "\t"
                    //           << meanTimeStamp << "\t"
                    //           << std::endl;
                    // std::cout << std::flush;
                  }
                  meanTimeStamp = meanTimeStamp/sumWeight;



                  double poliCorrCTR = meanTimeStamp;
                  crystal[iCry].poliCorrCTR->Fill(poliCorrCTR);
                  crystal[iCry].vPoli.push_back(poliCorrCTR);
                }






              }



              // end of temp commented
            }
          }
        }
        // delete Formula;
      }
    }

    counterAnalysis++;

    int perc = ((100*counterAnalysis)/neventAnalysis); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
  }
  std::cout << "Good events = " << goodEventsAnalysis << std::endl;
  std::sort(crystal.begin(), crystal.end(), compareByNumber);

  // TFile *outputCryFile = new TFile("testCry.root","RECREATE");
  //
  // outputCryFile->cd();
  // outputCryFile


  TH1F* noCorr = new TH1F("No Correction","No Correction",bins,minCTR,maxCTR);
  TH1F* centralCorr = new TH1F("Central Correction","Central Correction",bins,minCTR,maxCTR);
  TH1F* fullCorr = new TH1F("Full Correction","Full Correction",bins,minCTR,maxCTR);
  TH1F* likeCorr = new TH1F("Likelihood Correction","Likelihood Correction",bins,minCTR,maxCTR);
  TH1F* poliCorr = new TH1F("Polished Correction","Polished Correction",bins,minCTR,maxCTR);
  TH1F* hybridCorr = new TH1F("Hybrid Correction","Hybrid Correction",bins,minCTR,maxCTR);

  TH1F* unbinnednoCorr      = new TH1F("Unbinned No Correction","No Correction",bins,minCTR,maxCTR);
  TH1F* unbinnedcentralCorr = new TH1F("Unbinned Central Correction","Central Correction",bins,minCTR,maxCTR);
  TH1F* unbinnedfullCorr    = new TH1F("Unbinned Full Correction","Full Correction",bins,minCTR,maxCTR);
  TH1F* unbinnedpoliCorr    = new TH1F("Unbinned Polished Correction","Polished Correction",bins,minCTR,maxCTR);

  // std::vector<TH1F*> histograms;

  // do summary canvases for checking the fits

  int sqrtCrystals = ceil(sqrt( crystal.size() ) );

  TCanvas *cSumSimple  = new TCanvas("Summary Basic CTR","Summary Basic CTR",1200,1200);
  TCanvas *cSumCentral = new TCanvas("Summary Central CTR","Summary Central CTR",1200,1200);
  TCanvas *cSumAll     = new TCanvas("Summary Full CTR","Summary Full CTR",1200,1200);
  TCanvas *cPoliAll     = new TCanvas("Summary Polished CTR","Summary Polished CTR",1200,1200);
  cSumSimple ->Divide(sqrtCrystals,sqrtCrystals);
  cSumCentral->Divide(sqrtCrystals,sqrtCrystals);
  cSumAll->Divide(sqrtCrystals,sqrtCrystals);
  cPoliAll->Divide(sqrtCrystals,sqrtCrystals);

  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();

  std::cout << "Saving results..." << std::endl;

  for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  {

    if(likelihood)
    {
      if(crystal[iCry].accepted)
      {
        //create the 2d matrix of tgraphs
        int nLike = crystal[iCry].slice[0].tChannel.size();
        for(int iCorr = 0; iCorr < nLike ; iCorr++)
        {
          for(int jCorr = 0; jCorr < nLike ; jCorr++)
          {
            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->Write();
          }
        }
      }
    }

    // if(crystal[iCry].accepted)
    // {
    //   for(unsigned int iDet = 0 ;iDet < crystal[iCry].ctr_aligned.size(); iDet++)
    //   {
    //
    //     crystal[iCry].ctr_aligned[iDet].original_scatter->Write();
    //     crystal[iCry].ctr_aligned[iDet].aligned_scatter->Write();
    //
    //     if(crystal[iCry].ctr_aligned[iDet].delay_graph != NULL)
    //     {
    //       crystal[iCry].ctr_aligned[iDet].delay_graph->Write();
    //     }
    //
    //     if(crystal[iCry].ctr_aligned[iDet].delay_rms_graph != NULL)
    //     {
    //       crystal[iCry].ctr_aligned[iDet].delay_rms_graph->Write();
    //     }
    //
    //
    //
    //     for(unsigned int iSlice = 0; iSlice < crystal[iCry].ctr_aligned[iDet].slice.size();iSlice++)
    //     {
    //       crystal[iCry].ctr_aligned[iDet].slice[iSlice]->Write();
    //     }
    //
    //     crystal[iCry].ctr_aligned[iDet].ctr_aligned_graph->Write();
    //     crystal[iCry].ctr_aligned[iDet].ctr_aligned_rms_graph->Write();
    //
    //   }
    //
    // }

    std::stringstream sname;

    if(smooth)
    {
      if(crystal[iCry].simpleCTR)
      {
        crystal[iCry].simpleCTR ->Smooth(smooth);
      }
      if(crystal[iCry].centralCTR)
      {
        crystal[iCry].centralCTR ->Smooth(smooth);
      }
      if(crystal[iCry].allCTR)
      {
        crystal[iCry].allCTR ->Smooth(smooth);
      }
      if(crystal[iCry].poliCorrCTR)
      {
        crystal[iCry].poliCorrCTR ->Smooth(smooth);
      }
    }


    Float_t realBasicCTRfwhm,realBasicCTRfwtm ;
    Float_t realCentralCTRfwhm,realCentralCTRfwtm;
    Float_t realAllCTRfwhm,realAllCTRfwtm;
    Float_t poliCorrCTRfwhm,poliCorrCTRfwtm;
    Float_t reallikeCTRfwhm,reallikeCTRfwtm;
    Float_t realhybridCTRfwhm,realhybridCTRfwtm;
    double unbinnedSimpleCTR;
    double unbinnedCentralCTR;
    double unbinnedAllCTR;
    double unbinnedPoliCTR;
    double ret[2];
    double fitRes[3];
    // Int_t CTRentries;
    Float_t lightCentral;
    Float_t lightAll;

    // get data on entries and light collected
    // CTR entries

    crystal[iCry].basicCTRhisto->Write();

    // light central
    TF1 *gaussCentral = new TF1("gaussCentral","gaus");
    crystal[iCry].lightCentralHisto->Fit(gaussCentral,"Q");
    lightCentral = gaussCentral->GetParameter(1);
    crystal[iCry].lightCentralHisto->Write();

    // light central
    TF1 *gaussAll = new TF1("gaussAll","gaus");
    crystal[iCry].lightAllHisto->Fit(gaussAll,"Q");
    lightAll = gaussAll->GetParameter(1);
    crystal[iCry].lightAllHisto->Write();

    // std::cout << "BASIC CTRs --------------------" << std::endl;
    // std::cout << crystal[iCry]
    if(crystal[iCry].simpleCTR)
    {
      Int_t CTRentries = crystal[iCry].simpleCTR->GetEntries();
      crystal[iCry].simpleCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].simpleCTR->SetFillStyle(3001);
      crystal[iCry].simpleCTR->SetFillColor(kGreen);
      crystal[iCry].simpleCTR->SetLineColor(kGreen);
      crystal[iCry].simpleCTR->SetStats(0);
      crystal[iCry].simpleCTR_norm = (TH1F*) crystal[iCry].simpleCTR->Clone();
      if(func == 0)
      {

        extractCTR(crystal[iCry].simpleCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithEMG(crystal[iCry].simpleCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
          // extractWithGaussAndExp(crystal[iCry].simpleCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].simpleCTR,0.68,true,tagFwhm);
        }
      }

      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      std::cout << "No corr    - cry " << crystal[iCry].number << "\t"
                << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << "\t"
                << CTRentries << "\t"
                << lightCentral << "\t"
                << lightAll << "\t"
                << fitRes[0] << "\t"
                << fitRes[1] << "\t"
                << fitRes[2] << "\t"
                << std::endl;


      //
      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      textfile  << "No corr    - cry " << crystal[iCry].number << "\t"
                << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << "\t"
                << CTRentries << "\t"
                << lightCentral << "\t"
                << lightAll << "\t"
                << fitRes[0] << "\t"
                << fitRes[1] << "\t"
                << fitRes[2] << "\t"
                << std::endl;

      realBasicCTRfwhm = ret[0]*1e12;
      realBasicCTRfwtm = ret[1]*1e12;
      noCorr->Fill(ret[0]*1e12);

      crystal[iCry].simpleCTR->Write();
      cSumSimple->cd(iCry+1);
      crystal[iCry].simpleCTR->Draw();
      crystal[iCry].simpleCTR_norm->Scale(1.0/crystal[iCry].simpleCTR_norm->GetMaximum());

      // use unbinned method
      if(unbinned)
      {
        double mean,meanErr,min,max;
        double delta = FindSmallestInterval(mean,
                                            meanErr,
                                            min,
                                            max,
                                            crystal[iCry].vSimple,
                                            0.68,
                                            true);
        //now pass to fwhm
        double fwhm = 2.355 * (delta/2.0);
        unbinnedSimpleCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
        unbinnednoCorr->Fill(unbinnedSimpleCTR);

        std::cout << "Unbinned No corr    - cry " << crystal[iCry].number << "\t"
                  << unbinnedSimpleCTR << "\t"
                  << 0 << std::endl;

        textfile  << "Unbinned No corr    - cry " << crystal[iCry].number << "\t"
                  << unbinnedSimpleCTR << "\t"
                  << 0 << std::endl;
      }
    }

    if(crystal[iCry].centralCTR)
    {
      Int_t CTRentries = crystal[iCry].centralCTR->GetEntries();
      crystal[iCry].centralCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].centralCTR->SetFillStyle(3001);
      crystal[iCry].centralCTR->SetFillColor(kBlue);
      crystal[iCry].centralCTR->SetLineColor(kBlue);
      crystal[iCry].centralCTR->SetStats(0);
      crystal[iCry].centralCTR_norm = (TH1F*) crystal[iCry].centralCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithEMG(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
          // extractWithGaussAndExp(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].centralCTR,0.68,true,tagFwhm);
        }
      }


      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      std::cout << "Central    - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Central    - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      realCentralCTRfwhm = ret[0]*1e12;
      realCentralCTRfwtm = ret[1]*1e12;
      centralCorr->Fill(ret[0]*1e12);

      crystal[iCry].centralCTR->Write();
      cSumCentral->cd(iCry+1);
      crystal[iCry].centralCTR->Draw();
      // crystal[iCry].centralCTR->Scale(1.0/crystal[iCry].centralCTR->GetMaximum());
      crystal[iCry].centralCTR_norm->Scale(1.0/crystal[iCry].centralCTR_norm->GetMaximum());

      // use unbinned method
      if(unbinned)
      {
        double mean,meanErr,min,max;
        double delta = FindSmallestInterval(mean,
                                            meanErr,
                                            min,
                                            max,
                                            crystal[iCry].vCentral,
                                            0.68,
                                            true);
        //now pass to fwhm
        double fwhm = 2.355 * (delta/2.0);
        unbinnedCentralCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
        unbinnedcentralCorr->Fill(unbinnedCentralCTR);

        std::cout << "Unbinned Central    - cry " << crystal[iCry].number << "\t"
                  << unbinnedCentralCTR << "\t"
                  << 0 << std::endl;

        textfile  << "Unbinned Central    - cry " << crystal[iCry].number << "\t"
                  << unbinnedCentralCTR << "\t"
                  << 0 << std::endl;
      }
    }

    if(crystal[iCry].allCTR)
    {
      Int_t CTRentries = crystal[iCry].allCTR->GetEntries();
      crystal[iCry].allCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].allCTR->SetFillStyle(3001);
      crystal[iCry].allCTR->SetFillColor(kRed);
      crystal[iCry].allCTR->SetLineColor(kRed);
      crystal[iCry].allCTR->SetStats(0);

      crystal[iCry].allCTR_norm = (TH1F*) crystal[iCry].allCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithEMG(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
          // extractWithGaussAndExp(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].allCTR,0.68,true,tagFwhm);
        }
      }
      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      std::cout << "Full corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Full corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      realAllCTRfwhm = ret[0]*1e12;
      realAllCTRfwtm = ret[1]*1e12;
      fullCorr->Fill(ret[0]*1e12);
      crystal[iCry].allCTR->Write();
      cSumAll->cd(iCry+1);
      crystal[iCry].allCTR->Draw();
      crystal[iCry].allCTR_norm->Scale(1.0/crystal[iCry].allCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());

      // use unbinned method
      if(unbinned)
      {
        double mean,meanErr,min,max;
        double delta = FindSmallestInterval(mean,
                                            meanErr,
                                            min,
                                            max,
                                            crystal[iCry].vAll,
                                            0.68,
                                            true);
        //now pass to fwhm
        double fwhm = 2.355 * (delta/2.0);
        unbinnedAllCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
        unbinnedfullCorr->Fill(unbinnedAllCTR);

        std::cout << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
                  << unbinnedAllCTR << "\t"
                  << 0 << std::endl;

        textfile  << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
                  << unbinnedAllCTR << "\t"
                  << 0 << std::endl;
      }
    }





    if(crystal[iCry].likeCTR)
    {
      Int_t CTRentries = crystal[iCry].likeCTR->GetEntries();
      crystal[iCry].likeCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].likeCTR->SetFillStyle(3001);
      crystal[iCry].likeCTR->SetFillColor(kRed);
      crystal[iCry].likeCTR->SetLineColor(kRed);
      crystal[iCry].likeCTR->SetStats(0);

      crystal[iCry].likeCTR_norm = (TH1F*) crystal[iCry].likeCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].likeCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithEMG(crystal[iCry].likeCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
          // extractWithGaussAndExp(crystal[iCry].likeCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].likeCTR,0.68,true,tagFwhm);
        }
      }



      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      std::cout << "Like corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Full corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;


      reallikeCTRfwhm = ret[0]*1e12;
      reallikeCTRfwtm = ret[1]*1e12;
      likeCorr->Fill(ret[0]*1e12);
      crystal[iCry].likeCTR->Write();
      cSumAll->cd(iCry+1);
      crystal[iCry].likeCTR->Draw();
      crystal[iCry].likeCTR_norm->Scale(1.0/crystal[iCry].likeCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());

      // use unbinned method
      // if(unbinned)
      // {
      //   double mean,meanErr,min,max;
      //   double delta = FindSmallestInterval(mean,
      //                                       meanErr,
      //                                       min,
      //                                       max,
      //                                       crystal[iCry].vAll,
      //                                       0.68,
      //                                       true);
      //   //now pass to fwhm
      //   double fwhm = 2.355 * (delta/2.0);
      //   unbinnedAllCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
      //   unbinnedfullCorr->Fill(unbinnedAllCTR);
      //
      //   std::cout << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedAllCTR << "\t"
      //             << 0 << std::endl;
      //
      //   textfile  << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedAllCTR << "\t"
      //             << 0 << std::endl;
      // }
    }



    if(crystal[iCry].hybridCTR)
    {
      Int_t CTRentries = crystal[iCry].hybridCTR->GetEntries();
      crystal[iCry].hybridCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].hybridCTR->SetFillStyle(3001);
      crystal[iCry].hybridCTR->SetFillColor(kRed);
      crystal[iCry].hybridCTR->SetLineColor(kRed);
      crystal[iCry].hybridCTR->SetStats(0);

      crystal[iCry].hybridCTR_norm = (TH1F*) crystal[iCry].hybridCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].hybridCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithEMG(crystal[iCry].hybridCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
          // extractWithGaussAndExp(crystal[iCry].hybridCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].hybridCTR,0.68,true,tagFwhm);
        }
      }



      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      std::cout << "Hybr corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Full corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;


      realhybridCTRfwhm = ret[0]*1e12;
      realhybridCTRfwtm = ret[1]*1e12;
      hybridCorr->Fill(ret[0]*1e12);
      crystal[iCry].hybridCTR->Write();
      cSumAll->cd(iCry+1);
      crystal[iCry].hybridCTR->Draw();
      crystal[iCry].hybridCTR_norm->Scale(1.0/crystal[iCry].hybridCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());

      // use unbinned method
      // if(unbinned)
      // {
      //   double mean,meanErr,min,max;
      //   double delta = FindSmallestInterval(mean,
      //                                       meanErr,
      //                                       min,
      //                                       max,
      //                                       crystal[iCry].vAll,
      //                                       0.68,
      //                                       true);
      //   //now pass to fwhm
      //   double fwhm = 2.355 * (delta/2.0);
      //   unbinnedAllCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
      //   unbinnedfullCorr->Fill(unbinnedAllCTR);
      //
      //   std::cout << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedAllCTR << "\t"
      //             << 0 << std::endl;
      //
      //   textfile  << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedAllCTR << "\t"
      //             << 0 << std::endl;
      // }
    }








    if(crystal[iCry].poliCorrCTR)
    {
      Int_t CTRentries = crystal[iCry].poliCorrCTR->GetEntries();
      crystal[iCry].poliCorrCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].poliCorrCTR->SetFillStyle(3001);
      crystal[iCry].poliCorrCTR->SetFillColor(kBlack);
      crystal[iCry].poliCorrCTR->SetLineColor(kBlack);
      crystal[iCry].poliCorrCTR->SetStats(0);

      crystal[iCry].poliCorrCTR_norm = (TH1F*) crystal[iCry].poliCorrCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithEMG(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
          // extractWithGaussAndExp(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].poliCorrCTR,0.68,true,tagFwhm);
        }
      }


      // std::cout << "Polished corr. - cry " << crystal[iCry].number << "\t"
      //           << ret[0]*1e12 << "\t"
      //           << ret[1]*1e12 << std::endl;
      // //
      // std::cout << "Polished FIT   - cry " << crystal[iCry].number << "\t"
      //           << fitRes[0] << "\t"
      //           << fitRes[1] << "\t"
      //           << fitRes[2] << "\t"
      //           << std::endl;
      //
      // textfile  << "Polished corr. - cry " << crystal[iCry].number << "\t"
      //           << ret[0]*1e12 << "\t"
      //           << ret[1]*1e12 << std::endl;
      // //
      // textfile << "Polished FIT   - cry " << crystal[iCry].number << "\t"
      //           << fitRes[0] << "\t"
      //           << fitRes[1] << "\t"
      //           << fitRes[2] << "\t"
      //           << std::endl;


      //
      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      std::cout << "Polished corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Polished corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      poliCorrCTRfwhm = ret[0]*1e12;
      poliCorrCTRfwtm = ret[1]*1e12;
      poliCorr->Fill(ret[0]*1e12);
      crystal[iCry].poliCorrCTR->Write();
      cPoliAll->cd(iCry+1);
      crystal[iCry].poliCorrCTR->Draw();
      crystal[iCry].poliCorrCTR_norm->Scale(1.0/crystal[iCry].poliCorrCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());


      // use unbinned method
      if(unbinned)
      {
        double mean,meanErr,min,max;
        double delta = FindSmallestInterval(mean,
                                            meanErr,
                                            min,
                                            max,
                                            crystal[iCry].vPoli,
                                            0.68,
                                            true);
        //now pass to fwhm
        double fwhm = 2.355 * (delta/2.0);
        unbinnedPoliCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
        unbinnedpoliCorr->Fill(unbinnedPoliCTR);

        std::cout << "Unbinned Polished corr. - cry " << crystal[iCry].number << "\t"
                  << unbinnedPoliCTR << "\t"
                  << 0 << std::endl;

        textfile  << "Unbinned Polished corr. - cry " << crystal[iCry].number << "\t"
                  << unbinnedPoliCTR << "\t"
                  << 0 << std::endl;
      }
    }



    sname.str("");

    sname << "Summary - Crystal " << crystal[iCry].number;
    TCanvas* c_summary = new TCanvas(sname.str().c_str(),sname.str().c_str(),1200,800);
    c_summary->cd();
    THStack *hs = new THStack("hs","");
    hs->Add(crystal[iCry].simpleCTR_norm);
    hs->Add(crystal[iCry].centralCTR_norm);
    hs->Add(crystal[iCry].allCTR_norm);
    hs->Add(crystal[iCry].likeCTR_norm);
    hs->Add(crystal[iCry].poliCorrCTR_norm);
    hs->Add(crystal[iCry].hybridCTR_norm);


    // std::cout << "Crystal " << crystal[iCry].number << std::endl;
    // std::cout << "CTR FWHM [ps] " << std::endl;
    hs->Draw("hist nostack");
    sname.str("");
    sname << "CTR - Crystal " << crystal[iCry].number << " - width in FWHM";
    hs->SetTitle(sname.str().c_str());
    hs->GetXaxis()->SetTitle("Time [s]");
    hs->GetXaxis()->SetTitleOffset(1);
    hs->GetXaxis()->SetTitleSize(0.045);
    hs->GetXaxis()->SetLabelSize(0.045);
    hs->GetYaxis()->SetLabelSize(0.045);
    TLegend *legend = new TLegend(0.54,0.62,0.89,0.89,"");
    legend->SetFillStyle(0);
    if(crystal[iCry].simpleCTR)
    {
      sname.str("");
      sname << "No correction        = " << realBasicCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].simpleCTR,sname.str().c_str(),"f");
      // std::cout << "No correction       = "<< realBasicCTRfwhm   << " ps" << std::endl;
    }
    if(crystal[iCry].centralCTR)
    {
      sname.str("");
      sname << "Central correction = " << realCentralCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].centralCTR,sname.str().c_str(),"f");
      // std::cout << "Central correction  = "<< realCentralCTRfwhm << " ps" << std::endl;
    }
    if(crystal[iCry].allCTR)
    {
      sname.str("");
      sname << "Full correction       = " << realAllCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].allCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }
    if(crystal[iCry].likeCTR)
    {
      sname.str("");
      sname << "Likelihood correction = " << reallikeCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].likeCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }
    if(crystal[iCry].poliCorrCTR)
    {
      sname.str("");
      sname << "Polished correction       = " << poliCorrCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].poliCorrCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }
    if(crystal[iCry].hybridCTR)
    {
      sname.str("");
      sname << "Polished correction       = " << realhybridCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].hybridCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }


    sname.str("");
    legend->Draw();
    gStyle->SetOptTitle(0);
    TPaveLabel *title = new TPaveLabel(.11,.95,.35,.99,"new title","brndc");
    title->Draw();
    // std::cout << std::endl;

    c_summary->Write();


    TH1F* cloneBasic;
    TH1F* cloneCentral;
    TH1F* cloneAll;
    TH1F* cloneLike;
    TH1F* clonePoli;
    TH1F* cloneHybrid;
    THStack *cloneHs = (THStack*) hs->Clone();
    TLegend *legend1 = new TLegend(0.15,0.69,0.49,0.89,"");
    legend1->SetFillStyle(0);
    sname.str("");
    sname << "Multi - Crystal " << crystal[iCry].number;
    TCanvas* c_multi = new TCanvas(sname.str().c_str(),sname.str().c_str(),1800,1400);
    c_multi->Divide(2,2);

    if(crystal[iCry].simpleCTR_norm)
    {
      cloneBasic   = (TH1F*) crystal[iCry].simpleCTR->Clone();
      c_multi->cd(1);
      sname.str("");
      sname << "No correction        = " << realBasicCTRfwhm << "ps";
      legend1->AddEntry(cloneBasic,sname.str().c_str(),"f");
      cloneBasic->Draw();
      legend1->Draw();
    }
    if(crystal[iCry].centralCTR_norm)
    {
      cloneCentral = (TH1F*) crystal[iCry].centralCTR->Clone();
      c_multi->cd(2);
      TLegend *legend2 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend2->SetFillStyle(0);
      sname.str("");
      sname << "Central correction   = " << realCentralCTRfwhm << "ps";
      legend2->AddEntry(cloneCentral,sname.str().c_str(),"f");
      cloneCentral->Draw();
      legend2->Draw();
    }
    if(crystal[iCry].allCTR_norm)
    {
      cloneAll     = (TH1F*) crystal[iCry].allCTR->Clone();
      c_multi->cd(3);
      TLegend *legend3 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend3->SetFillStyle(0);
      sname.str("");
      sname << "Full correction      = " << realAllCTRfwhm << "ps";
      legend3->AddEntry(cloneAll,sname.str().c_str(),"f");
      cloneAll->Draw();
      legend3->Draw();
    }
    if(crystal[iCry].likeCTR_norm)
    {
      cloneLike     = (TH1F*) crystal[iCry].likeCTR->Clone();
      c_multi->cd(3);
      TLegend *legend4 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend4->SetFillStyle(0);
      sname.str("");
      sname << "Likelihood correction = " << reallikeCTRfwhm << "ps";
      legend4->AddEntry(cloneLike,sname.str().c_str(),"f");
      cloneLike->Draw();
      legend4->Draw();
    }
    if(crystal[iCry].hybridCTR_norm)
    {
      cloneHybrid     = (TH1F*) crystal[iCry].hybridCTR->Clone();
      c_multi->cd(3);
      TLegend *legend5 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend5->SetFillStyle(0);
      sname.str("");
      sname << "Hybrid correction = " << realhybridCTRfwhm << "ps";
      legend5->AddEntry(cloneHybrid,sname.str().c_str(),"f");
      cloneHybrid->Draw();
      legend5->Draw();
    }


    // if(crystal[iCry].poliCorrCTR_norm)
    // {
    //   clonePoli     = (TH1F*) crystal[iCry].poliCorrCTR->Clone();
    //   c_multi->cd(3);
    //   TLegend *legend3 = new TLegend(0.15,0.69,0.49,0.89,"");
    //   legend3->SetFillStyle(0);
    //   sname.str("");
    //   sname << "Polished correction      = " << poliCorrCTRfwhm << "ps";
    //   legend3->AddEntry(cloneAll,sname.str().c_str(),"f");
    //   cloneAll->Draw();
    //   legend3->Draw();
    // }


    c_multi->cd(4);
    c_multi->cd(4)->SetGrid();
    cloneHs->Draw("hist nostack");
    c_multi->Write();

  }
  noCorr->Write();
  centralCorr->Write();
  fullCorr->Write();
  poliCorr->Write();
  likeCorr->Write();
  hybridCorr->Write();

  unbinnednoCorr->Write();
  unbinnedcentralCorr->Write();
  unbinnedfullCorr->Write();
  unbinnedpoliCorr->Write();

  cSumSimple ->Write();
  cSumCentral->Write();
  cSumAll->Write();
  cPoliAll->Write();


  TNamed CommandNameD("Command",streamCommand.str().c_str());
  CommandNameD.Write();
  // treeFile->Close();

  calibrationFile->Close();
  outputFile->Close();
  textfile.close();
  std::cout << std::endl;
  std::cout << "Histograms saved in   " << outputFileName << std::endl;
  std::cout << "Text summary saved in " << textFileName << std::endl;
  return 0;
}
