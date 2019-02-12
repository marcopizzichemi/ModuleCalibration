
//---------------------------------------------------//
//                                                   //
// EXTRACT WITH EXPONENTIALLY MODIFIED GAUSSIAN      //
//                                                   //
//---------------------------------------------------//


// functions
// paper-like EMG
// TString gexp_string = "[0]/sqrt(2)*exp([2]^2/2/[3]^2-(x-[1])/[3])*(1-TMath::Erf(([1]-x+[2]^2/[3])/(sqrt(2*[2]^2))))";
// TString gexp_inv_string = "[0]/sqrt(2)*exp([2]^2/2/[3]^2-((-x)-(-[1]))/[3])*(1-TMath::Erf(((-[1])-(-x)+[2]^2/[3])/(sqrt(2*[2]^2))))";

// standard EMG
// par[0] = A
// par[1] = mu
// par[2] = sigma
// par[3] = tau
Double_t emg_function(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];
   Double_t f;
   Double_t norm = par[0] /(2.0 * par[3]);
   Double_t expo = TMath::Exp( (1.0/(2.0*par[3])) * ( 2.0*par[1] + ( TMath::Power(par[2],2)/par[3] ) - 2.0*xx ) );
   Double_t erf_c = TMath::Erfc( ( par[1] + ( TMath::Power(par[2],2)/par[3] ) - xx )/(TMath::Sqrt(2)  * par[2] )  ) ;
   f = norm * expo * erf_c;
   return f;
}

// standard EMG inverted. x -> -x, then mu will be -mu
Double_t emginv_function(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];
   Double_t m = par[1] + par[3];
   Double_t f;
   Double_t norm = par[0] /(2.0 * par[3]);
   Double_t expo = TMath::Exp( (1.0/(2.0*par[3])) * ( 2.0*par[1] + ( TMath::Power(par[2],2)/par[3] ) + 2.0*(xx) ) );
   Double_t erf_c = TMath::Erfc( ( par[1] + ( TMath::Power(par[2],2)/par[3] ) + (xx ) )/(TMath::Sqrt(2)  * par[2] )  ) ;
   f = norm * expo * erf_c;
   return f;
}

// forward declarations
template <class T>
void fitWithEMG(T* histo,float* res); //just fitting, ignoring fitres and tagfwhm. this is used in ModuleCalibration to fit the slices and find mean. relevant values are in *res
template <class T>
void extractCTRWithEMG_withRef(T* histo,int divs, float tagFwhm,float* res,float* fitRes);//fit and directly calculate the fwhm on the function
template <class T>
void fitWithEMG_core(T* histo,float &fitMin,float &fitMax,float &f1min,float &f1max,int &fType,double *fParams); // actual EMG fitting routine

// routines
template <class T>
void fitWithEMG(T* histo,float* res)
{
  //just fitting, ignoring fitres and tagfwhm
  // this is used in ModuleCalibration to fit the slices and find mean. relevant values are in *res
  float fitMin = 0;
  float fitMax = 0;
  float f1min = 0;
  float f1max = 0;
  int fType;
  double fParams[4];
  fitWithEMG_core(histo,fitMin,fitMax,f1min,f1max,fType,fParams);

  if(fType == 0)
  {
    res[0] = 0;
    res[1] = 0;
    res[2] = 0;
    res[3] = 0;
  }
  else
  {
    TF1 *fitFunction;
    if(fType == 1)
    {
      fitFunction = new TF1("gexp",emg_function,f1min,f1max,4);
    }
    else
    {
      fitFunction = new TF1("gexp_inv",emginv_function,f1min,f1max,4);
    }
    fitFunction->SetParameter(0,fParams[0]);
    fitFunction->SetParameter(1,fParams[1]);
    fitFunction->SetParameter(2,fParams[2]);
    fitFunction->SetParameter(3,fParams[3]);

    histo->Fit(fitFunction,"Q","",fitMin,fitMax); // re-fit just to store only the good one
    // write variables or it's gonna be a mess
    float mu = fitFunction->GetParameter(1);
    float e_mu = fitFunction->GetParError(1);
    float s = fitFunction->GetParameter(2);
    float e_s = fitFunction->GetParError(2);
    float tau = fitFunction->GetParameter(3);
    float e_tau = fitFunction->GetParError(3);

    float mean;
    if(fType == 1)
    {
      mean = mu + tau;
    }
    else
    {
      mean = -(mu + tau);
    }
    float sigma = TMath::Sqrt(TMath::Power(s,2) + TMath::Power(tau,2));
    float meanErr = TMath::Sqrt(TMath::Power(e_mu,2) + TMath::Power(e_tau,2));
    float sigmaErr = TMath::Sqrt( ( (TMath::Power(s,2))/(TMath::Power(s,2) + TMath::Power(tau,2))  )*TMath::Power(e_s,2) + ((TMath::Power(tau,2))/(TMath::Power(s,2) + TMath::Power(tau,2)))*TMath::Power(e_tau,2));

    res[0] = mean;
    res[1] = sigma;
    res[2] = meanErr;
    res[3] = sigmaErr;
  }
}

//fit and directly calculate the fwhm on the function
template <class T>
void extractCTRWithEMG_withRef(T* histo,int divs, float tagFwhm,float* res,float* fitRes)
{
  float fitMin = 0;
  float fitMax = 0;
  float f1min = 0;
  float f1max = 0;
  double fParams[4];
  int fType;
  fitWithEMG_core(histo,fitMin,fitMax,f1min,f1max,fType,fParams);
  fitRes[0] = 0;
  fitRes[1] = 0;
  fitRes[2] = 0;

  TF1 *fitFunction;
  if(fType == 0)
  {
    fitFunction = NULL;
    res[0] = 0;
    res[1] = 0;
    res[2] = 0;
    res[3] = 0;
  }
  else
  {
    if(fType == 1)
    {
      fitFunction = new TF1("gexp",emg_function,f1min,f1max,4);
    }
    else
    {
      fitFunction = new TF1("gexp_inv",emginv_function,f1min,f1max,4);
    }
    fitFunction->SetParameter(0,fParams[0]);
    fitFunction->SetParameter(1,fParams[1]);
    fitFunction->SetParameter(2,fParams[2]);
    fitFunction->SetParameter(3,fParams[3]);

    histo->Fit(fitFunction,"Q","",fitMin,fitMax); // re-fit just to store only the good one
    float min,max,min10,max10;
    float step = (f1max-f1min)/divs;
    float funcMax = fitFunction->GetMaximum(fitMin,fitMax);
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
    res[2] = 0;
    res[3] = 0;

    fitRes[0] = fitFunction->GetChisquare();
    fitRes[1] = fitFunction->GetNDF();
    fitRes[2] = fitFunction->GetProb();
  }
}
//---------------------------------------------------//


template <class T>
void fitWithEMG_core(T* histo,float &fitMin,float &fitMax,float &f1min,float &f1max,int &fType,double *fParams) // actual EMG fitting routine
{
  // actual EMG fitting routine
  TCanvas *cTemp  = new TCanvas("temp","temp");
  f1min = histo->GetXaxis()->GetXmin();
  f1max = histo->GetXaxis()->GetXmax();
  fitMin = f1min;
  fitMax = f1max;
  //define the standard gauss plus exp and the same but with inverted sign for x
  // TF1* gexp      = new TF1("gexp",gexp_string,f1min,f1max); //paper-like
  TF1* gexp      = new TF1("gexp",emg_function,f1min,f1max,4);
  gexp->SetLineColor(kGreen);
  gexp->SetParName(0,"N");
  gexp->SetParName(1,"mu");
  gexp->SetParName(2,"Sigma");
  gexp->SetParName(3,"tau");
  gexp->SetParameters(histo->GetEntries()*histo->GetXaxis()->GetBinWidth(1),
                      histo->GetMean(),
                      histo->GetRMS(),
                      histo->GetRMS());
  TFitResultPtr r_gexp = histo->Fit(gexp,"QS","",fitMin,fitMax);

  // inverted emg
  TF1* gexp_inv  = new TF1("gexp_inv",emginv_function,f1min,f1max,4);
  gexp_inv->SetLineColor(kBlue);
  gexp_inv->SetParName(0,"N");
  gexp_inv->SetParName(1,"mu");
  gexp_inv->SetParName(2,"Sigma");
  gexp_inv->SetParName(3,"tau");
  gexp_inv->SetParameters(histo->GetEntries()*histo->GetXaxis()->GetBinWidth(1),
                          -histo->GetMean(),
                          histo->GetRMS(),
                          histo->GetRMS());
  TFitResultPtr r_gexp_inv = histo->Fit(gexp_inv,"QS","",fitMin,fitMax);
  //try both, see what fits better
  Int_t fitStatusGexp = r_gexp;
  Int_t fitStatusGexp_inv = r_gexp_inv;
  float chi2gexp = 0;
  float chi2gexp_inv = 0;

  if(fitStatusGexp == 0) // if gexp worked
  {
    chi2gexp = r_gexp->Chi2();
  }
  if(fitStatusGexp_inv == 0)// if gexp_inv worked
  {
    chi2gexp_inv   = r_gexp_inv->Chi2();
  }

  if(fitStatusGexp == 0 && fitStatusGexp_inv != 0) //if gexp worked and gexp_inv didn't worked
  {
    fType = 1;
    fParams[0] = gexp->GetParameter(0);
    fParams[1] = gexp->GetParameter(1);
    fParams[2] = gexp->GetParameter(2);
    fParams[3] = gexp->GetParameter(3);
  }
  else
  {
    if(fitStatusGexp != 0 && fitStatusGexp_inv == 0) //if gexp didn't worked and gexp_inv worked
    {
      fType = 2;
      fParams[0] = gexp_inv->GetParameter(0);
      fParams[1] = gexp_inv->GetParameter(1);
      fParams[2] = gexp_inv->GetParameter(2);
      fParams[3] = gexp_inv->GetParameter(3);
    }
    else // both worked or nothing worked
    {
      if(fitStatusGexp == 0 && fitStatusGexp_inv == 0)
      {
        if(chi2gexp > chi2gexp_inv) // gexp_inv better than gexp
        {
          fType = 2;
          fParams[0] = gexp_inv->GetParameter(0);
          fParams[1] = gexp_inv->GetParameter(1);
          fParams[2] = gexp_inv->GetParameter(2);
          fParams[3] = gexp_inv->GetParameter(3);
        }
        else // gexp better than gexp_inv
        {
          fType = 1;
          fParams[0] = gexp->GetParameter(0);
          fParams[1] = gexp->GetParameter(1);
          fParams[2] = gexp->GetParameter(2);
          fParams[3] = gexp->GetParameter(3);
          // fitFunction = gexp;
        }
      }
      else  // nothing worked
      {
        // delete all func
        fType = 0;
        delete gexp_inv;
        delete gexp;
      }
    }
  }

  delete cTemp;
}



void extractCTR(TH1F* histo,float fitPercMin,float fitPercMax, int divs, float tagFwhm, float* res, float* fitRes)
{

  // preliminary gauss fit
  TCanvas *cTemp  = new TCanvas("temp","temp");
  TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  // resctrict the fitting range of gauss function

  gaussDummy->SetLineColor(kRed);
  float fitGaussMin = histo->GetMean()-2.0*histo->GetRMS();
  float fitGaussMax = histo->GetMean()+2.0*histo->GetRMS();
  float f1min = histo->GetXaxis()->GetXmin();
  float f1max = histo->GetXaxis()->GetXmax();
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

  float fitMin;
  float fitMax;
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

  float chi2gexp;
  float chi2cb;

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

      float min,max,min10,max10;
      // int divs = 3000;
      float step = (f1max-f1min)/divs;
      float funcMax = f1->GetMaximum(fitMin,fitMax);
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
float FindSmallestInterval(float& mean,
                            float& meanErr,
                            float& min,
                            float& max,
                            std::vector<float>& vals,
                            const float& fraction,
                            const bool& verbosity)
{
   if( verbosity )
     std::cout << ">>>>>> FindSmallestInterval" << std::endl;


   std::sort(vals.begin(),vals.end());

   unsigned int nPoints = vals.size();
   unsigned int maxPoints = (unsigned int)(fraction * nPoints);

   unsigned int minPoint = 0;
   unsigned int maxPoint = 0;
   float delta = 999999.;
   for(unsigned int point = 0; point < nPoints-maxPoints; ++point)
   {
     float tmpMin = vals.at(point);
     float tmpMax = vals.at(point+maxPoints-1);
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
void FindSmallestInterval(float* retValues, TH1F* histo, const float&
fraction, const bool& verbosity, float tagFwhm)
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
   float fwhm = 2.355* ((max-min) / 2.0);

   retValues[0] = sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
   retValues[1] = 0;
}
