//---------------------------------------------------//
//                                                   //
// EXTRACT WITH EXPONENTIALLY MODIFIED GAUSSIAN      //
//                                                   //
//---------------------------------------------------//
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
//---------------------------------------------------//




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
