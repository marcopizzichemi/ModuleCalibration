#include <iostream>

#include "Mppc.h"
#include "TSpectrum2.h"
#include "TF2.h"
#include "TMath.h"
#include "TEllipse.h"
Mppc::Mppc()
{
  //default contructor
  name = "Default MPPC";
  digitizerChannel = -1;
  canvasPosition = -1;
}

Mppc::Mppc(const Mppc &obj) 
{
  std::cout << "Copy ctor"  << std::endl;  
}

Mppc::~Mppc()
{
  std::cout << "destructing mppc " << name << std::endl;
}

void Mppc::SetModule(Module *amodule)
{
  parentModule = new Element;
  *parentModule = *((Element*)amodule);
}

void Mppc::SetCrystal(Crystal *pCrystal)
{
  Element *aElement = new Element;
  *aElement = *((Element*)pCrystal);
  vCrystal.push_back(aElement);
}

Crystal* Mppc::GetCrystal(int pi, int pj)
{
  return (Crystal*)vCrystal[pi*jChildren + pj];
}



int Mppc::Find2Dpeaks(int nofcrystals,TH2F* histogram2d)
{
//   int nofcrystals = 2;//FIXME for the moment hardcoded to find the line 
  TSpectrum2 *peak2d = new TSpectrum2(nofcrystals); //FIXME for the moment hardcoded to find the line 
  int nfound = peak2d->Search(histogram2d,2,"col",0.1);
  
//   std::cout << nfound << std::endl;
  
  //store what has been found in the params struct
  float *xpeaks = peak2d->GetPositionX();
  float *ypeaks = peak2d->GetPositionY();
  double *xsigma = new double[nofcrystals];
  double *ysigma = new double[nofcrystals];
  
  
  double *tempXsigma = new double[nofcrystals];
  double *tempYsigma = new double[nofcrystals];
  
//   fit2DmeanX = new double[nofcrystals];
//   fit2DmeanY = new double[nofcrystals];
//   fit2DsigmaX = new double[nofcrystals];
//   fit2DsigmaY = new double[nofcrystals];
  
  for(int i = 0 ; i < nofcrystals ; i++)
  {
    tempXsigma[i] = 0.2;
    tempYsigma[i] = 0.2;
    xsigma[i] = tempXsigma[i];
    ysigma[i] = tempYsigma[i];
  }
  
  
  //2d fit
  for ( int j = 0 ; j < nfound ; j++)
  {
    
    TF2 *f2 = new TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",xpeaks[j]-2.0*xsigma[j],xpeaks[j]+2.0*xsigma[j],ypeaks[j]-2.0*ysigma[j],ypeaks[j]+2.0*ysigma[j]);
    f2->SetParameters(14,xpeaks[j],xsigma[j],ypeaks[j],ysigma[j]);
    FloodMap2D.Fit("f2","QNR");
//     fit2DmeanX[j] = xpeaks[j];
//     fit2DmeanY[j] = ypeaks[j];
//     //we insert the 2.5 sigma limit here, it makes things easier
//     fit2DsigmaX[j] = 2.5*TMath::Abs(f2->GetParameter(2));
//     fit2DsigmaY[j] = 2.5*TMath::Abs(f2->GetParameter(4));
    fit2DmeanX.push_back(xpeaks[j]);
    fit2DmeanY.push_back(ypeaks[j]);
    //we insert the 2.5 sigma limit here, it makes things easier
    fit2DsigmaX.push_back(2.5*TMath::Abs(f2->GetParameter(2)));
    fit2DsigmaY.push_back(2.5*TMath::Abs(f2->GetParameter(4)));
    
  }
  
  
  for ( int j = 0 ; j < nfound ; j++)
  {
    float sx,sy;
    sx = fit2DsigmaX[j];
    sy = fit2DsigmaY[j];
    fit2DsigmaX[j] = (sx+sy)/2.0;
    fit2DsigmaY[j] = (sx+sy)/2.0;
  }
  
 
  
  
  
  
  
  //   now for each peak, check if any other peak is closer than the sum of the relative circles
  for ( int j = 0 ; j < nfound ; j++)// run on all peaks
  {
    for ( int jOther = 0 ; jOther < nfound ; jOther++)//run on all peaks again, but...
    {
      if (jOther != j) //...do it only if it's not the same peak
      {
	float distance = TMath::Sqrt( TMath::Power(fit2DmeanX[j]-fit2DmeanX[jOther],2) + TMath::Power(fit2DmeanY[j]-fit2DmeanY[jOther],2) );
	float sumOfRadii = (fit2DsigmaX[j]+fit2DsigmaX[jOther]);
	if ( distance < sumOfRadii )
	{
	  //std::cout << "WARNING: Peaks of Module " << k << " Channel " << i << " are overlapping!" << std::endl;
	  fit2DsigmaX[j] = distance * ( fit2DsigmaX[j] / sumOfRadii );
	  fit2DsigmaX[jOther] = distance * ( fit2DsigmaX[jOther] / sumOfRadii );
	  fit2DsigmaY[j] = fit2DsigmaX[j];
	  fit2DsigmaY[jOther] = fit2DsigmaX[jOther];
	}
      }
    }
  }
  
  // sort them so that first is lower y
  double swapFit2DmeanX,swapFit2DmeanY,swapFit2DsigmaX,swapFit2DsigmaY;
  if(fit2DmeanY[0] > fit2DmeanY[1])
  {
    swapFit2DmeanX  = fit2DmeanX[1];
    swapFit2DmeanY  = fit2DmeanY[1];
    swapFit2DsigmaX = fit2DsigmaX[1];
    swapFit2DsigmaY = fit2DsigmaY[1];
    
    fit2DmeanX[1] = fit2DmeanX[0];
    fit2DmeanY[1] = fit2DmeanY[0];
    fit2DsigmaX[1] = fit2DsigmaX[0];
    fit2DsigmaY[1] = fit2DsigmaY[0];
    
    fit2DmeanX[0] = swapFit2DmeanX  ;
    fit2DmeanY[0] = swapFit2DmeanY  ;
    fit2DsigmaX[0] =swapFit2DsigmaX ;
    fit2DsigmaY[0] =swapFit2DsigmaY ;
  }
  
//   for ( int j = 0 ; j < nfound ; j++)// run on all peaks
//   {
//     std::cout  <<  fit2DmeanX[j] << " " << fit2DmeanY[j]<< " " << fit2DsigmaX[j] << " " << fit2DsigmaY[j] << std::endl;
//   }
  
  //set the ellipses from here
//   for ( int j = 0 ; j < nfound ; j++) // runs max to j = 1
//   {
//     
//      // we have only two peaks. assign them to two crystals in a line in this mppc (which ones? doesn't matter for now.... FIXME)
//     this->GetCrystal(j,0)->SetCrystalOn(true);
//     this->GetCrystal(j,0)->SetCrystalData(fit2DmeanX[j],fit2DmeanY[j],fit2DsigmaX[j],fit2DsigmaY[j],0);
//     TEllipse *ellipse = new TEllipse(u,v,wu,wv,0,360,t);
//     this->GetCrystal(j,0)->SetGraphicalCut(*ellipse);
//     
//     //     ellipseCenterX[j] = fit2DmeanX[j];
//     //     ellipseCenterY[j] = fit2DmeanY[j];
//     //     ellipseWidthX[j] = fit2DsigmaX[j];
//     //     ellipseWidthY[j] = fit2DsigmaY[j];
//     //     TEllipse *ellipses;
//     //     ellipses = new TEllipse(ellipseCenterX[j],ellipseCenterY[j],ellipseWidthX[j],ellipseWidthY[j]);
//     //     //ellipses->SetFillColor(42);
//     //     ellipses->SetFillStyle(4001);
//     //     ellipses->SetLineColor(kRed);
//     //     ellipses->SetLineWidth(2);
//     //     ellipses->Draw();
//     std::cout  <<  fit2DmeanX[j] << " " << fit2DmeanY[j]<< " " << fit2DsigmaX[j] << " " << fit2DsigmaY[j] << std::endl;
//   }
  
  
  
  
  return nfound;
}






void Mppc::PrintSpecific()
{
  std::cout << "Label \t\t = " << label << std::endl;
  std::cout << "Digitizer Ch. \t = " << digitizerChannel << std::endl;
  std::cout << "Staus for Doi    \t = " << IsOnForDoi << std::endl;
}