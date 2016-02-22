#include <iostream>
#include "Mppc.h"
#include "TSpectrum2.h"
#include "TF2.h"
#include "TMath.h"
#include "TEllipse.h"
#include "TROOT.h"


Double_t g2(Double_t *x, Double_t *par) { //A Gaussian function
  Double_t r1 = Double_t((x[0]-par[1])/par[2]);
  Double_t r2 = Double_t((x[1]-par[3])/par[4]);
  return par[0]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}

//generic 2d gaussian
// par[] = {amplitude,mean_x,mean_y,sigma_x,sigma_y,theta}
// par[0] = amplitude
// par[1] = mean x 
// par[2] = mean y
// par[3] = sigma x
// par[4] = sigma y
// par[5] = theta

Double_t gauss2d(Double_t *x, Double_t *par) 
{
  Double_t xrot = x[0]*cos(par[5]) - x[1]*sin(par[5]);
  Double_t yrot = x[0]*sin(par[5]) + x[1]*cos(par[5]);
  Double_t x0rot = par[1]*cos(par[5]) - par[2]*sin(par[5]);
  Double_t y0rot = par[1]*sin(par[5]) + par[2]*cos(par[5]);
  
  
  Double_t result = par[0]* TMath::Exp( -( (xrot - x0rot )*(xrot - x0rot )/(2.0*par[3]*par[3]) + (yrot - y0rot )*(yrot - y0rot )/(2.0*par[4]*par[4]) ) );
  return result;
}


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
  /**Function to look for the 2D peaks
   * 
   */
  TSpectrum2 *peak2d = new TSpectrum2(nofcrystals); //FIXME for the moment hardcoded to find the line 
  int nfound = peak2d->Search(histogram2d,2,"col",0.1);
  
  //   std::cout << nfound << std::endl;
  
  //store what has been found in some array
  float *xpeaks = peak2d->GetPositionX();
  float *ypeaks = peak2d->GetPositionY();
  
  
  //   float *zpeaks = peak2d->GetPositionZ();
  
  double *xsigma = new double[nofcrystals];
  double *ysigma = new double[nofcrystals];
  //   double *tempXsigma = new double[nofcrystals];
  //   double *tempYsigma = new double[nofcrystals];
  
  
  // for each point, find the minimum fitting area
  for ( int j = 0 ; j < nfound ; j++) // cycle on all points found
  {
    double minDistance = INFINITY;
    for ( int i = 0 ; i < nfound ; i++) // cycle on all the other points
    {
      if(i != j)
      {
	//calculate the distance
	
	double distance = sqrt( pow((xpeaks[j] - xpeaks[i]),2) +   pow(ypeaks[j] - ypeaks[i],2) );
	if(distance < minDistance)
	  minDistance = distance;
	// 	std::cout << j << " " << i << " "  << xpeaks[j] << " "  << ypeaks[j]  << " "  << xpeaks[i]<< " "  << ypeaks[i] << " "  << distance << " " << minDistance << std::endl;
      }
    }
    xsigma[j] = ysigma[j] = minDistance / 2.0;
  }
  
  
  // manually set the sigma for the fitting area
  //   for(int i = 0 ; i < nofcrystals ; i++)
  //   {
  //     tempXsigma[i] = 0.2;
  //     tempYsigma[i] = 0.2;
  //     xsigma[i] = tempXsigma[i];
  //     ysigma[i] = tempYsigma[i];
  //   }
  
  // zero level. just use the peak positions and separation ranges from half distance
  //   for ( int j = 0 ; j < nfound ; j++)
  //   {
  //     fit2DmeanX.push_back(xpeaks[j]);
  //     fit2DmeanY.push_back(ypeaks[j]);
  //     fit2DsigmaX.push_back(xsigma[j]);
  //     fit2DsigmaY.push_back(ysigma[j]);
  //     fit2Dtheta.push_back(0);
  //   }
  
  
  
  
  //2d fit
  for ( int j = 0 ; j < nfound ; j++)
  {
    //find the peak height (to use it as amplitude parameter)
    Double_t zpeaks = FloodMap2D.GetBinContent(FloodMap2D.GetXaxis()->FindBin(xpeaks[j]) ,FloodMap2D.GetYaxis()->FindBin(ypeaks[j]));
    const Int_t npar = 6;
    Double_t f2params[npar] = {zpeaks,xpeaks[j],ypeaks[j],xsigma[j],ysigma[j],0};
    //  f2params[] = {amplitude,mean_x,mean_y,sigma_x,sigma_y,theta}
    //declare the fit function
    TF2 *f2 = new TF2("f2",gauss2d,xpeaks[j]-xsigma[j],xpeaks[j]+xsigma[j],ypeaks[j]-ysigma[j],ypeaks[j]+ysigma[j], npar);
    f2->SetParameters(f2params);  // set initial guesses
    
    //parameter limits
    //f2->SetParLimits(0,0,INFINITY);
    f2->SetParLimits(1,xpeaks[j]-0.5*xsigma[j],xpeaks[j]+0.5*xsigma[j]);
    f2->SetParLimits(2,ypeaks[j]-0.5*ysigma[j],ypeaks[j]+0.5*ysigma[j]);
    f2->SetParLimits(3,-xsigma[j],xsigma[j]);
    f2->SetParLimits(4,-ysigma[j],ysigma[j]);
    f2->SetParLimits(5,0,TMath::Pi());
    
    //     TF2 *f2 = new TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",xpeaks[j]-2.0*xsigma[j],xpeaks[j]+2.0*xsigma[j],ypeaks[j]-2.0*ysigma[j],ypeaks[j]+2.0*ysigma[j]);
    //     f2->SetParameters(14,xpeaks[j],xsigma[j],ypeaks[j],ysigma[j]);
    FloodMap2D.Fit("f2","QNR");
    fit2DmeanX.push_back(f2->GetParameter(1));
    fit2DmeanY.push_back(f2->GetParameter(2));
    //we insert the 2.5 sigma limit here, it makes things easier
    fit2DsigmaX.push_back(TMath::Abs(f2->GetParameter(3)));
    fit2DsigmaY.push_back(TMath::Abs(f2->GetParameter(4)));
    //careful: theta will set t, and t is in degrees, not in radiands!!
    fit2Dtheta.push_back(180.0*f2->GetParameter(5)/TMath::Pi()); 
  }  
  
  
  //forget about ellipses
  //   for ( int j = 0 ; j < nfound ; j++)
  //   {
  //     float sx,sy;
  //     sx = fit2DsigmaX[j];
  //     sy = fit2DsigmaY[j];
  //     fit2DsigmaX[j] = (sx+sy)/2.0;
  //     fit2DsigmaY[j] = (sx+sy)/2.0;
  //   }
  
  //   now for each peak, check if any other peak is closer than the sum of the relative circles
  //   for ( int j = 0 ; j < nfound ; j++)// run on all peaks
  //   {
  //     for ( int jOther = 0 ; jOther < nfound ; jOther++)//run on all peaks again, but...
  //     {
  //       if (jOther != j) //...do it only if it's not the same peak
  //       {
  // 	float distance = TMath::Sqrt( TMath::Power(fit2DmeanX[j]-fit2DmeanX[jOther],2) + TMath::Power(fit2DmeanY[j]-fit2DmeanY[jOther],2) );
  // 	float sumOfRadii = (fit2DsigmaX[j]+fit2DsigmaX[jOther]);
  // 	if ( distance < sumOfRadii )
  // 	{
  // 	  //std::cout << "WARNING: Peaks of Module " << k << " Channel " << i << " are overlapping!" << std::endl;
  // 	  fit2DsigmaX[j] = distance * ( fit2DsigmaX[j] / sumOfRadii );
  // 	  fit2DsigmaX[jOther] = distance * ( fit2DsigmaX[jOther] / sumOfRadii );
  // 	  fit2DsigmaY[j] = fit2DsigmaX[j];
  // 	  fit2DsigmaY[jOther] = fit2DsigmaX[jOther];
  // 	}
  //       }
  //     }
  //   }
  //   
  //   // sort them so that first is lower y
  //   double swapFit2DmeanX,swapFit2DmeanY,swapFit2DsigmaX,swapFit2DsigmaY;
  //   if(fit2DmeanY[0] > fit2DmeanY[1])
  //   {
  //     swapFit2DmeanX  = fit2DmeanX[1];
  //     swapFit2DmeanY  = fit2DmeanY[1];
  //     swapFit2DsigmaX = fit2DsigmaX[1];
  //     swapFit2DsigmaY = fit2DsigmaY[1];
  //     
  //     fit2DmeanX[1] = fit2DmeanX[0];
  //     fit2DmeanY[1] = fit2DmeanY[0];
  //     fit2DsigmaX[1] = fit2DsigmaX[0];
  //     fit2DsigmaY[1] = fit2DsigmaY[0];
  //     
  //     fit2DmeanX[0] = swapFit2DmeanX  ;
  //     fit2DmeanY[0] = swapFit2DmeanY  ;
  //     fit2DsigmaX[0] =swapFit2DsigmaX ;
  //     fit2DsigmaY[0] =swapFit2DsigmaY ;
  //   }
  
  return nfound;
}


void Mppc::FindProjectionPlane()
{
  TString name = FloodMap3D.GetName();
  //   std::cout << FloodMap3D.GetName() << std::endl;
  
  TString string_zy = name + "_zy";
  TString string_zx = name + "_zx";
  
  FloodMap3D.Project3D("zy");  // projection of u,v,w points on v,w plane 
  projection_zy = (TH2D*) gDirectory->Get(string_zy);
  
  FloodMap3D.Project3D("zx");  // projection of u,v,w points on u,w plane
  projection_zx = (TH2D*) gDirectory->Get(string_zx);
  
  //   projection_y = projection_zy->ProjectionY("spectrum3d_x"); 
  //   projection_x = projection_zx->ProjectionY("spectrum3d_y");
  
  profileX = projection_zx->ProfileY();  // profile along w of u,w plot. for each bin in w, average u is plotted with stdev 
  profileY = projection_zy->ProfileY();  // profile along w of u,w plot. for each bin in w, average v is plotted with stdev 
  
  lineX = new TF1("lineX",  "[0]*x + [1]",0,1);  // line becomes u(w) = m*w + c
  lineY = new TF1("lineY",  "[0]*x + [1]",0,1);  // line becomes v(w) = m*w + c
  
  profileX->Fit("lineX","Q");
  profileY->Fit("lineY","Q");
  
  //FIXME
  ThetaWU = TMath::ATan(lineX->GetParameter(0)); // --> ThetaWU = atan(-1/(m_w,u)) --> angle perpendicular to the angle from w to u in radiants
  ThetaWV = TMath::ATan(lineY->GetParameter(0)); // --> ThetaWV = atan(-1/(m_w,v)) --> angle perpendicular to the angle from w to v in radiants
  
  //equations become 
  // u" = u*cos(t_wu) + (v*sin(t_wv) + w*cos(t_wv))*sin(t_wu) 
  // v" = v * cos(t_wv) - w*sin(t_wv)
  
  
}


void Mppc::MakeRotatedFlood()
{
  //-------------------------------------------------------------------------------  
  // Flood histogram
  // Plot a different histogram depending on the position of the mppc
  std::stringstream varX,varY; // the variables of the following 2d plots will have to be build custom depending on the position of the mppc
  
  varX << "FloodX";
  varY << "FloodY";
  
  
  
  
  // the if statements below are an embarassing example of how poor my coding is. But hei, i'm in a rush for a conference..
//   spectrum2d = new TH2F("spectrum2d","spectrum2d",histo2DchannelBin,-7,7,histo2DchannelBin,-7,7);
//   if( ((iModule*nmppcx)+iMppc) > 0 && (((iModule*nmppcx)+iMppc) < nmppcx -1) && ((jModule*nmppcy)+jMppc) > 0 && (((jModule*nmppcy)+jMppc) < nmppcy -1 )) // central mppcs
//   {
//     // standard flood 2d
//     varX << "FloodX";
//     varY << "FloodY";
//     // 	    std::cout << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " CENTRAL"<< std::endl; 
//   }
//   else // we are not in the center
//   {
//     if( ((iModule*nmppcx)+iMppc) > 0 && (((iModule*nmppcx)+iMppc) < nmppcx -1) ) // we are in top or bottom lateral
//     {
//       if(((jModule*nmppcy)+jMppc) == 0 )   // bottom lateral crystals
//       {
// 	lateralQ2         = - base_lateralQ2;
// 	lateralDeltaV     = - base_lateralDeltaV;
// 	lateralRescaleTB  = + base_lateralRescaleTB;
// 	// 		std::cout << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " BOTTOM LATERAL"<< std::endl;
//       }
//       if(((jModule*nmppcy)+jMppc) == nmppcy -1 )   // top lateral crystals
//       {
// 	lateralQ2         = + base_lateralQ2;
// 	lateralDeltaV     = + base_lateralDeltaV;
// 	lateralRescaleTB  = + base_lateralRescaleTB;
// 	// 		std::cout << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " TOP LATERAL"<< std::endl;
//       }
//       // left and right crystals have the same var structure
//       varY << "("
//       << lateralDeltaV
//       << " + "
//       << lateralRescaleTB
//       << "* ((FloodY)*TMath::Cos(" 
//       << lateralQ2
//       << ") - (FloodZ)*TMath::Sin("
//       << lateralQ2
//       << ")))";
//       varX << "(FloodX)";  
//     }
//     else
//     {
//       if(((jModule*nmppcy)+jMppc) > 0 && (((jModule*nmppcy)+jMppc) < nmppcy -1 ) )  // we are in left or right lateral
//       {
// 	if(((iModule*nmppcx)+iMppc) == 0 )   // left lateral crystals
// 	{
// 	  lateralQ1         = + base_lateralQ1;
// 	  lateralDeltaU     = - base_lateralDeltaU;
// 	  lateralRescaleRL  = + base_lateralRescaleRL;
// 	  // 		  std::cout << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " LEFT LATERAL "<< std::endl;
// 	}
// 	if(((iModule*nmppcx)+iMppc) == nmppcx -1 )   // right lateral crystals
// 	{
// 	  lateralQ1         = - base_lateralQ1;
// 	  lateralDeltaU     = + base_lateralDeltaU;
// 	  lateralRescaleRL  = + base_lateralRescaleRL;
// 	  // 		  std::cout << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " RIGHT LATERAL" << std::endl;
// 	}
// 	// left and right crystals have the same var structure
// 	varY << "(FloodY)";  
// 	varX << "("
// 	<< lateralDeltaU
// 	<< " + "
// 	<< lateralRescaleRL
// 	<< "*(FloodZ*TMath::Sin("
// 	<< lateralQ1
// 	<< ")+ (FloodX)*TMath::Cos("
// 	<< lateralQ1
// 	<<")))";
//       }
//       else // only corner crystals remain..
//       {
// 	if(((iModule*nmppcx)+iMppc) == 0 &&  ((jModule*nmppcy)+jMppc) == 0 )   // bottom left crystals
// 	{
// 	  cornerDeltaU   = - base_cornerDeltaU  ;
// 	  cornerDeltaV   = - base_cornerDeltaV  ; 
// 	  cornerQ1       = + base_cornerQ1      ;
// 	  cornerQ2       = - base_cornerQ2      ;
// 	  cornerRescale  =   base_cornerRescale ;
// 	  // 		  std::cout << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " BOTTOM LEFT"<< std::endl;
// 	}
// 	if(((iModule*nmppcx)+iMppc) == nmppcx -1 &&  ((jModule*nmppcy)+jMppc) == 0 )   // bottom right crystals
// 	{
// 	  cornerDeltaU   = + base_cornerDeltaU  ;
// 	  cornerDeltaV   = - base_cornerDeltaV  ;
// 	  cornerQ1       = - base_cornerQ1      ;
// 	  cornerQ2       = - base_cornerQ2      ;
// 	  cornerRescale  =   base_cornerRescale ;
// 	  // 		  std::cout << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " BOTTOM RIGTH"<< std::endl;
// 	}
// 	if(((iModule*nmppcx)+iMppc) == 0 &&  ((jModule*nmppcy)+jMppc) == nmppcy -1 )   // top left crystals
// 	{
// 	  cornerDeltaU   = - base_cornerDeltaU  ;
// 	  cornerDeltaV   = + base_cornerDeltaV  ;
// 	  cornerQ1       = - base_cornerQ1      ;
// 	  cornerQ2       = + base_cornerQ2      ;
// 	  cornerRescale  =   base_cornerRescale ;
// 	  // 		  std::cout << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " TOP LEFT"<< std::endl;
// 	}
// 	if(((iModule*nmppcx)+iMppc) == nmppcx -1 &&  ((jModule*nmppcy)+jMppc) == nmppcy -1)   // top right crystals
// 	{
// 	  cornerDeltaU   = + base_cornerDeltaU  ;
// 	  cornerDeltaV   = + base_cornerDeltaV  ;
// 	  cornerQ1       = + base_cornerQ1      ;
// 	  cornerQ2       = + base_cornerQ2      ;
// 	  cornerRescale  =   base_cornerRescale ;
// 	  // 		  std::cout << mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel() << " TOP RIGHT"<< std::endl;
// 	}
// 	varY 
// 	<< "( "
// 	<< cornerDeltaV
// 	<< " + " 
// 	<< cornerRescale
// 	<<"*((FloodX*TMath::Sin(" 
// 	<< cornerQ1 
// 	<< ") + FloodY*TMath::Cos( "
// 	<< cornerQ1 
// 	<< ")) *TMath::Cos("
// 	<< cornerQ2
// 	<< ") - FloodZ*TMath::Sin("
// 	<< cornerQ2 
// 	<< ")))" ;
// 	varX << " ( " 
// 	<< cornerDeltaU
// 	<< " + (FloodX*TMath::Cos("
// 	<< cornerQ1
// 	<< ") - FloodY*TMath::Sin("
// 	<< cornerQ1
// 	<< ")))";
// 	
// 	// 		std::cout << varY.str() << std::endl;
// 	// 		std::cout << varX.str() << std::endl;
//       }
//     }	    
//   }
//   var << varY.str() << ":" << varX.str() << " >> spectrum2d";
//   std::cout << var << std::endl;
//   tree->Draw(var.str().c_str(),CutXYZ+CutTrigger,"COLZ");
//   name = "Flood Histogram 2D - MPPC " + mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->GetLabel();
//   spectrum2d->SetName(name); 
//   spectrum2d->SetTitle(name);
//   spectrum2d->GetXaxis()->SetTitle("U");
//   spectrum2d->GetYaxis()->SetTitle("V");
//   mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetFloodMap2D(*spectrum2d);
//   mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetXvariable(varX.str());
//   mppc[(iModule*nmppcx)+iMppc][(jModule*nmppcy)+jMppc]->SetYvariable(varY.str());
//   /*TSpectrum2 *peaks2D = new TSpectrum2(ncrystalsx*ncrystalsy,1);
//    *  int nfound2D = peaks2D->Search(spectrum2d,1,"col",0.3);*/	
//   module[iModule][jModule]->GetFloodMap2DSeparated()->Add(spectrum2d);
//   // 	  module[iModule][jModule]->SetFloodMap2DSeparated(*spectrum2d);
//   varX.str("");
//   varY.str("");
//   var.str("");
//   delete spectrum2d; 
  //-------------------------------------------------------------------------------
  
}




void Mppc::PrintSpecific()
{
  std::cout << "Label \t\t = " << label << std::endl;
  std::cout << "Digitizer Ch. \t = " << digitizerChannel << std::endl;
  std::cout << "Staus for Doi    \t = " << IsOnForDoi << std::endl;
}