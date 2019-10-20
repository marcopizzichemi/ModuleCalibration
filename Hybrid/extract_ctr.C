

void extract_ctr(std::string fileName = "resB3.root", int crystal = 9 ,int rebin = 20, double z_fwhm = 2.0)
{
  TFile *_file0 = TFile::Open(fileName.c_str());

  double z_sigma = z_fwhm/2.355;
  // basic CTR
  std::stringstream sname;
  sname << "Basic CTR vs. Z - Crystal " << crystal ;
  TH2F *basicCTRvsZ = (TH2F*) gDirectory->Get(sname.str().c_str());
  sname.str("");
  basicCTRvsZ->RebinX(rebin);

  basicCTRvsZ->FitSlicesY(0, 0, -1, 0, "QNR");
  sname << basicCTRvsZ->GetName() << "_1";
  TH1D *spectrum2d_1 = (TH1D*)gDirectory->Get(sname.str().c_str()); // _1 is the TH1D automatically created by ROOT when FitSlicesY is called, holding the TH1F of the mean values
  sname.str("");
  new TCanvas;
  spectrum2d_1->Draw();

  sname << basicCTRvsZ->GetName() << "_2";
  TH1D *spectrum2d_2 = (TH1D*)gDirectory->Get(sname.str().c_str()); // _2 is the TH1D automatically created by ROOT when FitSlicesY is called, holding the TH1F of the sigma values
  sname.str("");

  new TCanvas;
  spectrum2d_2->Draw();

  int i = 0;
  std::vector<double> x,y,ex,ey;
  for(i = 0; i < spectrum2d_2->GetNbinsX(); i++)
  {
    x.push_back(spectrum2d_2->GetBinCenter(i+1));
    y.push_back(spectrum2d_2->GetBinContent(i+1));
    ex.push_back(z_sigma);
    ey.push_back(2.5);
  }

  for(i = 0; i < x.size(); i++)
  {
    y[i] = 1e12*sqrt(2)*sqrt(pow(y[i]*2.355,2)-pow(90e-12,2));
  }

  TGraphErrors *g = new TGraphErrors(x.size(),&x[0],&y[0],&ex[0],&ey[0]);
  new TCanvas;
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlue);

  g->GetXaxis()->SetTitle("Interaction Position in crystal (0 = MPPC side) [mm]");
  g->GetYaxis()->SetTitle("CTR [ps]");
  g->GetXaxis()->SetRangeUser(0,15);
  g->GetYaxis()->SetRangeUser(140,250);
  g->Draw("AP");


  // full CTR
  sname << "Full CTR vs. Z - Crystal " << crystal;
  TH2F *fullCTRvsZ = (TH2F*) gDirectory->Get(sname.str().c_str());
  sname.str("");
  fullCTRvsZ->RebinX(rebin);
  // std::stringstream sname;
  fullCTRvsZ->FitSlicesY(0, 0, -1, 0, "QNR");
  sname << fullCTRvsZ->GetName() << "_1";
  TH1D *fspectrum2d_1 = (TH1D*)gDirectory->Get(sname.str().c_str()); // _1 is the TH1D automatically created by ROOT when FitSlicesY is called, holding the TH1F of the mean values
  sname.str("");
  new TCanvas;
  fspectrum2d_1->Draw();

  sname << fullCTRvsZ->GetName() << "_2";
  TH1D *fspectrum2d_2 = (TH1D*)gDirectory->Get(sname.str().c_str()); // _2 is the TH1D automatically created by ROOT when FitSlicesY is called, holding the TH1F of the sigma values
  sname.str("");

  new TCanvas;
  fspectrum2d_2->Draw();

  // int i = 0;
  std::vector<double> xf,yf,exf,eyf;
  for(i = 0; i < fspectrum2d_2->GetNbinsX(); i++)
  {
    xf.push_back(fspectrum2d_2->GetBinCenter(i+1));
    yf.push_back(fspectrum2d_2->GetBinContent(i+1));
    exf.push_back(z_sigma);
    eyf.push_back(2.5);
  }

  for(i = 0; i < xf.size(); i++)
  {
    yf[i] = 1e12*sqrt(2)*sqrt(pow(yf[i]*2.355,2)-pow(90e-12,2));
  }

  TGraphErrors *gf = new TGraphErrors(xf.size(),&xf[0],&yf[0],&exf[0],&eyf[0]);
  new TCanvas;
  gf->SetMarkerStyle(20);
  gf->SetMarkerColor(kRed);
  gf->GetXaxis()->SetTitle("Interaction Position in crystal (0 = MPPC side) [mm]");
  gf->GetYaxis()->SetTitle("CTR [ps]");
  gf->GetXaxis()->SetRangeUser(0,15);
  gf->GetYaxis()->SetRangeUser(140,250);
  gf->Draw("AP");

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(g);
  mg->Add(gf);

  new TCanvas;
  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("Interaction Position in crystal (0 = MPPC side) [mm]");
  mg->GetYaxis()->SetTitle("CTR FWHM [ps]");
  mg->GetXaxis()->SetRangeUser(0,15);
  mg->GetYaxis()->SetRangeUser(130,240);

  TLegend *legend = new TLegend(0.11,0.65,0.5,0.89,"");
  legend->SetFillStyle(0);
  legend->AddEntry(g,"Standard CTR","p");
  legend->AddEntry(gf,"Fully corrected CTR","p");
  legend->Draw();
  // mg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  // mg->GetYaxis()->SetTitle("Coefficients");
  gPad->Modified();



}
