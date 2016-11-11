//line 416

nameModule = "Total Charge no Cuts " + module[iModule][jModule]->GetName();
varModule << "(";
for(int k=0; k<ncrystalsy*ncrystalsx*nmppcy*nmppcx*nmodulex*nmoduley-1; k++)
{
varModule << "TotalCryEnergy[" << k << "]+";
}
varModule << "TotalCryEnergy[" << ncrystalsy*ncrystalsx*nmppcy*nmppcx*nmodulex*nmoduley-1 << "])" << ">>" << nameModule;
TH1F *TotalNoCuts = new TH1F(nameModule, nameModule, histo1Dbins, 0, 1);
tree->Draw(varModule.str().c_str());
TotalNoCuts->GetXaxis()->SetTitle("Charge Deposited [MeV]");
TotalNoCuts->GetYaxis()->SetTitle("N");
std::cout <<  TotalNoCuts->Integral() << std::endl;
module[iModule][jModule]->SetTotalNoCuts(TotalNoCuts);
varModule.str("");


