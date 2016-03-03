#include <iostream>
#include "Element.h"



Element::Element()
{
  /**Default constructor
     Set a standard name ("Defaut Element") 
     Set id = i = j = 0
     Set x = y = z = 0
  */
  name = "Default Element";
  id = i = j = 0;
  x = y = z = 0;
//   GlobalTag = "x.x-x.x-x.x";
//   parentID = 0;
  
}

Element::Element(const Element &obj) 
{
  /**Copy constructor
   */
}

// int Element::MakeFloodMap2D(TTree* tree)
// {
//   std::stringstream var,cut,sname; 
//   sname << "Flood Histogram 2D - " << this->GetName();
//   var << "FloodY:FloodX >> " << sname.str();
//   int histo2DglobalBins = 1000;
//   TH2F *spectrum2d = new TH2F(sname.str().c_str(),sname.str().c_str(),histo2DglobalBins,-7,7,histo2DglobalBins,-7,7);
//   tree->Draw(var.str().c_str() ,"","COLZ");
// //   nameModule = "Flood Histogram 2D - " + module[iModule][jModule]->GetName();
//   spectrum2d->SetName(sname.str().c_str()); 
//   spectrum2d->SetTitle(sname.str().c_str());
//   spectrum2d->GetXaxis()->SetTitle("U");
//   spectrum2d->GetYaxis()->SetTitle("V");
//   this->SetFloodMap2D(*spectrum2d);
// //   delete spectrum2dModule;
//   return 0;
// }


void Element::PrintGlobal()
{
  /**Prints Element info to terminal
   * Variables printed here are common to all families of elements
   */
  std::cout << "Name \t\t = "    << name << std::endl;
  std::cout << "ID \t\t = "      << id << std::endl;
  std::cout << "Extended ID \t = "      << extendedID << std::endl;
  std::cout << "i = \t\t = "      << i << std::endl;
  std::cout << "j = \t\t = "      << j << std::endl;
  std::cout << "Position \t = (" << x << "," << y << "," << z << ")" << std::endl;
//   std::cout << "GlobalTag \t = " << GlobalTag << std::endl;
  std::cout << "Parent \t\t = " << parentName << std::endl;
  for(int k = 0; k < childrenName.size() ; k++)
  {
    std::cout << "Child \t\t = " << childrenName[k] << std::endl;
  }
}

void Element::PrintSpecific()
{
  /**Prints Element info to terminal
   * Variables printed here are specific to this type of element
   */
//   std::cout << "Abstract Element, nothing specific to print" << std::endl;
}