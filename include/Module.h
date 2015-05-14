#ifndef MODULE_H
#define MODULE_H

#include <string>
#include "Element.h"

class Module : public Element
{
// private: 
//   Element*** childrenMppc; // matrix of children mppcs
  
public:
  Module(); // default constructor
  //Module(std::string aname, int pid, float px, float py, float pz); //constructor
  Module(const Module &obj); // copy constructor
  ~Module(); // destructor
    
//   Element*               GetChildrenMppc(int i, int j){return childrenMppc[i][j]};
//   void                   MakeMppcPointers(int i, int j);
//   void                   SetChildrenMppc(int i, int j, Mppc* amppc){childrenMppc[i][j] = amppc;}
  
  Element* GetMppc(int i , int j){return GetChild(i,j);};
  Element* GetMppc(int digiCh);
  Element* GetMppc(std::string alabel);
  Element* GetCrystal(int iMppc , int jMppc , int iCrystal , int jCrystal){return GetChild(iMppc,jMppc)->GetChild(iCrystal,jCrystal);};
  
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // MODULE_H