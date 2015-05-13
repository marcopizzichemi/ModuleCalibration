#ifndef MPPC_H
#define MPPC_H
#include "Element.h"
#include "Module.h"
// #include <string>

class Mppc : public Element
{
private:
  Module*    parentModule    = NULL;           // one pointer for its parent, it's only one
//   Crystal*** childrenCrystal = NULL;           // and array of pointers for its children
  
public:
  Mppc(); // default constructor
  //Mppc(std::string aname, int pid, float px, float py, float pz); //constructor
  Mppc(const Mppc &obj); // copy constructor
  ~Mppc(); // destructor
  
  Module*                GetModule(){return parentModule;};
  void                   SetModule(Module* amodule){parentModule = amodule;}; 
//   void                   MakeCrystalsPointers(int i, int j);
//   Crystal*               GetChildrenCrystal(int i, int j){return childrenCrystal[i][j];};
//   void                   SetChildrenCrystal(int i, int j, Crystal* acrystal){childrenCrystal[i][j] = acrystal;}
  
  void PrintGlobal();
  void PrintSpecific();
};




#endif  // MPPC_H