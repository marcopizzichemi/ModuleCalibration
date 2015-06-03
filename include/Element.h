// parent class of the elements in the detector
// module, mppc and crystal inheritate from it

#ifndef ELEMENT_H
#define ELEMENT_H
#include <sstream>
#include <iostream>
#include <vector>
#include "TH2F.h"
#include "TH3F.h"

class Module;
class Mppc;
class Crystal;

class Element
{
protected:
  std::string          name;              // name or label
  std::string          label; 
  std::string          parentName;
  std::vector<std::string> childrenName;
  int                  id;                // id number
  int                  i,j;               // i and j IDs
  float                x,y,z;             // position coordinate space
//   std::string          GlobalTag;         // tag identifing the element
                                          // structure is module-mppc-crystal
                                          // by default is set to x-x.x-x.x
  int                  iChildren;         // number of children on "i"
  int                  jChildren;         // number of children on "j"
  
  //2d histos
  TH2F                 FloodMap2D;
  TH2F                 SphericalMap;
  TH2F                 CylindricalXMap;
  TH2F                 CylindricalYMap;
  //3d histos
  TH3F                 FloodMap3D;
  
  
  
public:
  
  Element(); // default constructor
  Element(const Element &obj); // copy constructor
  ~Element(){}; // destructor
  
  std::string          GetName()                                 {return name;};
  std::string          GetLabel()                                {return label;};
  int                  GetID()                                   {return id;};
  int                  GetI()                                    {return i;};
  int                  GetJ()                                    {return j;};
  float                GetX()                                    {return x;};
  float                GetY()                                    {return y;};
  float                GetZ()                                    {return z;};
  int                  GetChildrenI()                            {return iChildren;};
  int                  GetChildrenJ()                            {return jChildren;};
  void                 SetName(std::string aname)                {name = aname;};
  void                 SetLabel(std::string aname)               {label = aname;};
  void                 SetID(int pid)                            {id = pid;};
  void                 SetI(int pi)                              {i = pi;};
  void                 SetJ(int pj)                              {j = pj;};
  void                 SetPosition(float px, float py, float pz) {x = px; y = py; z = pz;};
  void                 SetChildrenI(int pi)                      {iChildren = pi;};
  void                 SetChildrenJ(int pj)                      {jChildren = pj;};
  
//   void                 SetGlobalTag(int module, int mppcx, int mppcy, int cryx , int cryy);
  void                 PrintGlobal();
  virtual void         PrintSpecific();
  void                 Print(){std::cout<<std::endl;PrintGlobal(); PrintSpecific();std::cout<<std::endl;};
  
  void                 SetParentName(std::string aName){parentName = aName;};
  std::string          GetParentName(){return parentName;};
  void                 AddChild(std::string aName){childrenName.push_back(aName);};
  std::vector<std::string> GetChildren(){return childrenName;};
  
  
  TH2F*                GetFloodMap2D(){return &FloodMap2D;};
  TH2F*                GetSphericalMap(){return &SphericalMap;};
  TH2F*                GetCylindricalXMap(){return &CylindricalXMap;};
  TH2F*                GetCylindricalYMap(){return &CylindricalYMap;};
  
  TH3F*                GetFloodMap3D(){return &FloodMap3D;};
  
  
  void                 SetFloodMap2D(TH2F aHisto){FloodMap2D = aHisto;};
  void                 SetSphericalMap(TH2F aHisto){SphericalMap = aHisto;};
  void                 SetCylindricalXMap(TH2F aHisto){CylindricalXMap = aHisto;};
  void                 SetCylindricalYMap(TH2F aHisto){CylindricalYMap = aHisto;};
  
  void                 SetFloodMap3D(TH3F aHisto){FloodMap3D = aHisto;};
  
  
  
  
  
};








#endif  // ELEMENT_H