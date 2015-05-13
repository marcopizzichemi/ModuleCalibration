// parent class of the elements in the detector
// module, mppc and crystal inheritate from it

#ifndef ELEMENT_H
#define ELEMENT_H
#include <sstream>
// #include <string>

class Element
{
protected:
  std::string          name;              // name or label
  int                  id;                // id number
  int                  i,j;               // i and j IDs
  float                x,y,z;             // position coordinate space
  std::string          GlobalTag;         // tag identifing the element
                                          // structure is module-mppc-crystal
                                          // by default is set to x-x.x-x.x
//   int                  parentID;          // ID of parent element
  Element*             pParent = NULL;
  Element***           pChild;
  
public:
  
  Element(); // default constructor
//   Element(std::string aname, int pid, float px, float py, float pz); //constructor
  Element(const Element &obj); // copy constructor
  ~Element(); // destructor
  
  std::string          GetName()                                 {return name;};
  std::string          GetGlobalTag()                            {return GlobalTag;};
  int                  GetID()                                   {return id;};
  int                  GetI()                                    {return i;};
  int                  GetJ()                                    {return j;};
  float                GetX()                                    {return x;};
  float                GetY()                                    {return y;};
  float                GetZ()                                    {return z;};
//   int                  GetParentID()                             {return parentID;};
  void                 SetName(std::string aname)                {name = aname;};
  void                 SetID(int pid)                            {id = pid;};
  void                 SetI(int pi)                              {i = pi;};
  void                 SetJ(int pj)                              {j = pj;};
  void                 SetPosition(float px, float py, float pz) {x = px; y = py; z = pz;};
//   void                 SetParentID(int pid)                      {parentID = pid;};
  
  void                 MakeChildrenPointers(int i, int j);
  Element*             GetChild(int i, int j){return pChild[i][j];};
  void                 SetChild(int i, int j, Element* achild){pChild[i][j] = achild;}
  
  void                 SetGlobalTag(int module, int mppcx, int mppcy, int cryx , int cryy);
  void                 PrintGlobal();
  virtual void         PrintSpecific();
  void                 Print(){PrintGlobal(); PrintSpecific();};
  
};




#endif  // ELEMENT_H