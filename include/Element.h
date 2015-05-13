// parent class of the elements in the detector
// module, mppc and crystal inheritate from it

#ifndef ELEMENT_H
#define ELEMENT_H

// #include <string>

class Element
{
protected:
  std::string          name;        // name or label
  int                  id;          // id number
  float                x,y,z;       // position coordinate space
  
  
public:
  
  Element(); // default constructor
  Element(std::string aname, int pid, float px, float py, float pz); //constructor
  Element(const Element &obj); // copy constructor
  ~Element(); // destructor
  
  std::string          GetName(){return name;};
  int                  GetID()  {return id;};
  float                GetX()   {return x;};
  float                GetY()   {return y;};
  float                GetZ()   {return z;};
  
  void                 SetName(std::string aname)               {name = aname;};
  void                 SetID(int pid)                           {id = pid;};
  void                 SetPosition(float px, float py, float pz){x = px; y = py; z = pz;};
  
};




#endif  // ELEMENT_H