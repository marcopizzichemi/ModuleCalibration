#include <iostream>
#include "Element.h"


Element::Element()
{
  //default constructor
  name = "Default Element";
  id = 0;
  x = y = z = 0;
}

Element::Element(std::string aname, int pid, float px, float py, float pz)
{
  //constructor
  name = aname;
  id = pid;
  x = px;
  y = py;
  z = pz;
}

Element::Element(const Element &obj) 
{
  //copy constructor
}

Element::~Element()
{
  //destructor
}

