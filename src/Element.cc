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