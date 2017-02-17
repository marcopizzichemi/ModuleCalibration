#ifndef ELEMENT_H
#define ELEMENT_H
#include <sstream>
#include <iostream>
#include <vector>
#include "TH2F.h"
#include "TH3I.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TProfile.h"
#include "TTree.h"

class Module;
class Mppc;
class Crystal;


/** @brief Generic element
 *
 *  Parent class of the elements in the detector
 *  Module, mppc and crystal inheritate from here
 *
 *  @author M. Pizzichemi
 *  @date Nov 2015
 */

class Element
{
protected:
  // general variables of an element

  std::string          name;              ///< name
  std::string          label;             ///< label
  std::string          parentName;        ///< name of parent element
  std::vector<std::string> childrenName;  ///< name of children elements
  std::string          extendedID;        ///< sort of unique identifier --> moduleI.moduleJ.mppcI.mppcJ.crystalI.crystalJ
  std::string          variableX;         ///< x variable for the specific roto-translated plot of this element
  std::string          variableY;         ///< y variable for the specific roto-translated plot of this element
  int                  id;                ///< id number
  int                  i;                 ///< i ID
  int                  j;                 ///< i ID
  float                x;                 ///< x position coordinate space [mm]
  float                y;                 ///< y position coordinate space [mm]
  float                z;                 ///< z position coordinate space [mm]
  float                dx;                ///< x dimensions of the element [mm]
  float                dy;                ///< y dimensions of the element [mm]
  float                dz;                ///< z dimensions of the element [mm]
  int                  iChildren;         ///< number of children on "i"
  int                  jChildren;         ///< number of children on "j"
  bool                 isOnForModular;
  TH1D*                SlicesSigma;           ///< histogram of fitted sigma values of profiles from TH2F EnDepSumCharge
  TGraph*              EnDepSumChargeGraph;
  //2d histos
  TH2F*                FloodMap2D;          ///< u,v map for this element
  TH2F*                FloodMap2DSingleCrystalHit; ///< Simulation dataset: u,v map for this element, accepting only events where energy was deposited in just one crystal
  TH2F*                ADCversusW;          ///< 2d histo of "energy" versus w coordinate
  TH2F*                ADCversusWComplete;  ///< 2d histo of "energy" versus w coordinate - without cuts
  //3d histos
  TH3I*                FloodMap3D;          ///< u,v,w map for this element
  TH3I*                FloodMap3DSeparation;

public:

  Element();                                ///< default constructor
  Element(const Element &obj);              ///< copy constructor
  ~Element(){};                             ///< destructor

  // methods to get and set the private variables. Names should be self explanatory
  std::string          GetName()                                 {return name;};
  std::string          GetLabel()                                {return label;};
  int                  GetID()                                   {return id;};
  std::string          GetExtendedID()                           {return extendedID;};
  int                  GetI()                                    {return i;};
  int                  GetJ()                                    {return j;};
  float                GetX()                                    {return x;};
  float                GetY()                                    {return y;};
  float                GetZ()                                    {return z;};
  float                GetDimensionX()                           {return dx;};
  float                GetDimensionY()                           {return dy;};
  float                GetDimensionZ()                           {return dz;};
  int                  GetChildrenI()                            {return iChildren;};
  int                  GetChildrenJ()                            {return jChildren;};
  std::string          GetParentName()                           {return parentName;};
  TH2F*                GetFloodMap2D()                           {return FloodMap2D;};
  TH2F*                GetFloodMap2DSingleCrystalHit()           {return FloodMap2DSingleCrystalHit;};
  TH2F*                GetADCversusW()                           {return ADCversusW;};
  TH2F*                GetADCversusWComplete()                   {return ADCversusWComplete;};
  //   TH2F*                GetFloodMap2DSeparated()                  {return &FloodMap2DSeparated;};
  //   TH2F*                GetSphericalMap()                         {return &SphericalMap;};
  //   TH2F*                GetCylindricalXMap()                      {return &CylindricalXMap;};
  //   TH2F*                GetCylindricalYMap()                      {return &CylindricalYMap;};
  //   TH2F*                GetLateralMap()                           {return &LateralMap;};
  //   TH2F*                GetCornerMap()                            {return &CornerMap;};
  //   TH2F*                GetCentralMap()                           {return &CentralMap;};
  TH3I*                GetFloodMap3D()                           {return FloodMap3D;};
  TH3I*                GetFloodMap3DSeparation()                           {return FloodMap3DSeparation;};
  std::string          GetXvariable()                            {return variableX;};
  std::string          GetYvariable()                            {return variableY;};
  bool                 GetIsOnForModular()                       {return isOnForModular;};
  //   TGraph*              GetGraphFlood2D()                         {return &GraphFlood2D;};
  //   TGraph2D*            GetGraphFlood3D()                         {return &GraphFlood3D;};
  TH1D*                GetSlicesSigma()                          {return SlicesSigma;};
  TGraph*              GetEnDepSumChargeGraph()                  {return EnDepSumChargeGraph;};

  void                 SetName(std::string aname)                {name = aname;};
  void                 SetLabel(std::string aname)               {label = aname;};
  void                 SetID(int pid)                            {id = pid;};
  void                 SetExtendedID(std::string pid)            {extendedID = pid;};
  void                 SetI(int pi)                              {i = pi;};
  void                 SetJ(int pj)                              {j = pj;};
  void                 SetPosition(float px, float py, float pz) {x = px; y = py; z = pz;};
  void                 SetDimension(float px, float py, float pz){dx = px; dy = py; dz = pz;};
  void                 SetChildrenI(int pi)                      {iChildren = pi;};
  void                 SetChildrenJ(int pj)                      {jChildren = pj;};
  void                 SetParentName(std::string aName)          {parentName = aName;};
  void                 SetFloodMap2D(TH2F* aHisto)                {FloodMap2D = aHisto;};
  void                 SetFloodMap2DSingleCrystalHit(TH2F* aHisto)  {FloodMap2DSingleCrystalHit = aHisto;};
  void                 SetADCversusW(TH2F* aHisto)                {ADCversusW = aHisto;};
  void                 SetADCversusWComplete(TH2F* aHisto)        {ADCversusWComplete = aHisto;};
  void                 SetIsOnForModular(bool aBool)              {isOnForModular = aBool;};
  //   void                 SetFloodMap2DSeparated(TH2F aHisto)       {FloodMap2DSeparated = aHisto;};
  //   void                 SetSphericalMap(TH2F aHisto)              {SphericalMap = aHisto;};
  //   void                 SetCylindricalXMap(TH2F aHisto)           {CylindricalXMap = aHisto;};
  //   void                 SetCylindricalYMap(TH2F aHisto)           {CylindricalYMap = aHisto;};
  void                 SetFloodMap3D(TH3I* aHisto)                {FloodMap3D = aHisto;};
  void                 SetFloodMap3DSeparation(TH3I* aHisto)                {FloodMap3DSeparation = aHisto;};
  //   void                 SetLateralMap(TH2F aHisto)                {LateralMap = aHisto;};
  //   void                 SetCornerMap(TH2F aHisto)                 {CornerMap = aHisto;};
  //   void                 SetCentralMap(TH2F aHisto)                {CentralMap = aHisto;};
  void                 SetXvariable(std::string astring)         {variableX = astring;};
  void                 SetYvariable(std::string astring)         {variableY = astring;};
  //   void                 SetGraphFlood2D(TGraph aGraph)            {GraphFlood2D = aGraph;};
  //   void                 SetGraphFlood3D(TGraph2D aGraph)          {GraphFlood3D = aGraph;};

  //methods to add and return children elements
  void                 AddChild(std::string aName)               {childrenName.push_back(aName);};
  void                 SetSlicesSigma(TH1D* aHisto)              {SlicesSigma = aHisto;};
  void                 SetEnDepSumChargeGraph(TGraph* aGraph)    {EnDepSumChargeGraph = aGraph;};


  std::vector<std::string> GetChildren()                         {return childrenName;};


  //   void SetTempMap2D(TH2F aHisto){TempMap2D = aHisto;};
  //   TH2F* GetTempMap2D(){return &TempMap2D;};

  int MakeFloodMap2D(TTree* tree);

  //   TH2F*                Generate2Dhistogram(TString var,TCut);

  // methods to print element information
  void                 PrintGlobal();        ///< prints global info
  virtual void         PrintSpecific();      ///< prints specific info of this element. polimorphic implementation in the specific class of each elements
  void                 Print()                                   {std::cout<<std::endl;PrintGlobal(); PrintSpecific();std::cout<<std::endl;}; ///< Prints info about this element


};

#endif  // ELEMENT_H
