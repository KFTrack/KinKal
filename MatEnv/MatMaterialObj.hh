//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatMaterialObj.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatMaterialObj (Transient version)
//      Header file
// 
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby     (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   June 18, 1998 - Talby : created
//-----------------------------------------------------------------------------
#ifndef MATMATERIALOBJ_HH
#define MATMATERIALOBJ_HH 1

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
// Mandatory include file

#include <string>
#include <vector>
#include <math.h>
namespace MatEnv {

  class MatMaterialObj
  {

    public:   

      // Constructor to create Material from a combination of Elements or Materials
      MatMaterialObj();
      MatMaterialObj(const std::string& name, double density,  
	  double Zeff, double Aeff, int nbrComp,
	  std::vector<int>* iflg, 
	  std::vector<double>* weight, 
	  std::vector<std::string>* compName, 
	  double radLength, double intLength, double refIndex,
	  double temperature, double pressure, 
	  const std::string& state, double Tcut=0.0);

      // Destructor
      virtual ~MatMaterialObj();

      // print constructor
      void print();

      // Copy constructor
      MatMaterialObj(const MatMaterialObj&);

      // Operators (Assignment op)
      MatMaterialObj& operator= ( const MatMaterialObj& );   
      int operator==( const MatMaterialObj& ) const;

      // Selectors

      std::string getName() const { return std::string(_matName); };              
      double getDensity() const { return _matDensity; };           
      double getZeff() const { return _matZeff; };  
      double getAeff() const { return _matAeff; };       
      int getNbrComp() const { return _nbrComp; };           
      int getIflg(size_t index) const { return _iflg[index]; }
      double getWeight(size_t index) const { return _weight[index]; };
      std::string getCompName(size_t index) const { return std::string(_compName[index]); };

      // These are Optionals parameters
      double getRadLength() const { return _radLength; };           
      double getIntLength() const {return _intLength; }            
      double getRefIndex() const { return _refIndex; };  
      double getTemperature() const { return _matTemperature; };           
      double getPressure() const { return _matPressure; };  
      std::string getState() const { return std::string(_matState); };              
      double getTcut() const { return _matTcut; };

      void setName(const std::string Name) {_matName=Name; };
      void setDensity(double Density) {_matDensity=Density; };
      void setZeff(double Zeff) {_matZeff=Zeff; };  
      void setAeff(double Aeff) {_matAeff=Aeff; };  
      void setNbrComp(int NbrComp) {_nbrComp=NbrComp; };
      void setIflg(int Iflg) {_iflg.push_back(Iflg); };
      void setWeight(double Weight) {_weight.push_back(Weight); };
      void setCompName(const std::string CompName) {_compName.push_back(CompName); };

      // These are Optionals 
      void setRadLength(double RadLength) {_radLength=RadLength; };
      void setIntLength(double IntLength) {_intLength=IntLength; };
      void setRefIndex(double RefIndex) {_refIndex=RefIndex; };
      void setTemperature(double Temperature) {_matTemperature=Temperature; };
      void setPressure(double Pressure) {_matPressure=Pressure; };
      void setState(const std::string State) {_matState=State; };
      void setTcut(double Tcut) {_matTcut=Tcut; };


    private:

      // Basic data members (define a Material)

      std::string _matName;      // Material name
      double _matDensity;        // Material density
      double _matZeff;           // Material effective Z (to allow building material from scratch)
      double _matAeff;           // Material effective A (to allow building material from scratch)
      int _nbrComp;              // number of components < 0 : if number of atomes
      //                      > 0 : if fractions
      std::vector<int> _iflg;             // iflg = 0 : Material made of Elements
      //      = 1 : Material made of other Materials
      std::vector<double> _weight;        // weight of each component
      std::vector<std::string> _compName; // name of each component

      // These are Optionals 
      double _radLength;            // Radiation length of the Material 
      double _intLength;            // Interaction length of the Material
      double _refIndex;             // refraction index of the Material
      double _matTemperature;       // Material Temperature
      double _matPressure;          // Material Pressure
      std::string _matState;        // Material State (Gas, liquid, ...)
      double _matTcut;              // Maximum energy transfer allowed per interaction step (this means that for the dE/dx is used the restricted energy loss rate parameterization instead of the Bethe-Bloch)

      friend bool testCdb(const MatMaterialObj*, const MatMaterialObj*);
  };
}
#endif
