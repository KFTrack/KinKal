//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MtrPropObj.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MtrPropObj 
//      this class computes the materials properties (radiation length 
//      and ionisation parameters) from the basic elements data stored 
//      in the condition database. 
//      The functions AddElement(), AddMaterial(), 
//      ComputeDerivedQuantities(), ComputeRadiationLength(), 
//      ComputeIonisationParam() are taken from Geant4 
//      (the 4 first are from G4Material and the last one is from 
//      G4IonisParamMat). 
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
#ifndef MTRPROPOBJ_HH
#define MTRPROPOBJ_HH 1

//-------------
// C Headers --
//-------------
#include <sys/types.h>

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/DetMtrProp.hh"
#include <string>
#include <vector>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
class ElmPropObj;
class MatMaterialObj;

class MtrPropObj : public DetMtrProp
{
public:   

  enum { numShellV = 3 };

  // Default constructor 
  MtrPropObj(); 

  // Constructor to create Material from a combination of Elements or Materials
  MtrPropObj( MatMaterialObj* theMaterial ); 

  // Copy constructor
  MtrPropObj(const MtrPropObj&);

  // Destructor
  virtual ~MtrPropObj();

  // Operators (Assignment op)
  MtrPropObj& operator=( const MtrPropObj& );   
  bool operator==( const MtrPropObj& ) const;

  // Add and element, giving number of atoms
  void AddElement( ElmPropObj* element, int nAtoms );

  // Add and element or material, giving fraction of mass
  void AddElement ( ElmPropObj* element , double fraction );
  void AddMaterial( MtrPropObj* material, double fraction );

  // print constructor
  void print();

  // Selectors
  virtual const std::string& getName() const;
  virtual double getDensity() const { return _matDensity; }
  virtual double getCdensity() const { return _cdensity; }
  virtual double getMdensity() const { return _mdensity; }
  virtual double getAdensity() const { return _adensity; }
  virtual double getX0density() const { return _x0density; }
  virtual double getX1density() const { return _x1density; }
  virtual double getTaul() const { return _taul; }
  virtual double getDEdxFactor() const { return _dEdxFactor; }
  virtual double getMeanExciEnergy() const { return _meanExciEnergy; }
  virtual const std::vector< double >& getShellCorrectionVector() const;
  virtual double getEnergyTcut() const { return _energyTcut; }

  virtual double getRadLength() const { return _radLength; }
  virtual double getIntLength() const { return _intLength; }

  virtual std::vector< ElmPropObj* >* getElementVector()  const {return _theElementVector; }
  virtual size_t getNumberOfElements()       const {return _numberOfElements; }
  virtual const std::vector< double >& getFractionVector()        const;
  virtual const std::vector< int >& getAtomsVector()              const;
  virtual const std::vector< double >& getVecNbOfAtomsPerVolume() const;
  virtual const std::vector< double >& getVecTau0() const;
  virtual const std::vector< double >& getVecAlow() const;
  virtual const std::vector< double >& getVecBlow() const;
  virtual const std::vector< double >& getVecClow() const;
  virtual const std::vector< double >& getVecZ()    const;

  virtual double getZ() const;           
  virtual double getA() const;

  virtual const std::string& getState() const;

private:

// Compute derived quantities
  void ComputeDerivedQuantities();
  void ComputeRadiationLength();
  void ComputeIonisationParam();
  void ComputeOtherParams();
  void ComputeInteractionLength();

private:

// Basic data members (define a Material)

  std::string* _matName;     
  double _matDensity;     
  double _cdensity;
  double _mdensity;
  double _adensity;
  double _x0density;
  double _x1density;
  double _taul;
  double _radLength;    
  double _intLength;

  double _dEdxFactor;
  double _meanExciEnergy;
  std::vector< double >* _shellCorrectionVector;
  double _energyTcut;

  double _zeff;
  double _aeff;
  std::string* _state;
  double _temp;
  double _pressure;
  size_t _maxNbComponents;
  size_t _numberOfComponents;
  size_t _numberOfElements;

  std::vector< ElmPropObj* >* _theElementVector;
  std::vector< double >* _massFractionVector;
  std::vector< int >* _atomsVector;
  std::vector< double >* _theTau0Vector;
  std::vector< double >* _theAlowVector;
  std::vector< double >* _theBlowVector;
  std::vector< double >* _theClowVector;
  std::vector< double >* _theZVector;

// Derived data members (computed from the basic data members)

  std::vector< double >* _vecNbOfAtomsPerVolume;
                                       // vector of nb of atoms per volume
  double  _totNbOfAtomsPerVolume;      // total nb of atoms per volume 
  double  _totNbOfElectPerVolume;      // total nb of electrons per volume 

};

#endif
