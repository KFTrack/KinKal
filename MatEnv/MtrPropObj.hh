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
#ifndef MatEnv_MtrPropObj_hh
#define MatEnv_MtrPropObj_hh

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
#include <string>
#include <vector>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
namespace MatEnv {
  class ElmPropObj;
  class MatMaterialObj;

  class MtrPropObj
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
      ~MtrPropObj();

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
      const std::string& getName() const;
      double getDensity() const { return _matDensity; }
      double getCdensity() const { return _cdensity; }
      double getMdensity() const { return _mdensity; }
      double getAdensity() const { return _adensity; }
      double getX0density() const { return _x0density; }
      double getX1density() const { return _x1density; }
      double getTaul() const { return _taul; }
      double getDEdxFactor() const { return _dEdxFactor; }
      double getMeanExciEnergy() const { return _meanExciEnergy; }
      const std::vector< double >& getShellCorrectionVector() const;
      double getEnergyTcut() const { return _energyTcut; }

      double getRadLength() const { return _radLength; }
      double getIntLength() const { return _intLength; }

      std::vector< ElmPropObj* >* getElementVector()  const {return _theElementVector; }
      size_t getNumberOfElements()       const {return _numberOfElements; }
      const std::vector< double >& getFractionVector()        const;
      const std::vector< int >& getAtomsVector()              const;
      const std::vector< double >& getVecNbOfAtomsPerVolume() const;
      const std::vector< double >& getVecTau0() const;
      const std::vector< double >& getVecAlow() const;
      const std::vector< double >& getVecBlow() const;
      const std::vector< double >& getVecClow() const;
      const std::vector< double >& getVecZ()    const;

      double getZ() const;           
      double getA() const;

      const std::string& getState() const;

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
}
#endif
