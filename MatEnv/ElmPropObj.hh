//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: ElmPropObj.hh 671 2010-09-20 11:36:17Z stroili $
//
// Description:
//      Class ElmPropObj 
//      this class computes the elements properties (radiation length 
//      and ionisation parameters) from the basic elements data stored 
//      in the condition database. 
//      The functions ComputeDerivedQuantities(), ComputeCoulombFactor(), 
//      ComputeRadTsai(), ComputeIonisationParam() are taken from Geant4 
//      (the 3 first are from G4Element and the last one is from 
//      G4IonisParamElm). 
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby     (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   October 30, 1998 - Talby : created
//-----------------------------------------------------------------------------
#ifndef ELMPROPOBJ_HH
#define ELMPROPOBJ_HH 1

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include <string>
#include <vector>
namespace MatEnv {

  //------------------------------------
  // Collaborating Class Declarations --
  //------------------------------------
  class MatElementObj;
  class ElmPropObj
  {
    public:   
      enum { numInShellV = 3 };

      // Constructor to create Elements from Cond. DB stored info 
      // either directly or from a combination of isotopes
      ElmPropObj();
      ElmPropObj(MatElementObj* theElement);

      // Constructor to create an Element directly from scratch
      // no reference to Cond. DB
      ElmPropObj(const std::string& name, double zeff, double aeff);

      // Destructor
      virtual ~ElmPropObj();

      // print constructor
      void print() const;  

      // Copy constructor
      ElmPropObj(const ElmPropObj&);

      // Operators (Assignment op)
      ElmPropObj& operator= ( const ElmPropObj& );   
      bool operator==( const ElmPropObj& ) const;

      // Selectors
      const std::string& getName() const;
      double getZ() const { return _elmZ; }
      double getA() const { return _elmA; }
      double getRadTsai() const { return _radTsai; }
      double getTau0() const { return _tau0; }
      double getTaul() const { return _taul; }
      double getAlow() const { return _alow; }
      double getBlow() const { return _blow; }
      double getClow() const { return _clow; }
      double getMeanExciEnergy() const { return _meanExciEnergy; }
      const std::vector< double >& getShellCorrectionVector() const;

    private:

      // Compute derived quantities
      void ComputeDerivedQuantities();
      void ComputeCoulombFactor();
      void ComputeRadTsai();
      void ComputeIonisationParam();

      // Basic data members (define an Element)

      std::string* _elmName;
      double _elmZ;
      double _elmA;
      double _radTsai;
      double _tau0;
      double _taul;
      double _alow;
      double _blow;
      double _clow;
      double _meanExciEnergy;
      std::vector< double >* _shellCorrectionVector;
      double _coulomb;
  };
}
#endif





