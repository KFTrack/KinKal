//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: ElmPropObj.cc 671 2010-09-20 11:36:17Z stroili $
//
// Description:
//      Class ElmPropObj 
//      Source file (see.hh file for more details) 
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

//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/ElmPropObj.hh"

#include <iostream>
#include <math.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Units/PhysicalConstants.h"
#include "MatEnv/MatElementObj.hh"
#include <string>
#include <vector>
using namespace CLHEP;
namespace MatEnv {

  // Constructor to create Elements from a combination of Isotopes

  ElmPropObj::ElmPropObj() 
    : _elmName( new std::string ),
    _elmZ(0),
    _elmA(0),
    _radTsai(0),
    _tau0(0),
    _taul(0),
    _alow(0),
    _blow(0),
    _clow(0),
    _meanExciEnergy(0),
    _shellCorrectionVector( new std::vector< double >( numInShellV ) ),
    _coulomb(0)
  {
  }

  ElmPropObj::ElmPropObj( const std::string& name, double zeff, double aeff )
    : _elmName( new std::string( name ) ),
    _elmZ( zeff ),
    _elmA( aeff ),
    _shellCorrectionVector( new std::vector< double >( numInShellV  ) )
  {
    ComputeDerivedQuantities();
  }

  ElmPropObj::ElmPropObj( MatElementObj* theElement )
    : _elmName( new std::string( theElement->getName() ) ),
    _elmZ( theElement->getZeff() ),
    _elmA( theElement->getAeff() ),
    _shellCorrectionVector( new std::vector< double >( numInShellV ) )
  {
    ComputeDerivedQuantities();
  }

  ElmPropObj::~ElmPropObj()
  {
    delete _elmName;
    delete _shellCorrectionVector;
  }

  ElmPropObj::ElmPropObj( const ElmPropObj& elmcp )
    : _elmName( new std::string( elmcp.getName() ) ),
    _elmZ( elmcp.getZ() ),
    _elmA( elmcp.getA() ),
    _radTsai( elmcp.getRadTsai() ),
    _tau0( elmcp.getTau0() ),
    _taul( elmcp.getTaul() ),
    _alow( elmcp.getAlow() ),
    _blow( elmcp.getBlow() ),
    _clow( elmcp.getClow() ),
    _meanExciEnergy( elmcp.getMeanExciEnergy() ),
    _shellCorrectionVector( new std::vector< double >
	( elmcp.getShellCorrectionVector() ) )
  {
  }

  //------------
  // Operator --
  //------------
  ElmPropObj&
    ElmPropObj::operator=( const ElmPropObj& elmt )
    {
      if (&elmt == this) return *this;

      *_elmName = elmt.getName();
      _elmZ = elmt.getZ();
      _elmA = elmt.getA();
      _radTsai = elmt.getRadTsai();
      _tau0 = elmt.getTau0();
      _taul = elmt.getTaul();
      _alow = elmt.getAlow();
      _blow = elmt.getBlow();
      _clow = elmt.getClow();
      _meanExciEnergy = elmt.getMeanExciEnergy();
      *_shellCorrectionVector = elmt.getShellCorrectionVector();

      return *this;
    }

  bool
    ElmPropObj::operator==( const ElmPropObj& other ) const
    {
      bool equal = true;
      if ( *_elmName != other.getName() || _elmZ != other.getZ() || 
	  _elmA != other.getA() || _radTsai != other.getRadTsai() || 
	  _taul != other.getTaul() || _tau0 != other.getTau0() ||
	  _alow != other.getAlow() || _blow != other.getBlow() ||
	  _clow != other.getClow() || 
	  _meanExciEnergy != other.getMeanExciEnergy()) equal = false;
      return equal;
    }

  const std::string&
    ElmPropObj::getName() const
    {
      return *_elmName;
    }

  const std::vector< double >&
    ElmPropObj::getShellCorrectionVector() const
    {
      return *_shellCorrectionVector;
    }

  void
    ElmPropObj::ComputeDerivedQuantities()
    {
      // some basic functions of the atomic number

      // Radiation Length
      ComputeCoulombFactor();
      ComputeRadTsai(); 

      // parameters for energy loss by ionisation   
      ComputeIonisationParam();
    }

  void
    ElmPropObj::ComputeCoulombFactor()
    {
      //  Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)

      const double k1 = 0.0083 , k2 = 0.20206 ,k3 = 0.0020 , k4 = 0.0369 ;

      double az2 = (fine_structure_const*_elmZ)*(fine_structure_const*_elmZ);
      double az4 = az2 * az2;

      _coulomb = (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
    }

  void
    ElmPropObj::ComputeRadTsai()
    {
      //  Compute Tsai's Expression for the Radiation Length
      //  (Phys Rev. D50 3-1 (1994) page 1254)

      const double Lrad_light[]  = {5.31  , 4.79  , 4.74 ,  4.71} ;
      const double Lprad_light[] = {6.144 , 5.621 , 5.805 , 5.924} ;

      const double logZ3 = log((double)_elmZ)/3.;

      double Lrad, Lprad;
      int iz = (int)(_elmZ+0.5) - 1 ;
      if (iz <= 3) { Lrad = Lrad_light[iz] ;  Lprad = Lprad_light[iz] ; }
      else { Lrad = log(184.15) - logZ3 ; Lprad = log(1194.) - 2*logZ3 ; }

      _radTsai = 4*alpha_rcl2*_elmZ*(_elmZ*(Lrad-_coulomb) + Lprad); 
    }

  void ElmPropObj::ComputeIonisationParam()
  {
    double Z3 = pow(_elmZ,1./3.);
    //  double ZZ3 = pow(_elmZ*(_elmZ+1.),1./3.);
    //  double logZ3 = log((double)_elmZ)/3.;

    // parameters for energy loss by ionisation
    _meanExciEnergy = 16.*pow(_elmZ,0.9)*eV;
    _tau0 = 0.1*Z3*MeV/proton_mass_c2;
    _taul = 2.*MeV/proton_mass_c2;
    double rate = _meanExciEnergy/electron_mass_c2;
    double w = _taul*(_taul+2.);
    double BetheBlochLow = (_taul+1.)*(_taul+1.)*log(2.*w/rate)/w - 1.;
    BetheBlochLow = 2.*_elmZ*twopi_mc2_rcl2*BetheBlochLow ; 

    _clow = sqrt(_taul)*BetheBlochLow;
    _alow = 6.458040 *_clow/_tau0;
    double Taum = 0.035*Z3*MeV/proton_mass_c2;
    _blow =-3.229020*_clow/(_tau0*sqrt(Taum));

    // Shell correction factors
    double rate2 = rate*rate;

    (*_shellCorrectionVector)[0] = ( 1.10289e5 + 5.14781e8*rate)*rate2 ;
    (*_shellCorrectionVector)[1] = ( 7.93805e3 - 2.22565e7*rate)*rate2 ;
    (*_shellCorrectionVector)[2] = (-9.92256e1 + 2.10823e5*rate)*rate2 ;
  }

  void
    ElmPropObj::print() const
    {
      std::cout << "Name: " << getName() << "  Z: " << getZ() << " A: " 
	<< getA() 
	<< " RadTsai: " << getRadTsai() << "  Taul: " << getTaul() 
	<< " MeanExciEnergy: " << getMeanExciEnergy() << std::endl;
    }
}
