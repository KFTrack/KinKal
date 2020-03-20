//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MtrPropObj.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MtrPropObj 
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
//   3 Nov 2003 - M. Kelsey -- move ideal gas info to ErrMsg(trace)
//-----------------------------------------------------------------------------


//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/MtrPropObj.hh"

#include <iostream>
#include <math.h>
#include <cfloat>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "MatEnv/ErrLog.hh"
#include "MatEnv/ElmPropObj.hh"
#include "MatEnv/MatMaterialObj.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include <string>
#include <vector>
using std::endl;
using namespace CLHEP;
namespace MatEnv {

  static const double GasThreshold    = 10.*mg/cm3;

  // Constructor to create Material from scratch 

  MtrPropObj::MtrPropObj()
    : _matDensity(0),
    _cdensity(0),
    _mdensity(0),
    _adensity(0),
    _x0density(0),
    _x1density(0),
    _taul(0),
    _radLength(0),
    _intLength(0),
    _dEdxFactor(0),
    _meanExciEnergy(0),
    _energyTcut(0),
    _zeff( 0 ),
    _aeff( 0 ),
    _temp( 0 ),
    _pressure( 0 ),
    _maxNbComponents( 0 ),
    _numberOfComponents( 0 ),
    _numberOfElements( 0 ),
    _theElementVector( 0 ),
    _massFractionVector( 0 ),
    _atomsVector( 0 ),
    _theTau0Vector( 0 ),
    _theAlowVector( 0 ),
    _theBlowVector( 0 ),
    _theClowVector( 0 ),
    _theZVector( 0 ),
    _vecNbOfAtomsPerVolume( 0 ),
    _totNbOfAtomsPerVolume( 0 ),
    _totNbOfElectPerVolume( 0 )
  {
    _matName = new std::string;
    _shellCorrectionVector = new std::vector< double >( numShellV );
    _state = new std::string;
  }

  // Constructor to create Material from a combination of Elements or Materials

  MtrPropObj::MtrPropObj( MatMaterialObj* theMaterial )
    : _matDensity( theMaterial->getDensity() ),
    _cdensity(0),
    _mdensity(0),
    _adensity(0),
    _x0density(0),
    _x1density(0),
    _taul(0),
    _radLength( theMaterial->getRadLength() ),
    _intLength( theMaterial->getIntLength() ),
    _energyTcut( theMaterial->getTcut() ),
    _zeff( theMaterial->getZeff() ),
    _aeff( theMaterial->getAeff() ),
    _temp( theMaterial->getTemperature()+STP_Temperature ),
    _pressure( theMaterial->getPressure()*atmosphere ),
    _maxNbComponents( 0 ),
    _numberOfComponents( 0 ),
    _numberOfElements( 0 ),
    _theElementVector( 0 ),
    _massFractionVector( 0 ),
    _atomsVector( 0 ),
    _theTau0Vector( 0 ),
    _theAlowVector( 0 ),
    _theBlowVector( 0 ),
    _theClowVector( 0 ),
    _theZVector( 0 ),
    _vecNbOfAtomsPerVolume( 0 ),
    _totNbOfAtomsPerVolume( 0 ),
    _totNbOfElectPerVolume( 0 )
  {
    _matName = new std::string( theMaterial->getName() );
    _shellCorrectionVector = new std::vector< double >( numShellV );
    _state = new std::string( theMaterial->getState() );

    if (_matDensity*g/cm3 < universe_mean_density) { 
      ErrMsg(error)
	<< " Warning: define a material with density=0 is not allowed. \n"
	<< " The material " << *_matName << " will be constructed with the"
	<< " default minimal density: " << universe_mean_density/(g/cm3) 
	<< "g/cm3" << endmsg;
      _matDensity = universe_mean_density/(g/cm3);
    } 

    int ncomp = theMaterial->getNbrComp();
    if (ncomp == 0) 
    {
      // Initialize theElementVector allocating one
      // element corresponding to this material
      _maxNbComponents = _numberOfComponents = _numberOfElements = 1;
      _theElementVector    = new std::vector< ElmPropObj* >( 1 );
      (*_theElementVector)[0] = new ElmPropObj( *_matName, _zeff, _aeff );
      _massFractionVector = new std::vector< double >( 1 );
      (*_massFractionVector)[0] = 1. ;

      ComputeDerivedQuantities();

    } else {
      _maxNbComponents = abs(ncomp);
      _numberOfComponents = _numberOfElements = 0;
      _theElementVector = new std::vector< ElmPropObj* >(_maxNbComponents);
    }

    if ( *_state == "undefined")
    {
      if (_matDensity*g/cm3 > GasThreshold) *_state = "solid";
      else                                  *_state = "gas";
    }

  }

  MtrPropObj::~MtrPropObj() 
  {
    delete _matName;
    delete _state;

    delete _theElementVector;
    delete _massFractionVector;
    delete _shellCorrectionVector;
    delete _atomsVector;
    delete _vecNbOfAtomsPerVolume;
    delete _theTau0Vector;
    delete _theAlowVector;
    delete _theBlowVector;
    delete _theClowVector;
    delete _theZVector;
  }

  //------------
  // Operator --
  //------------
  MtrPropObj::MtrPropObj(const MtrPropObj& matcp)
  {
    _matName = new std::string( matcp.getName() );
    _matDensity = matcp.getDensity();
    _cdensity = matcp.getCdensity();
    _mdensity = matcp.getMdensity();
    _adensity = matcp.getAdensity();
    _x0density = matcp.getX0density();
    _x1density = matcp.getX1density();
    _taul = matcp.getTaul();
    _dEdxFactor = matcp.getDEdxFactor();
    _meanExciEnergy = matcp.getMeanExciEnergy();
    _shellCorrectionVector
      = new std::vector< double >( matcp.getShellCorrectionVector() );
    _radLength = matcp.getRadLength();
    _intLength = matcp.getIntLength();
    _energyTcut = matcp.getEnergyTcut();

  }

  MtrPropObj&
    MtrPropObj::operator= (const MtrPropObj& matrl)
    {
      if (&matrl == this) return *this;

      *_matName = matrl.getName();
      _matDensity = matrl.getDensity();
      _cdensity = matrl.getCdensity();
      _mdensity = matrl.getMdensity();
      _adensity = matrl.getAdensity();
      _x0density = matrl.getX0density();
      _x1density = matrl.getX1density();
      _taul = matrl.getTaul();
      _dEdxFactor = matrl.getDEdxFactor();
      _meanExciEnergy = matrl.getMeanExciEnergy();
      *_shellCorrectionVector = matrl.getShellCorrectionVector();
      _radLength = matrl.getRadLength();
      _intLength = matrl.getIntLength();
      _energyTcut = matrl.getEnergyTcut();
      *_state = matrl.getState();

      return *this;
    }

  bool
    MtrPropObj::operator==(const MtrPropObj& other) const
    {
      bool equal = true;
      if ( *_matName != other.getName() || 
	  _matDensity != other.getDensity() || 
	  _cdensity != other.getCdensity() ||
	  _mdensity != other.getMdensity() ||
	  _adensity != other.getAdensity() ||
	  _x0density != other.getX0density() ||
	  _x1density != other.getX1density() ||
	  _taul != other.getTaul() ||
	  _dEdxFactor != other.getDEdxFactor() ||
	  _meanExciEnergy != other.getMeanExciEnergy() ||
	  _radLength != other.getRadLength() ||
	  _intLength != other.getIntLength() 
	 ) equal=false;

      const std::vector< double >& myV = *_shellCorrectionVector;
      const std::vector< double >& otherV = other.getShellCorrectionVector();
      size_t i=0, length = myV.size();
      if ( length != otherV.size() )
      {
	equal = false;
      }
      else
      {
	for ( i=0; i<length; i++ ) if ( myV[i] != otherV[i] ) equal = false;
      }

      return equal;
    }

  const std::string&
    MtrPropObj::getName() const
    {
      return *_matName;
    }

  const std::vector< double >&
    MtrPropObj::getShellCorrectionVector() const
    {
      return *_shellCorrectionVector;
    }

  const std::vector< double >&
    MtrPropObj::getFractionVector() const
    {
      return *_massFractionVector;
    }

  const std::vector< int >&
    MtrPropObj::getAtomsVector() const
    {
      return *_atomsVector;
    }

  const std::vector< double >&
    MtrPropObj::getVecNbOfAtomsPerVolume() const
    {
      return *_vecNbOfAtomsPerVolume;
    }

  const std::vector< double >&
    MtrPropObj::getVecTau0() const
    {
      return *_theTau0Vector;
    }

  const std::vector< double >&
    MtrPropObj::getVecAlow() const
    {
      return *_theAlowVector;
    }

  const std::vector< double >&
    MtrPropObj::getVecBlow() const
    {
      return *_theBlowVector;
    }

  const std::vector< double >&
    MtrPropObj::getVecClow() const
    {
      return *_theClowVector;
    }

  const std::vector< double >&
    MtrPropObj::getVecZ()    const
    {
      return *_theZVector;
    }

  const std::string&
    MtrPropObj::getState() const
    {
      return *_state;
    }

  // AddElement -- composition by atom count

  void
    MtrPropObj::AddElement(ElmPropObj* element, int nAtoms)
    {   
      // initialization
      if ( _numberOfElements == 0 ) {
	_atomsVector        = new std::vector< int >( _maxNbComponents );
	_massFractionVector = new std::vector< double >( _maxNbComponents );
      }

      // filling ...
      if ( _numberOfElements < _maxNbComponents ) {
	(*_theElementVector)[_numberOfElements] = element;
	(*_atomsVector)     [_numberOfElements] = nAtoms;
	_numberOfComponents = ++_numberOfElements;
      }
      else
	ErrMsg(error)
	  << "Attempt to add more than the declared number of elements." 
	  << endmsg; 

      // filled.
      if ( _numberOfElements == _maxNbComponents ) {     
	// compute proportion by mass
	size_t i=0;
	double Zmol(0.), Amol(0.);
	for (i=0;i<_numberOfElements;i++) {
	  Zmol +=  (*_atomsVector)[i]*(*_theElementVector)[i]->getZ();
	  Amol +=  (*_atomsVector)[i]*(*_theElementVector)[i]->getA();
	}
	for (i=0;i<_numberOfElements;i++) {
	  (*_massFractionVector)[i] = (*_atomsVector)[i]*(*_theElementVector)[i]->getA()/Amol;
	}
	ComputeDerivedQuantities();
      }
    }

  // AddElement -- composition by fraction of mass

  void
    MtrPropObj::AddElement(ElmPropObj* element, double fraction)
    {
      // if atomsVector is non-NULL, complain. Apples and oranges.  $$$
      if ( _atomsVector != 0 ) {
	ErrMsg(error)
	  << "This material is already being defined via elements by"
	  << "atoms." << endmsg;
      }

      // initialization
      if ( _numberOfComponents == 0)
      {
	_massFractionVector = new std::vector< double >( 100 );
	// zero all the elements before doing +=
	for ( size_t zeroer = 0; zeroer<100; zeroer++ )
	  (*_massFractionVector)[zeroer]=0;
      }


      // filling ...
      if (_numberOfComponents < _maxNbComponents) {
	size_t el = 0;
	while ((el<_numberOfElements)&&(element!=(*_theElementVector)[el])) el++;
	if (el<_numberOfElements) (*_massFractionVector)[el] += fraction;
	else {
	  if(el>=_theElementVector->size()) _theElementVector->resize(el+1);
	  (*_theElementVector)[el] = element;
	  (*_massFractionVector)[el] = fraction;
	  _numberOfElements ++;
	}
	_numberOfComponents++;  
      }
      else
	ErrMsg(error)
	  << "Attempt to add more than the declared number of components."
	  << endmsg;

      // filled.
      if (_numberOfComponents == _maxNbComponents) {
	// check sum of weights -- OK?
	size_t i;
	double wtSum(0.0);
	for (i=0;i<_numberOfElements;i++) { wtSum +=  (*_massFractionVector)[i]; }
	if (fabs(1.-wtSum) > 0.001) {
	  ErrMsg(error)
	    << "WARNING !! - Fractional masses do not sum to 1: Delta is > 0.001"
	    << "( the weights are NOT renormalized; the results may be wrong)" 
	    << endmsg;
	}

	ComputeDerivedQuantities();

      }
    }

  // AddMaterial -- composition by fraction of mass

  void
    MtrPropObj::AddMaterial(MtrPropObj* material, double fraction)
    {
      // if atomsVector is non-NULL, complain. Apples and oranges.  $$$
      if ( _atomsVector != 0 ) {
	ErrMsg(error)
	  << "This material is already being defined via elements by"
	  << "atoms." << endmsg;
      }
      // initialization
      if (_numberOfComponents == 0) {
	_massFractionVector = new std::vector< double >( 100 );
	// zero all the elements before doing +=
	for ( size_t zeroer = 0; zeroer<100; zeroer++ )
	  (*_massFractionVector)[zeroer]=0;
      }


      // filling ...
      if (_numberOfComponents < _maxNbComponents) {
	for (size_t elm=0; elm < material->getNumberOfElements(); elm++)
	{ 
	  ElmPropObj* element = (*(material->getElementVector()))[elm];
	  size_t el = 0;
	  while ((el<_numberOfElements)&&(element!=(*_theElementVector)[el])) el++;
	  if (el<_numberOfElements) (*_massFractionVector)[el] += fraction
	    *(material->getFractionVector())[elm];
	  else {
	    if(el>=_theElementVector->size()) _theElementVector->resize(el+1);
	    (*_theElementVector)[el] = element;
	    (*_massFractionVector)[el] = fraction*(material->getFractionVector())[elm];
	    _numberOfElements ++;
	  }
	} 
	_numberOfComponents++;  
      }
      else
	ErrMsg(error)
	  << "Attempt to add more than the declared number of components."
	  << endmsg;
      // filled.
      if (_numberOfComponents == _maxNbComponents) {
	// check sum of weights -- OK?
	size_t i;
	double wtSum(0.0);
	for (i=0;i<_numberOfElements;i++) { wtSum +=  (*_massFractionVector)[i]; }
	if (fabs(1.-wtSum) > 0.001) {
	  ErrMsg(error)
	    << "WARNING !! - Fractional masses do not sum to 1: Delta is > 0.001"
	    << "( the weights are NOT renormalized; the results may be wrong)" 
	    << endmsg;
	}

	ComputeDerivedQuantities();
      }
    }

  void
    MtrPropObj::ComputeDerivedQuantities()
    {
      // Header routine to compute various properties of material.
      // Atoms density vector, Electrons density, effective Z and A for the material
      double Zi, Ai;
      double density = _matDensity*g/cm3;
      _theTau0Vector = new std::vector< double >( _numberOfElements );
      _theAlowVector = new std::vector< double >( _numberOfElements );
      _theBlowVector = new std::vector< double >( _numberOfElements );
      _theClowVector = new std::vector< double >( _numberOfElements );
      _theZVector = new std::vector< double >( _numberOfElements );
      _totNbOfAtomsPerVolume = 0.;
      _vecNbOfAtomsPerVolume = new std::vector< double >( _numberOfElements );
      _totNbOfElectPerVolume = 0.;
      for (size_t i=0;i<_numberOfElements;i++) 
      {
	Zi = (*_theElementVector)[i]->getZ();
	Ai = (*_theElementVector)[i]->getA();
	Ai *= g/mole;

	ElmPropObj* element = (*_theElementVector)[i];
	(*_theTau0Vector)[i] = element->getTau0(); 
	(*_theAlowVector)[i] = element->getAlow(); 
	(*_theBlowVector)[i] = element->getBlow(); 
	(*_theClowVector)[i] = element->getClow(); 
	(*_theZVector)[i]    = element->getZ(); 

	(*_vecNbOfAtomsPerVolume)[i] = Avogadro*density*(*_massFractionVector)[i]/Ai;
	_totNbOfAtomsPerVolume += (*_vecNbOfAtomsPerVolume)[i];
	_totNbOfElectPerVolume += (*_vecNbOfAtomsPerVolume)[i]*Zi;

      }

      ComputeOtherParams();
      if(_radLength < 0) { ComputeRadiationLength(); }
      ComputeIonisationParam();
      if(_intLength < 0) { ComputeInteractionLength(); }
    }

  void MtrPropObj::ComputeRadiationLength()
  {
    double radinv = 0.0 ;
    for (size_t i=0;i<_numberOfElements;i++) {
      radinv += (*_vecNbOfAtomsPerVolume)[i]*((*_theElementVector)[i]->getRadTsai()); 
    }
    _radLength = (radinv <= 0.0 ? DBL_MAX : 1./radinv);

    // convert _radLength from mm to g/cm2 (_matDensity is in g/cm3)
    _radLength *= _matDensity/10.;

  }

  void MtrPropObj::ComputeInteractionLength()
  {
    double A;
    double sigmaSum = 0.0;
    for (size_t i=0;i<_numberOfElements;i++) {
      A = (*_theElementVector)[i]->getA();
      sigmaSum += (*_vecNbOfAtomsPerVolume)[i] * (0.0245264 * pow(A, 0.7098081));
    }
    sigmaSum /= 6.022045e23;
    _intLength = (sigmaSum > 0.0 ? 1./sigmaSum : DBL_MAX);
    _intLength *= _matDensity/1000.;
  }

  void
    MtrPropObj::ComputeOtherParams()
    {

      // _dEdxFactor (in Mev/(g/cm2)) is exactly equal to 0.153536*Z/A 
      // (see R.M. Sternheimer et al. in Atomic and Nuclear Data Tables 
      // vol. 30, N02, March 1984, page 267.)

      _dEdxFactor = twopi_mc2_rcl2*_totNbOfElectPerVolume;

      // Convert _dEdxFactor from MeV/mm --> MeV/(g/cm2)

      _dEdxFactor = 10*_dEdxFactor/_matDensity;

    }

  void
    MtrPropObj::ComputeIonisationParam()
    {

      // Compute mean excitation energy and shell correction vector
      //

      _taul = (*_theElementVector)[0]->getTaul();
      double logMeanExciEnergy = 0.;

      for (size_t i=0; i<_numberOfElements; i++)
      {
	logMeanExciEnergy += (*_vecNbOfAtomsPerVolume)[i]
	  *((*_theElementVector)[i]->getZ())
	  *log((*_theElementVector)[i]->getMeanExciEnergy());
      }
      logMeanExciEnergy /= _totNbOfElectPerVolume;
      _meanExciEnergy = exp(logMeanExciEnergy);

      for (int j=0; j<=2; j++)
      {
	(*_shellCorrectionVector)[j] = 0.;
	for (size_t k=0; k<_numberOfElements; k++)
	{
	  (*_shellCorrectionVector)[j] += (*_vecNbOfAtomsPerVolume)[k] 
	    *((*_theElementVector)[k]->getShellCorrectionVector()[j]);
	}
	(*_shellCorrectionVector)[j] /= _totNbOfElectPerVolume;
      } 

      // Compute parameters for the density effect correction in DE/Dx formula.
      // The parametrization is from R.M. Sternheimer, Phys. Rev.B,3:3681 (1971)
      //
      const double Cd2 = 4*pi*hbarc_squared*classic_electr_radius;
      const double twoln10 = 2.*log(10.);
      int icase;
      _cdensity = 1. + log(_meanExciEnergy*_meanExciEnergy/(Cd2*_totNbOfElectPerVolume));

      //
      // condensed materials  
      //
      if ((*_state == "solid")||(*_state == "liquid")) {
	const double E100keV  = 100.*keV; 
	const double ClimiS[] = {3.681 , 5.215 };
	const double X0valS[] = {1.0   , 1.5   };
	const double X1valS[] = {2.0   , 3.0   };

	if(_meanExciEnergy < E100keV) icase = 0 ;
	else                          icase = 1 ;

	if(_cdensity < ClimiS[icase]) _x0density = 0.2;
	else                          _x0density = 0.326*_cdensity-X0valS[icase];

	_x1density = X1valS[icase] ; 
	_mdensity = 3.0;

	//special: Hydrogen
	if ((_numberOfElements==1)&&(getZ()==1)) {
	  _x0density = 0.425; _x1density = 2.0; _mdensity = 5.949;
	}
      }

      //
      // gases
      //
      if (*_state == "gas") { 

	const double ClimiG[] = { 10. , 10.5 , 11. , 11.5 , 12.25 , 13.804};
	const double X0valG[] = { 1.6 , 1.7 ,  1.8 ,  1.9 , 2.0   ,  2.0 };
	const double X1valG[] = { 4.0 , 4.0 ,  4.0 ,  4.0 , 4.0   ,  5.0 };

	icase = 5;
	_x0density = 0.326*_cdensity-2.5 ; _x1density = 5.0 ; _mdensity = 3. ; 
	while((icase > 0)&&(_cdensity < ClimiG[icase])) icase-- ;
	_x0density = X0valG[icase]  ; _x1density = X1valG[icase] ;

	//special: Hydrogen
	if ((_numberOfElements==1)&&(getZ()==1.)) {
	  _x0density = 1.837; _x1density = 3.0; _mdensity = 4.754;
	}

	//special: Helium
	if ((_numberOfElements==1)&&(getZ()==2.)) {
	  _x0density = 2.191; _x1density = 3.0; _mdensity = 3.297;
	}

	// change parameters if the gas is not in STP.
	// For the correction the density(STP) is needed. 
	// Density(STP) is calculated here : 

	double DensitySTP = _matDensity*STP_Pressure*_temp/(_pressure*STP_Temperature);

	double ParCorr = log(_matDensity/DensitySTP) ;

	_cdensity  -= ParCorr;
	_x0density -= ParCorr/twoln10 ;
	_x1density -= ParCorr/twoln10 ;
      }

      double Xa = _cdensity/twoln10 ;
      _adensity = twoln10*(Xa-_x0density)/pow((_x1density-_x0density),_mdensity);
    }

  void
    MtrPropObj::print()
    {
      ErrMsg( routine ) << " Name: " << getName()
	<< "  Density: " << getDensity() << endl
	<< " _cdensity: "  << getCdensity()  
	<< " _mdensity: "  << getMdensity()  
	<< " _adensity: "  << getAdensity()  
	<< " _x0density: " << getX0density()  
	<< " _x1density: " << getX1density() << endl
	<< " _taul: "      << getTaul()
	<< " _meanExciEnergy: " << getMeanExciEnergy()
	<< " _radLength: " << getRadLength() 
	<< " _intLength: " << getIntLength()
	<< endmsg;
    }

  double
    MtrPropObj::getZ() const
    { 
      if (_numberOfElements > 1) {
	//    ErrMsg(error)
	//  << "WARNING in getZ. The material: " << *_matName << " is a mixture." 
	//  <<" the Atomic number is not well defined." << endmsg; 
	double Zsum = 0.0;
	for (size_t i = 0; i<_numberOfElements; i++)
	  Zsum += (*_vecNbOfAtomsPerVolume)[i]/_totNbOfAtomsPerVolume * 
	    ((*_theElementVector)[i]->getZ());
	return Zsum;

      } 
      return (*_theElementVector)[0]->getZ();      
    }

  double
    MtrPropObj::getA() const
    { 
      if (_numberOfElements > 1) { 
	// ErrMsg(error)
	//  << "WARNING in getA. The material: " << *_matName << " is a mixture." 
	//  <<" the Atomic mass is not well defined." << endmsg; 
	double Asum = 0.0;
	for (size_t i = 0; i<_numberOfElements; i++)
	  Asum += (*_vecNbOfAtomsPerVolume)[i]/_totNbOfAtomsPerVolume * 
	    ((*_theElementVector)[i]->getA());
	return Asum;
      } 
      return  (*_theElementVector)[0]->getA();      
    }
}

