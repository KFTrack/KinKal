//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatIsotopeObj.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatIsotopeObj (transient version)  
//      Source file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   June 18, 1998 - Talby : created
//-----------------------------------------------------------------------------


//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/MatIsotopeObj.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include <iostream>
using std::cout;
using std::endl;
namespace MatEnv {

  // Constructor to create an Isotope

  MatIsotopeObj::MatIsotopeObj()
    : _isoName(""),
    _isoZ(0),
    _isoN(0),
    _isoA(0)
  {
  }

  MatIsotopeObj::MatIsotopeObj(const std::string& name, int z, int n, double a)
  {
    _isoName = name;
    _isoZ = z;
    _isoN = n;
    _isoA = a;

  }

  MatIsotopeObj::~MatIsotopeObj() {}

  //------------
  // Operator --
  //------------
  MatIsotopeObj::MatIsotopeObj(const MatIsotopeObj& isocp)
  {
    _isoName = isocp.getName();
    _isoZ = isocp.getZ();
    _isoN = isocp.getN();
    _isoA = isocp.getA();
  }

  MatIsotopeObj& MatIsotopeObj::operator= (const MatIsotopeObj& istp)
  {
    if (&istp == this) return *this;

    _isoName = istp.getName();
    _isoZ = istp.getZ();
    _isoN = istp.getN();
    _isoA = istp.getA();

    return *this;
  }

  int MatIsotopeObj::operator==(const MatIsotopeObj &other) const 
  {
    if (_isoName == other.getName() && _isoZ == other.getZ() && 
	_isoN == other.getN() && _isoA == other.getA()) return 1;
    return 0;
  }

  void MatIsotopeObj::print()
  {
    cout << "Name: " << getName() << "  Z: " << getZ() 
      << "  N: " << getN() << "  A: " << getA() << endl; 
  }
}

