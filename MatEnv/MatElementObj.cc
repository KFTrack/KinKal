//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatElementObj.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatElementObj (transient version)
//      Source file
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
#include <iostream>
#include "MatEnv/MatElementObj.hh"
using std::cout;
using std::endl;
namespace MatEnv {

  //-------------------------------
  // Collaborating Class Headers --
  //-------------------------------
  // Constructor to create Elements from a combination of Isotopes

  MatElementObj::MatElementObj() 
    : _elmName(" "),
    _elmSymbol(" "),
    _elmZeff(0),
    _elmAeff(0),
    _nbrIsotope(0),
    _weight(0),
    _isotopeName(0)
  {
  }

  MatElementObj::MatElementObj(const std::string& name, const std::string& symbol, 
      int Zeff, double Aeff, int nbrIsotope, 
      std::vector<double>* weight, 
      std::vector<std::string>* isotopeName)
  {
    _elmName = name;
    _elmSymbol = symbol;
    _elmZeff = Zeff;
    _elmAeff = Aeff;
    _nbrIsotope = nbrIsotope;
    for (int i=0; i<nbrIsotope; i++) {
      _weight.push_back((*weight)[i]);
      _isotopeName.push_back((*isotopeName)[i]);
    }
  }

  MatElementObj::~MatElementObj() {;}

  //------------
  // Operator --
  //------------
  MatElementObj::MatElementObj(const MatElementObj& elmcp)
  {
    _elmName = elmcp.getName();
    _elmSymbol = elmcp.getSymbol();
    _elmZeff = elmcp.getZeff();
    _elmAeff = elmcp.getAeff();
    _nbrIsotope = elmcp.getNbrIsotope();
    for(int index=0; index<_nbrIsotope; index++)
    {
      _weight.push_back(elmcp.getWeight(index));   
      _isotopeName.push_back(elmcp.getIsotopeName(index));
    }
  }

  MatElementObj& MatElementObj::operator= (const MatElementObj& elmt)
  {
    if (&elmt == this) return *this;

    _elmName = elmt.getName();
    _elmSymbol = elmt.getSymbol();
    _elmZeff = elmt.getZeff();
    _elmAeff = elmt.getAeff();
    _nbrIsotope = elmt.getNbrIsotope();
    for(int index=0; index<_nbrIsotope; index++)
    {
      _weight.push_back(elmt.getWeight(index));   
      _isotopeName.push_back(elmt.getIsotopeName(index));
    }

    return *this;
  }

  int MatElementObj::operator==(const MatElementObj& other) const
  {
    int equal = 1;
    if (_elmName != other.getName() || _elmSymbol != other.getSymbol() || 
	_elmZeff != other.getZeff() || _elmAeff != other.getAeff() || 
	_nbrIsotope || other.getNbrIsotope()) equal = 0;
    for(int index=0; index<_nbrIsotope; index++)
    {       
      if (_weight[index] != other.getWeight(index) ||  
	  _isotopeName[index] != other.getIsotopeName(index)) equal = 0;
    }
    return equal;
  }

  void MatElementObj::print()
  {
    cout << "Name: " << getName() << "  Symbol: " << getSymbol() 
      << "  Z: " << getZeff() << "  Aeff: " << getAeff() 
      << "  Nisotope: " << getNbrIsotope() << endl; 
    for(int idx=0; idx<getNbrIsotope(); idx++) {
      cout << "IsoName: " << getIsotopeName(idx) << "  Weight: " 
	<< getWeight(idx) << endl;
    } 
  }
}





