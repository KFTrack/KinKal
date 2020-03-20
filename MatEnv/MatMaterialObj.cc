//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatMaterialObj.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatMaterialObj (Transient version) 
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
#include <stdlib.h>
#include "MatEnv/MatMaterialObj.hh"
using std::cout;
using std::endl;
namespace MatEnv {

  //-------------------------------
  // Collaborating Class Headers --
  //-------------------------------

  // Constructor to create Material from a combination of Elements or Materials

  MatMaterialObj::MatMaterialObj()
    : _matName(" "),
    _matDensity(0),
    _matZeff(0),
    _matAeff(0),
    _nbrComp(0),
    _iflg(0),
    _weight(0),
    _compName(0),
    _radLength(0),
    _intLength(0),
    _refIndex(0),
    _matTemperature(0),
    _matPressure(0),
    _matState(" "),
    _matTcut(0)
  {
  }

  MatMaterialObj::MatMaterialObj(const std::string& name, double density, 
      double Zeff, double Aeff, int nbrComp, 
      std::vector<int>* iflg, 
      std::vector<double>* weight, 
      std::vector<std::string>* compName, 
      double radLength, double intLength, double refIndex,
      double temperature, double pressure, 
      const std::string& state, double Tcut)
  {
    _matName = name;
    _matDensity = density;
    _matZeff = Zeff;
    _matAeff = Aeff;
    _nbrComp = nbrComp;
    _radLength = radLength;
    _intLength = intLength;
    _refIndex = refIndex;
    _matTemperature = temperature;
    _matPressure = pressure;
    _matState = state;
    _matTcut = Tcut;

    for (int i=0; i<abs(nbrComp); i++) {
      _iflg.push_back((*iflg)[i]); 
      _weight.push_back((*weight)[i]);
      _compName.push_back((*compName)[i]);
    }

  }

  MatMaterialObj::~MatMaterialObj() {;}

  //------------
  // Operator --
  //------------
  MatMaterialObj::MatMaterialObj(const MatMaterialObj& matcp)
  {
    _matName = matcp.getName();
    _matDensity = matcp.getDensity();
    _matZeff = matcp.getZeff(); 
    _matAeff = matcp.getAeff(); 
    _nbrComp = matcp.getNbrComp();
    for (int index=0; index<abs(_nbrComp); index++)
    {
      _iflg.push_back(matcp.getIflg(index));   
      _weight.push_back(matcp.getWeight(index));   
      _compName.push_back(matcp.getCompName(index));
    }    
    _radLength = matcp.getRadLength();
    _intLength = matcp.getIntLength();
    _refIndex = matcp.getRefIndex();
    _matTemperature = matcp.getTemperature();
    _matPressure = matcp.getPressure();
    _matState = matcp.getState();
    _matTcut = matcp.getTcut();
  }

  MatMaterialObj& MatMaterialObj::operator= (const MatMaterialObj& matrl)
  {
    if (&matrl == this) return *this;

    _matName = matrl.getName();
    _matDensity = matrl.getDensity();
    _matZeff = matrl.getZeff(); 
    _matAeff = matrl.getAeff(); 
    _nbrComp = matrl.getNbrComp();
    for (int index=0; index<abs(_nbrComp); index++)
    {
      _iflg.push_back(matrl.getIflg(index));   
      _weight.push_back(matrl.getWeight(index));   
      _compName.push_back(matrl.getCompName(index));
    }    
    _radLength = matrl.getRadLength();
    _intLength = matrl.getIntLength();
    _refIndex = matrl.getRefIndex();
    _matTemperature = matrl.getTemperature();
    _matPressure = matrl.getPressure();
    _matState = matrl.getState();            
    _matTcut = matrl.getTcut();

    return *this;
  }

  int MatMaterialObj::operator==(const MatMaterialObj& other) const
  {
    size_t equal=1;
    if (_matName != other.getName() || 
	_matDensity != other.getDensity() || 
	_matZeff != other.getZeff() || 
	_matAeff != other.getAeff() || 
	_nbrComp != other.getNbrComp() ) equal=0; 
    for (int index=0; index<abs(_nbrComp); index++)
    {
      if(_iflg[index] != other.getIflg(index) ||  
	  _weight[index] != other.getWeight(index) ||  
	  _compName[index] != other.getCompName(index) ) equal=0; 
    }
    if(_radLength != other.getRadLength() || 
	_intLength != other.getIntLength() || 
	_refIndex != other.getRefIndex() ||
	_matTemperature != other.getTemperature() || 
	_matPressure != other.getPressure() || 
	_matState != other.getState() ) equal=0;

    return equal;
  }

  void MatMaterialObj::print()
  {
    cout << "Name: " << getName() << "  Density: " << getDensity() 
      << " Zeff: " << getZeff() << " Aeff: " << getAeff() 
      << "  NbrComp: " << getNbrComp() << endl; 
    for (int idx=0; idx<abs(getNbrComp()); idx++) {
      cout << "Iflg: " << getIflg(idx) << "  Weight: " 
	<< getWeight(idx) << "  CompName: " << getCompName(idx)  << endl;    
    } 
    cout << "RadLength: " << getRadLength() << "  IntLength: " << getIntLength() 
      << " RefIndex: " << getRefIndex()  << "  Temperature: " << getTemperature() 
      << "  Pressure: " << getPressure() << "  State: " << getState();
    if (_matTcut>0.0) { cout<< " Tcut: "<<getTcut(); }
    cout << endl;
  }

}




