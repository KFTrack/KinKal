//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: RecoMatFactory.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class RecoMatFactory 
//      Source file (see.hh file for more details)
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   October 21, 1998 - Talby : created
//-----------------------------------------------------------------------------


//----------------------
// This Class Headers --
//----------------------
#include "MatEnv/RecoMatFactory.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "MatEnv/ErrLog.hh" 
#include "MatEnv/ElmPropObj.hh"
#include "MatEnv/MtrPropObj.hh" 
#include "MatEnv/MatElmDictionary.hh"
#include "MatEnv/MatMtrDictionary.hh"
#include "MatEnv/BbrCollectionUtils.hh"
#include <cstdlib>
#include <string>
#include <map>
namespace MatEnv {

  // Constructors 

  RecoMatFactory::RecoMatFactory()
    : _theElmDict( new MatElmDictionary ),
    _theMtrDict( new MatMtrDictionary ),
    _theElmPropDict( new std::map< std::string*, ElmPropObj*, PtrLess >),
    _theMtrPropDict( new std::map< std::string*, MtrPropObj*, PtrLess >)
  {
  }

  RecoMatFactory*
    RecoMatFactory::getInstance()
    {
      static RecoMatFactory theFactory;

      return &theFactory;
    }

  ElmPropObj*
    RecoMatFactory::GetElmProperties( const std::string& name )
    {    
      std::map< std::string*, ElmPropObj*, PtrLess >::iterator elmPos;
      if ((elmPos = _theElmPropDict->find((std::string*)&name)) != _theElmPropDict->end()) {
	//    cout << " the ElmPropObj " << name << " is already built ! " << endl;
	return (elmPos->second);
      } else {
	ElmPropObj* theElmProp=0;
	MatElmDictionary::iterator matElmPos;
	if ( (matElmPos = _theElmDict->find((std::string*)&name)) != _theElmDict->end() ) {
	  MatElementObj* theElement=0;
	  theElement = matElmPos->second;
	  theElmProp = new ElmPropObj(theElement);

	  // Store infos in the dictionary
	  std::string* elmKey = new std::string(name);
	  (*_theElmPropDict)[elmKey] = theElmProp;
	} else {
	  ErrMsg( warning ) << "RecoMatFactory - the element: " << name 
	    << "does not exist in CondDB" << endmsg; 
	}
	return theElmProp;
      } 
    }

  MtrPropObj*
    RecoMatFactory::GetMtrProperties(const std::string& name) 
    {    
      std::map< std::string*, MtrPropObj*, PtrLess >::iterator mtrPos;
      if ((mtrPos = _theMtrPropDict->find((std::string*)&name)) != _theMtrPropDict->end()) {
	//    cout << " the MtrPropObj " << name << " is already built ! " << endl;
	return (mtrPos->second);
      } else {
	MtrPropObj* theMtrProp(0);
	MatMtrDictionary::iterator matMtrPos;
	if ( (matMtrPos = _theMtrDict->find((std::string*)&name)) != _theMtrDict->end() ) {
	  MatMaterialObj* theMaterial;
	  theMaterial = matMtrPos->second;
	  int ncomp = theMaterial->getNbrComp();
	  theMtrProp = new MtrPropObj(theMaterial);
	  int iflg;
	  int nAtomes;
	  double fraction;
	  std::string cmpName;
	  for (int i=0; i<abs(ncomp); i++) {
	    iflg = theMaterial->getIflg(i);
	    cmpName = theMaterial->getCompName(i);
	    if(iflg == 0) {
	      if(ncomp > 0) {
		fraction = theMaterial->getWeight(i);
		theMtrProp->AddElement(GetElmProperties(cmpName),fraction);
	      } else {
		nAtomes = int(theMaterial->getWeight(i));
		theMtrProp->AddElement(GetElmProperties(cmpName),nAtomes);
	      }
	    } else if(iflg == 1) {
	      fraction = theMaterial->getWeight(i);
	      theMtrProp->AddMaterial(GetMtrProperties(cmpName),fraction);
	    }	
	  }

	  // Store infos in the dictionary
	  std::string* mtrKey = new std::string(name);
	  (*_theMtrPropDict)[mtrKey] = theMtrProp;
	} else {
	  ErrMsg ( warning )
	    << "RecoMatFactory - the material: " << name 
	    << "does not exist in MatMtrDictionary " << endmsg; 
	}
	return theMtrProp;
      }
    } 

  //--------------
  // Destructor --
  //--------------
  RecoMatFactory::~RecoMatFactory()
  {
    std::map< std::string*, ElmPropObj*, PtrLess >::iterator 
      iterElm = _theElmPropDict->begin();
    for (; iterElm != _theElmPropDict->end(); ++iterElm) {
      delete iterElm->first;
      delete iterElm->second;
    }
    _theElmPropDict->clear();
    delete _theElmPropDict;

    std::map< std::string*, MtrPropObj*, PtrLess >::iterator 
      iterMtr = _theMtrPropDict->begin();
    for (; iterMtr != _theMtrPropDict->end(); ++iterMtr) {
      delete iterMtr->first;
      delete iterMtr->second;
    }
    _theMtrPropDict->clear();
    delete _theMtrPropDict;

    delete _theElmDict;
    delete _theMtrDict;
  }

}





