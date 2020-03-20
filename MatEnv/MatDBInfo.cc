//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: MatDBInfo.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//	Class MatDBInfo.  
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Dave Brown                      LBL
//
// Copyright Information:
//	Copyright (C) 1999		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------
#include "MatEnv/MatDBInfo.hh"
#include "MatEnv/DetMaterial.hh"

#include <string>
#include <map>
namespace MatEnv {

  MatDBInfo::MatDBInfo() :
    _genMatFactory(0)
  {;}

  MatDBInfo::~MatDBInfo() {
    // delete the materials
    std::map< std::string*, DetMaterial*, PtrLess >::iterator 
      iter = _matList.begin();
    for (; iter != _matList.end(); ++iter) {
      delete iter->first;
      delete iter->second;
    }
    _matList.clear();
  }

  void 
    MatDBInfo::declareMaterial( const std::string& db_name,
	const std::string& detMatName )
    {
      _matNameMap[detMatName] = db_name;
      materialNames().push_back( detMatName );
      return;
    }

  const DetMaterial*
    MatDBInfo::findDetMaterial( const std::string& matName ) const
    {
      if (_genMatFactory == 0)
	that()->_genMatFactory = RecoMatFactory::getInstance();

      DetMaterial* theMat;
      std::map< std::string*, DetMaterial*, PtrLess >::const_iterator pos;
      if ((pos = _matList.find((std::string*)&matName)) != _matList.end()) {
	theMat = pos->second;
      } else {
	// first, look for aliases
	std::string theName;
	std::map< std::string, std::string >::const_iterator matNamePos;
	if ((matNamePos = _matNameMap.find(matName)) != _matNameMap.end()) {
	  theName = matNamePos->second;
	  theMat = createMaterial<DetMaterial>( theName, matName);
	} else {
	  //then , try to find the material name directly
	  theMat = createMaterial<DetMaterial>( matName, matName);
	  // if we created a new material directly, add it to the list
	  if(theMat != 0)that()->declareMaterial(matName,matName);
	}
      }
      if(theMat == 0){
	ErrMsg( error ) << "MatDBInfo: Cannot find requested material " << matName
	  << "." << endmsg;
      }
      return theMat;
    }
}
