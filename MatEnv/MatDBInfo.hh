//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: MatDBInfo.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//	Class MatDBInfo.  Implementation of MaterialInfo interface
//      using the database.
//
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

#ifndef MATDBINFO_HH
#define MATDBINFO_HH

#include "MatEnv/MaterialInfo.hh"
#include "MatEnv/RecoMatFactory.hh"
#include "MatEnv/MtrPropObj.hh"
#include "MatEnv/ErrLog.hh"
#include <string>
#include <map>

namespace MatEnv {

  class DetMaterial;
  class RecoMatFactory;
  class MatBuildEnv;

  class MatDBInfo : public MaterialInfo {
    public:
      MatDBInfo();
      virtual ~MatDBInfo();
      //  Find the material, given the name
      virtual const DetMaterial* findDetMaterial( const std::string& matName ) const;
      template <class T> const T* findDetMaterial( const std::string& matName ) const;
      // utility functions
    private:
      template <class T> T* createMaterial( const std::string& dbName,
	  const std::string& detMatName ) const;
      void declareMaterial( const std::string& dbName, 
	  const std::string& detMatName );
      // Cache of RecoMatFactory pointer
      RecoMatFactory* _genMatFactory;
      // Cache of list of materials for DetectorModel
      std::map< std::string*, DetMaterial*, PtrLess > _matList;
      // Map for reco- and DB material names
      std::map< std::string, std::string > _matNameMap; 
      // function to cast-off const
      MatDBInfo* that() const {
	return const_cast<MatDBInfo*>(this);
      }
      // allow MatBuildEnv to mess with me
      friend class MatBuildEnv;
      friend class MatBuildCoreEnv;
  };

  template <class T> T*
    MatDBInfo::createMaterial( const std::string& db_name,
	const std::string& detMatName ) const
    {
      MtrPropObj* genMtrProp;
      T* theMat;

      if (_genMatFactory == 0)
	that()->_genMatFactory = RecoMatFactory::getInstance();

      genMtrProp = _genMatFactory->GetMtrProperties(db_name);
      if(genMtrProp != 0){
	theMat = new T( detMatName.c_str(), genMtrProp ) ;
	that()->_matList[new std::string( detMatName )] = theMat;
	return theMat;
      } else {
	return 0;
      }
    }

  template <class T> const T*
    MatDBInfo::findDetMaterial( const std::string& matName ) const
    {
      if (_genMatFactory == 0)
	that()->_genMatFactory = RecoMatFactory::getInstance();

      T* theMat;
      std::map< std::string*, DetMaterial*, PtrLess >::const_iterator pos;
      if ((pos = _matList.find((std::string*)&matName)) != _matList.end()) {
	theMat = (T*) pos->second;
      } else {
	// first, look for aliases
	std::string theName;
	std::map< std::string, std::string >::const_iterator matNamePos;
	if ((matNamePos = _matNameMap.find(matName)) != _matNameMap.end()) {
	  theName = matNamePos->second;
	  theMat = createMaterial<T>( theName, matName);
	} else {
	  //then , try to find the material name directly
	  theMat = createMaterial<T>( matName, matName);
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
#endif
