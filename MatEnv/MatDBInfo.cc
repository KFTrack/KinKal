//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: MatDBInfo.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatDBInfo.
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Dave Brown                      LBL
//
// Copyright Information:
//      Copyright (C) 1999              Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------
#include "KinKal/MatEnv/MatDBInfo.hh"

#include <string>
#include <map>
namespace MatEnv {

  MatDBInfo::MatDBInfo(FileFinderInterface const& interface, DetMaterialConfig const& dmconf ) :
    _genMatFactory(RecoMatFactory::getInstance(interface)), _dmconf(dmconf)
  {;}

  MatDBInfo::MatDBInfo(FileFinderInterface const& interface) : MatDBInfo(interface, DetMaterialConfig()){}


  MatDBInfo::~MatDBInfo() {;}

  void
    MatDBInfo::declareMaterial( const std::string& db_name,
        const std::string& detMatName )
    {
      _matNameMap[detMatName] = db_name;
      materialNames().push_back( detMatName );
      return;
    }

  const std::shared_ptr<DetMaterial> MatDBInfo::findDetMaterial( const std::string& matName ) const
  {
    std::shared_ptr<DetMaterial> theMat;
    std::map< std::string, std::shared_ptr<DetMaterial> >::const_iterator pos;
    if ((pos = _matList.find(matName)) != _matList.end()) {
      theMat = pos->second;
    } else {
      // first, look for aliases
      std::string theName;
      std::map< std::string, std::string >::const_iterator matNamePos;
      if ((matNamePos = _matNameMap.find(matName)) != _matNameMap.end()) {
        theName = matNamePos->second;
        theMat = createDetMaterial( theName, matName);
      } else {
        //then , try to find the material name directly
        theMat = createDetMaterial( matName, matName);
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

  std::shared_ptr<DetMaterial> MatDBInfo::createDetMaterial( const std::string& db_name,
      const std::string& detMatName ) const
  {
    MtrPropObj* genMtrProp;
    std::shared_ptr<DetMaterial> theMat;

    genMtrProp = _genMatFactory->GetMtrProperties(db_name);
    if(genMtrProp != 0){
      theMat = std::make_shared<DetMaterial> ( detMatName.c_str(), genMtrProp, _dmconf ) ;
      that()->_matList[ detMatName ] = theMat;
      return theMat;
    } else {
      return nullptr;
    }
  }
}

