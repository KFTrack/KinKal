//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: MatDBInfo.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatDBInfo.  Implementation of MaterialInfo interface
//      using the database.
//
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

#ifndef MATDBINFO_HH
#define MATDBINFO_HH

#include "KinKal/MatEnv/MaterialInfo.hh"
#include "KinKal/MatEnv/RecoMatFactory.hh"
#include "KinKal/MatEnv/MtrPropObj.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/MatEnv/ErrLog.hh"
#include "KinKal/MatEnv/FileFinderInterface.hh"
#include <string>
#include <map>
#include <memory>

namespace MatEnv {

  class RecoMatFactory;
  class MatBuildEnv;

  class MatDBInfo : public MaterialInfo {
    public:
      MatDBInfo(FileFinderInterface const& interface, DetMaterial::energylossmode elossmode);
      virtual ~MatDBInfo();
      //  Find the material, given the name
      const std::shared_ptr<DetMaterial> findDetMaterial( const std::string& matName ) const override;
      // utility functions
    private:
      std::shared_ptr<DetMaterial> createDetMaterial( const std::string& dbName,
          const std::string& detMatName ) const;
      void declareMaterial( const std::string& dbName,
          const std::string& detMatName );
      // Cache of RecoMatFactory pointer
      RecoMatFactory* _genMatFactory;
      // Cache of list of materials for DetectorModel
      std::map< std::string, std::shared_ptr<DetMaterial> > _matList;
      // Map for reco- and DB material names
      std::map< std::string, std::string > _matNameMap;
      // function to cast-off const
      MatDBInfo* that() const {
        return const_cast<MatDBInfo*>(this);
      }
      DetMaterial::energylossmode _elossmode;
  };

}
#endif
