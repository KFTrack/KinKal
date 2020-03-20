//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: MaterialInfo.hh 510 2010-01-14 15:28:09Z stroili $
//
// Description:
//	Class MaterialInfo.  Interface for accessing material information
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

#ifndef MATERIALINFO_HH
#define MATERIALINFO_HH

#include <string>
#include <vector>
namespace MatEnv {
  class DetMaterial;

  class MaterialInfo {
    public:
      MaterialInfo(){;}
      virtual ~MaterialInfo(){;}
      //  Find the material, given the name
      virtual const DetMaterial* findDetMaterial( const std::string& matName ) const = 0;
      const std::vector< std::string >& materialNames() const {
	return _matNameList; }
      std::vector< std::string >& materialNames() {
	return _matNameList; }
    private:
      std::vector< std::string > _matNameList;
  };
}
#endif
