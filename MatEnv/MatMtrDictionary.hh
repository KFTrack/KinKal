//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatMtrDictionary.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatMtrDictionary (Material Dictionary)
//      Header file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   October 15, 1998 - Talby : created
//-----------------------------------------------------------------------------

#ifndef MATMTRDICTIONARY_HH
#define MATMTRDICTIONARY_HH

//----------------------
// Base Class Headers --
//----------------------
#include <string>
#include <map>

#include "MatEnv/BbrCollectionUtils.hh"
#include "MatEnv/MatMaterialObj.hh"
#include "MatEnv/MatMaterialList.hh"
#include "MatEnv/FileFinderInterface.hh"
//-------------------------------
// Collaborating Class Headers --
//-------------------------------
namespace MatEnv {

  class MatMtrDictionary : public std::map<std::string*, MatMaterialObj*, PtrLess>
  {

    public:

      // Constructor 
      MatMtrDictionary(FileFinderInterface const& interface =SimpleFileFinder());

      // Destructor
      virtual ~MatMtrDictionary();

      void FillMtrDict(MatMaterialList* mtrlist);
    private:
      FileFinderInterface const& fileFinder_;

  };
}
#endif /* MATMTRDICTIONARY_HH */




