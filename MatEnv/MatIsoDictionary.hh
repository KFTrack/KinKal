//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatIsoDictionary.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatIsoDictionary (Isotope Dictionary)
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

#ifndef MATISODICTIONARY_HH
#define MATISODICTIONARY_HH

//----------------------
// Base Class Headers --
//----------------------
#include <string>
#include <map>

#include "MatEnv/MatIsotopeObj.hh"
#include "MatEnv/MatIsotopeList.hh"
#include "MatEnv/BbrCollectionUtils.hh"
#include "MatEnv/FileFinderInterface.hh"
//-------------------------------
// Collaborating Class Headers --
//-------------------------------
namespace MatEnv {

  class MatIsoDictionary : public std::map<std::string*, MatIsotopeObj*, PtrLess>
  {

    public:

      // Constructor 
      MatIsoDictionary(FileFinderInterface const& interface =SimpleFileFinder());

      // Destructor
      virtual ~MatIsoDictionary();

      void FillIsoDict(MatIsotopeList* isolist);
    private:
      FileFinderInterface const& fileFinder_;
  };
}
#endif /* MATISODICTIONARY_HH */

