//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatElmDictionary.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatElmDictionary (Element Dictionary)
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

#ifndef MATELMDICTIONARY_HH
#define MATELMDICTIONARY_HH

//----------------------
// Base Class Headers --
//----------------------
#include <string>
#include <map>

#include "MatEnv/MatElementObj.hh"
#include "MatEnv/MatElementList.hh"
#include "MatEnv/FileFinderInterface.hh"
//-------------------------------
// Collaborating Class Headers --
//-------------------------------
namespace MatEnv {

  class MatElmDictionary : public std::map<std::string*, MatElementObj*, PtrLess>
  {

    public:

      // Constructor 
      MatElmDictionary(FileFinderInterface const& interface =SimpleFileFinder());

      // Destructor
      virtual ~MatElmDictionary();

      void FillElmDict(MatElementList* elmlist);
    private:
      FileFinderInterface const& fileFinder_;
  };
}
#endif /* MATELMDICTIONARY_HH */

