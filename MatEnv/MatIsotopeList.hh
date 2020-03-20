//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatIsotopeList.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatIsotopeList (transient version)
//      Header file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   June 18, 1998 - Talby : created
//-----------------------------------------------------------------------------

#ifndef MATISOTOPELIST_HH
#define MATISOTOPELIST_HH

//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/MatIsotopeObj.hh"

#include <string>
#include <vector>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
namespace MatEnv {

  class MatIsotopeList
  {

    public:

      // Constructor to create an Isotope

      MatIsotopeList();
      MatIsotopeList(const std::string& isotopesFile); 
      MatIsotopeList(const std::vector<MatIsotopeObj*>& vector);

      // Destructor
      virtual ~MatIsotopeList();

      inline std::vector<MatIsotopeObj*>* getIsotopeVector();
      inline const std::vector<MatIsotopeObj*>* getVector() const;

    private:

      std::vector<MatIsotopeObj*> _vector;

      friend bool testCdb(const MatIsotopeList*, const MatIsotopeList*);
  };

  inline std::vector<MatIsotopeObj*>* MatIsotopeList::getIsotopeVector()
  {
    return &_vector;
  }

  inline const std::vector<MatIsotopeObj*>* MatIsotopeList::getVector() const
  {
    return &_vector;
  }
}
#endif
