//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatElementList.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatElementList (transient version)
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

#ifndef ELEMENTLIST_HH
#define ELEMENTLIST_HH

//----------------------
// Base Class Headers --
//----------------------
#include <string>
#include <vector>

#include "MatEnv/MatElementObj.hh"
//-------------------------------
// Collaborating Class Headers --
//-------------------------------

namespace MatEnv {
  class MatElementList
  {

    public:

      // Constructor to create an Element

      MatElementList();
      MatElementList(const std::string& ElementsFile); 
      MatElementList(const std::vector<MatElementObj*>&);

      // Destructor
      virtual ~MatElementList();

      inline std::vector<MatElementObj*>* getElementVector();
      inline const std::vector<MatElementObj*>* getVector() const;

    private:

      std::vector<MatElementObj*> _vector;

      friend bool testCdb(const MatElementList*, const MatElementList*);

  };

  inline std::vector<MatElementObj*>* MatElementList::getElementVector()
  {
    return &_vector;
  }

  inline const std::vector<MatElementObj*>* MatElementList::getVector() const
  {
    return &_vector;
  }
}
#endif
