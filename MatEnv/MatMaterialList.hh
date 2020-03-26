//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatMaterialList.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatMaterialList (transient version)
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

#ifndef MATMATERIALLIST_HH
#define MATMATERIALLIST_HH

//----------------------
// Base Class Headers --
//----------------------
#include <string>
#include <vector>
#include <iostream>

#include "MatEnv/MatMaterialObj.hh"
//-------------------------------
// Collaborating Class Headers --
//-------------------------------
namespace MatEnv {

  class MatMaterialList
  {

    public:

      // Constructor to create an Material

      MatMaterialList();
      MatMaterialList(const std::string& MaterialsFile); 
      MatMaterialList(const std::vector<MatMaterialObj*>& vector);

      void print(std::ostream& out = std::cout);

      // Destructor
      virtual ~MatMaterialList();

      inline std::vector<MatMaterialObj*>* getMaterialVector();
      inline const std::vector<MatMaterialObj*>* getVector() const;

    private:

      std::vector<MatMaterialObj*> _vector;

      friend bool testCdb(const MatMaterialList*, const MatMaterialList*);

  };

  inline std::vector<MatMaterialObj*>* MatMaterialList::getMaterialVector()
  {
    return &_vector;
  }

  inline const std::vector<MatMaterialObj*>* MatMaterialList::getVector() const
  {
    return &_vector;
  }
}
#endif
