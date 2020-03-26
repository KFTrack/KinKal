//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatIsoDictionary.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatIsoDictionary (transient version)
//      Source file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   October 15, 1998 - Talby : created
//-----------------------------------------------------------------------------


//----------------------
// C++ Headers --
//----------------------
#include <fstream>
#include <assert.h>

#include <vector>
#include <algorithm>

//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/MatIsoDictionary.hh"
#include "MatEnv/ErrLog.hh"
using std::fstream;
namespace MatEnv {

  // Create Constructor

  MatIsoDictionary::MatIsoDictionary(FileFinderInterface const& fileFinder ) : fileFinder_(fileFinder)
  {
    std::string fullPath = fileFinder_.matIsoDictionaryFileName();
    MatIsotopeList* mtrList = new MatIsotopeList(fullPath);
    FillIsoDict(mtrList);
  }

  void MatIsoDictionary::FillIsoDict(MatIsotopeList* isoList)
  {
    std::vector<MatIsotopeObj*>* isoVec = isoList->getIsotopeVector();
    size_t nisotope = isoVec->size();
    for (size_t is=0; is<nisotope; is++)
    {
      // copy the object into the dictionary. The disctionary now has
      // ownership of the copied objects. Use clearAndDestroy in dtor to remove
      // copies from memory.
      MatIsotopeObj* Obj = new MatIsotopeObj(*(*isoVec)[is]);
      std::string* key = new std::string(Obj->getName());
      (*this)[key] = Obj;
    }
  }

  MatIsoDictionary::~MatIsoDictionary()
  {
    std::map<std::string*, MatIsotopeObj*, PtrLess>::iterator iter;
    for (iter = begin(); iter != end(); ++iter) {
      delete iter->first;
      delete iter->second;
    }
    clear();
  }
}





