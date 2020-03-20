//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatElmDictionary.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatElmDictionary (transient version)
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

#include <string>
#include <vector>
#include <algorithm>
#include "MatEnv/BbrCollectionUtils.hh"

//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/MatElmDictionary.hh"
#include "MatEnv/ErrLog.hh"
using std::fstream;

namespace MatEnv {
  // Create Constructor

  MatElmDictionary::MatElmDictionary(FileFinderInterface const& fileFinder ) : fileFinder_(fileFinder)
  {
    std::string fullPath = fileFinder_.matElmDictionaryFileName();
    MatElementList* elmList = new MatElementList(fullPath);
    FillElmDict(elmList);
  }

  void MatElmDictionary::FillElmDict(MatElementList* elmList)
  {
    std::vector<MatElementObj*>* elmVec = elmList->getElementVector();
    size_t nelement = elmVec->size();
    for (size_t ie=0; ie<nelement; ie++) {
      //
      // copy the object into the dictionary. The disctionary now has
      // ownership of the copied objects. Use clearAndDestroy in dtor to remove
      // copies from memory.
      MatElementObj* Obj = new MatElementObj(*(*elmVec)[ie]);
      std::string* key = new std::string(Obj->getName());
      (*this)[key] = Obj;
    }
  }
  MatElmDictionary::~MatElmDictionary()
  {
    std::map<std::string*, MatElementObj*, PtrLess>::iterator iter;
    for (iter = begin(); iter != end(); ++iter) {
      delete iter->first;
      delete iter->second;
    }
    clear();
  }
}





