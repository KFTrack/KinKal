//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatMtrDictionary.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatMtrDictionary (transient version)
//      Source file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//      MatMaterialObj* theMaterial;

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
#include <map>

//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/MatMtrDictionary.hh"
#include "MatEnv/ErrLog.hh"

using std::fstream;
namespace MatEnv {

  // Create Constructor

  MatMtrDictionary::MatMtrDictionary(FileFinderInterface const& fileFinder) : fileFinder_(fileFinder)
  {
    std::string fullPath = fileFinder_.matMtrDictionaryFileName();
    MatMaterialList* mtrList = new MatMaterialList(fullPath);
    FillMtrDict(mtrList);
  }

  void MatMtrDictionary::FillMtrDict(MatMaterialList* mtrList)
  {
    std::vector<MatMaterialObj*>* mtrVec = mtrList->getMaterialVector();
    size_t nmaterial = mtrVec->size();
    for (size_t im=0; im<nmaterial; im++){
      //
      // copy the object into the dictionary. The disctionary now has
      // ownership of the copied objects. Use clearAndDestroy in dtor to remove
      // copies from memory.
      MatMaterialObj* Obj = new MatMaterialObj(*(*mtrVec)[im]);
      std::string* key = new std::string(Obj->getName());
      (*this)[key] = Obj;
      //ErrMsg(routine) << "MatMtrDictionary: Inserted Material " << *key << endmsg;
    }
  }

  MatMtrDictionary::~MatMtrDictionary()
  {
    std::map<std::string*, MatMaterialObj*, PtrLess>::iterator
      iter = begin();
    for (; iter != end(); ++iter) {
      delete iter->first;
      delete iter->second;
    }
    clear();
  }
}




