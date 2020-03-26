//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: RecoMatFactory.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class RecoMatFactory 
//      This class is a singleton and is designed in such a way that once the
//      static method getInstance() is invoked the first time, all isotopes,
//      elements and materials are loaded from the Condition DB in memory
//      through the MatEnv package. All transient objects are stored into
//      hash dictionaries and their keys are represented by their names.
//      The public methods GetElmProperties() and GetMtrProperties()
//      store the ElmPropObj and MtrPropObj objects (respectively) 
//      in a dictionary (if not already existing) and return it back 
//      to the caller.
//      
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   October 21, 1998 - Talby : created
//-----------------------------------------------------------------------------

#ifndef RECOMATFACTORY_HH
#define RECOMATFACTORY_HH

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "MatEnv/MatMtrDictionary.hh"
#include "MatEnv/MatElmDictionary.hh"


#include <string>
#include <map>

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
namespace MatEnv {
  class ElmPropObj;
  class MtrPropObj;
  class MatMtrObj;
  class MatMaterialObj;
  class MatElementObj;

  class RecoMatFactory {
    public:

      static RecoMatFactory* getInstance();

      // Destructor
      virtual ~RecoMatFactory();

      const std::map< std::string*, ElmPropObj*, PtrLess>* GetElmPropDict() const
      {return _theElmPropDict;}
      const std::map< std::string*, MtrPropObj*, PtrLess >* GetMtrPropDict() const
      {return _theMtrPropDict;}

      ElmPropObj* GetElmProperties( const std::string& );
      MtrPropObj* GetMtrProperties( const std::string& );

      std::map<std::string*, MatMaterialObj*, PtrLess>* materialDictionary() const 
      { return _theMtrDict; }
      std::map<std::string*, MatElementObj*, PtrLess>* elementDictionary() const 
      { return _theElmDict; }

    private:

      // Singleton: constructor private
      RecoMatFactory(); 

      // Data members
      MatElmDictionary* _theElmDict;
      MatMtrDictionary* _theMtrDict;
      std::map< std::string*, ElmPropObj*, PtrLess >* _theElmPropDict;
      std::map< std::string*, MtrPropObj*, PtrLess >* _theMtrPropDict; 
  };
}
#endif // RECOMATFACTORY_HH


