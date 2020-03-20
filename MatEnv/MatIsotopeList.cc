//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatIsotopeList.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatIsotopeList (transient version)
//      Source file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   June 18, 1998 - Talby : created
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
#include "MatEnv/MatIsotopeList.hh"
#include "MatEnv/ErrLog.hh"
using std::fstream;
using std::ifstream;
namespace MatEnv {

  //-------------------------------
  // Collaborating Class Headers --
  //-------------------------------


  // Constructor to create an Isotope

  MatIsotopeList::MatIsotopeList()
    : _vector(0)
  {
  }

  MatIsotopeList::MatIsotopeList(const std::vector<MatIsotopeObj*>& vector)
    : _vector(vector)
  {
  }

  MatIsotopeList::MatIsotopeList(const std::string& isotopesFile) 
  {
    // open input file isotopesFile to read isotopes one by one 
    ifstream isotopes( isotopesFile.c_str() );
    assert( isotopes.good() );
    if (isotopes.eof()) {
      ErrMsg(fatal) << "MatEnv/MatIsotopesList.data file empty!" << endmsg; 
    }

    std::string chline;
    std::string tagname;
    std::string fline;

    //  Read, skipping comments
    tagname = "Isotopes_list";
    bool tag = false;
    while(!tag){
      do {  
	getline(isotopes, fline);
      } while (fline == "" && !isotopes.eof());
      tag = ( fline.find(tagname) != std::string::npos );
    }
    assert(tag);
    // read first line (Name, Z, N, A)
    for( size_t i=0; i<4; i++) 
    { 
      isotopes >> chline;
    }
    // read Isotopes data (IsotopeName, Z, N, A)
    std::string name;
    int z = 0;
    int n = 0;
    double a=0.;
    isotopes >> name;
    while( !isotopes.eof() )
    {
      isotopes >> z >> n >> a;
      MatIsotopeObj* isoObj = new MatIsotopeObj();
      isoObj->setName(name);
      isoObj->setZ(z);
      isoObj->setN(n);
      isoObj->setA(a);
      _vector.push_back(isoObj);
      //       _vector.push_back(new MatIsotopeObj(name, z, n, a));
      isotopes >> name;
    } 
  }

  MatIsotopeList::~MatIsotopeList() 
  {
    std::for_each(_vector.begin(), _vector.end(), DeleteObject());
    _vector.clear();   
  }
}





