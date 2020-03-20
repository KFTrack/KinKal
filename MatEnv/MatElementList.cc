//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatElementList.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatElementList (transient version)
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
#include <vector>
#include "MatEnv/BbrCollectionUtils.hh"
#include <algorithm>
#include <fstream>
#include <assert.h>

//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/MatElementList.hh"
#include "MatEnv/ErrLog.hh"
using std::fstream;
using std::ifstream;

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

// Constructor to create an Element
namespace MatEnv {

  MatElementList::MatElementList()
    : _vector(0)
  {
  }

  MatElementList::MatElementList(const std::vector<MatElementObj*>& vector)
    : _vector(vector)
  {
  }

  MatElementList::MatElementList(const std::string& elementsFile) 

  {

    // open input file elementsFile to read elements one by one 
    ifstream elements( elementsFile.c_str() );
    assert( elements.good() );
    if (elements.eof()) {
      ErrMsg(fatal) << "MatEnv/MatElementsList.data file empty!" << endmsg; 
    }

    std::string tagname;
    std::string fline;

    //  Read, skipping comments
    tagname = "Elements_list";
    bool tag = false;
    while(!tag){
      do {  
	getline(elements, fline);
      } while (fline == "" && !elements.eof());
      tag = ( fline.find(tagname) != std::string::npos );
    }
    assert(tag);

    // read Elements data (Name, Symbol Z, Aeff, nisotopes, ...)
    std::string name;
    std::string symbol;
    std::string isname;
    std::vector<double> Isoabundance;
    std::vector<std::string> Isoname;
    int isoentry = 0;
    int zeff = 0;
    int nisotopes = 0;
    double aeff = 0.;
    double abundance = 0.;
    elements >> name;
    while( !elements.eof())
    {
      elements >> symbol >> zeff >> aeff >> nisotopes;
      if(nisotopes != 0) {
	for(int i=0; i<nisotopes; i++) 
	{ 
	  elements >> abundance;
	  Isoabundance.push_back(abundance);
	}
	for(int j=0; j<nisotopes; j++) 
	{ 
	  elements >> isname;
	  Isoname.push_back(isname);
	}

      } else {
	Isoabundance.push_back(0);
	Isoname.push_back(" ");
      }            
      MatElementObj* elmObj = new MatElementObj();
      elmObj->setName(name);
      elmObj->setSymbol(symbol);
      elmObj->setZeff(zeff);
      elmObj->setAeff(aeff);
      elmObj->setNbrIsotope(nisotopes);
      isoentry = Isoabundance.size();
      for (int idx=isoentry-nisotopes; idx<isoentry; idx++)
      {
	elmObj->setWeight(Isoabundance[idx]);
	elmObj->setIsotopeName(Isoname[idx]);
      }
      _vector.push_back(elmObj);
      //      _vector.push_back(new MatElementObj(name, symbol, zeff, aeff, Isoabundance, Isoname));
      elements >> name;
    }
  }

  MatElementList::~MatElementList() 
  {
    std::for_each(_vector.begin(), _vector.end(), DeleteObject());
    _vector.clear();   
  }

}




