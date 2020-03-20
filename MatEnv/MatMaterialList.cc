//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatMaterialList.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatMaterialList (transient version)
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
#include <sstream>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#include "MatEnv/BbrCollectionUtils.hh"

//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/MatMaterialList.hh"
#include "MatEnv/ErrLog.hh"
using std::fstream;
using std::ifstream;
namespace MatEnv {

  //-------------------------------
  // Collaborating Class Headers --
  //-------------------------------


  // Constructor to create an Material

  MatMaterialList::MatMaterialList()
    : _vector(0)
  {
  }

  MatMaterialList::MatMaterialList(const std::vector<MatMaterialObj*>& vector)
    : _vector(vector)
  {
  }

  MatMaterialList::MatMaterialList(const std::string& materialsFile) 

  {

    // open input file materialsFile to read materials one by one 
    ifstream materials( materialsFile.c_str() );
    assert( materials.good() );
    if (materials.eof()) {
      ErrMsg(fatal) << "MatEnv/MatMaterialsList.data file empty!" << endmsg; 
    }

    std::string tagname;
    std::string fline;

    //  Read, skipping comments
    tagname = "Materials_list";
    bool tag = false;
    while(!tag){
      do {  
	getline(materials, fline);
      } while (fline == "" && !materials.eof());
      tag = ( fline.find(tagname) != std::string::npos );
    }
    assert(tag);

    // read Materials data (Name, Zeff, Aeff, nComp ...)
    std::string name;
    std::string cpname;
    std::string state;
    std::vector<double> Compweight;
    std::vector<std::string> Compname;
    std::vector<int> Compflag;
    size_t matentry = 0;
    int nbrcomp = 0;
    int iflag = 0;
    double density = 0.;
    double zeff = 0.;
    double aeff = 0.;
    double weight = 0.;
    double radlen = 0.;
    double intlen = 0.;
    double refindex = 0.;
    double temperature = 0.;
    double pressure = 0.;
    double Tcut = 0.;
    materials >> name;
    while( !materials.eof())
    {
      if(name.find("#") != std::string::npos) {
	do {  
	  getline(materials, fline);
	} while (fline == "" && !materials.eof());
	materials >> name;
      } else {
	materials >> density >> zeff >> aeff >> nbrcomp;
	if(nbrcomp != 0) {
	  for (int i=0; i<abs(nbrcomp); i++) 
	  {
	    materials >> weight >> cpname >> iflag;
	    Compflag.push_back(iflag);
	    Compweight.push_back(weight);
	    Compname.push_back(cpname);
	  }

	} else {
	  Compflag.push_back(0);
	  Compweight.push_back(0);
	  Compname.push_back(" ");
	}            
	materials >> radlen >> intlen >> refindex >> temperature 
	  >> pressure >> state;

	getline(materials, fline);
	std::istringstream iss(fline);

	MatMaterialObj* matObj = new MatMaterialObj();
	matObj->setName(name);
	matObj->setDensity(density);
	matObj->setZeff(zeff);
	matObj->setAeff(aeff);
	matObj->setNbrComp(nbrcomp);
	matentry = Compflag.size();
	for (size_t idx=matentry-abs(nbrcomp); idx<matentry; idx++)
	{
	  matObj->setIflg(Compflag[idx]);
	  matObj->setWeight(Compweight[idx]);
	  matObj->setCompName(Compname[idx]); 
	}
	matObj->setRadLength(radlen);
	matObj->setIntLength(intlen);
	matObj->setRefIndex(refindex);
	matObj->setTemperature(temperature);
	matObj->setPressure(pressure);
	matObj->setState(state);
	if (iss>>Tcut) {
	  matObj->setTcut(Tcut);
	}

	_vector.push_back(matObj);

	//      _vector.push_back(new MatMaterialObj(name, density, zeff, 
	//                         aeff, nbrcomp, Compflag, Compweight, 
	//                         Compname, radlen, intlen, refindex, 
	//                         temperature, pressure, state));
	materials >> name;
      }
    }  

  }

  void MatMaterialList::print(std::ostream& out) 
  {

    std::string tagname;

    //  print tag name
    tagname = "Materials_list";
    out<<tagname<<std::endl;

    // write Materials data (Name, Zeff, Aeff, nComp ...)
    std::string name;
    std::string cpname;
    std::string state;
    std::vector<double> Compweight;
    std::vector<std::string> Compname;
    std::vector<int> Compflag;
    int nbrcomp = 0;
    double density = 0.;
    double zeff = 0.;
    double aeff = 0.;
    double radlen = 0.;
    double intlen = 0.;
    double refindex = 0.;
    double temperature = 0.;
    double pressure = 0.;
    double energyCut = 0.;
    std::vector<MatMaterialObj*>::iterator it=_vector.begin();
    for(;it!=_vector.end();it++)
    {
      MatMaterialObj* matObj = (*it);
      name     = matObj->getName();
      density  = matObj->getDensity();
      zeff     = matObj->getZeff();
      aeff     = matObj->getAeff();
      nbrcomp  = matObj->getNbrComp();
      Compflag.resize(abs(nbrcomp));
      Compweight.resize(abs(nbrcomp));
      Compname.resize(abs(nbrcomp));
      for (int idx=0; idx<abs(nbrcomp); idx++)
      {
	Compflag[idx]  =matObj->getIflg(idx);
	Compweight[idx]=matObj->getWeight(idx);
	Compname[idx]  =matObj->getCompName(idx); 
      }

      radlen   = matObj->getRadLength();
      intlen   = matObj->getIntLength();
      refindex = matObj->getRefIndex();
      temperature = matObj->getTemperature();
      pressure = matObj->getPressure();
      state    = matObj->getState();
      energyCut = matObj->getTcut();


      out<< name;
      out<<" "<< density <<" "<< zeff <<" "<< aeff <<" "<< nbrcomp;
      for (int i=0; i<abs(nbrcomp); i++) 
      { 
	out <<" "<< Compweight[i] <<" "<< Compname[i] <<" "<< Compflag[i];
      }

      out <<" "<< radlen <<" "<< intlen <<" "<< refindex <<" "<< temperature 
	<<" "<< pressure <<" "<< state;
      if (energyCut>0.0) { out <<" "<< energyCut; }

      out << std::endl;
    }
  }  



  MatMaterialList::~MatMaterialList() 
  {
    std::for_each(_vector.begin(), _vector.end(), DeleteObject());
    _vector.clear();   
  }
}





