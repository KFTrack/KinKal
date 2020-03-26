//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatElementObj.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatElementObj (transient version) 
//      Header file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby     (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   June 18, 1998 - Talby : created
//-----------------------------------------------------------------------------
#ifndef MATELEMENTOBJ_HH
#define MATELEMENTOBJ_HH 1

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
// Mandatory include file

#include <string>
#include <vector>

namespace MatEnv {
  class MatElementObj
  {

    public:   

      // Constructor to create Elements either directly or from a combination 
      // of isotopes
      MatElementObj();
      MatElementObj(const std::string& name, const std::string& symbol, 
	  int Zeff, double Aeff, int nbrIsotope, 
	  std::vector<double>* weight, 
	  std::vector<std::string>* isotopeName);
      // Destructor
      virtual ~MatElementObj();

      // print constructor
      void print();  

      // Copy constructor
      MatElementObj(const MatElementObj&);

      // Operators (Assignment op)
      MatElementObj& operator= ( const MatElementObj& );   
      int operator==( const MatElementObj& ) const;

      // Selectors
      std::string getName() const { return std::string(_elmName); };         
      std::string getSymbol() const { return std::string(_elmSymbol); };
      int getZeff() const { return _elmZeff; };           
      double getAeff() const { return _elmAeff; };             
      int getNbrIsotope() const { return _nbrIsotope; };
      double getWeight(size_t index) const { return _weight[index]; };
      std::string getIsotopeName(size_t index) const { return std::string(_isotopeName[index]); };

      void setName(const std::string Name) {_elmName=Name; };
      void setSymbol(const std::string Symbol) {_elmSymbol=Symbol; };
      void setZeff(int Zeff) {_elmZeff=Zeff; };
      void setAeff(double Aeff) {_elmAeff=Aeff; };
      void setNbrIsotope(int NbrIsotope) {_nbrIsotope=NbrIsotope; };
      void setWeight(double Weight) {_weight.push_back(Weight); };
      void setIsotopeName(const std::string IsotopeName) {_isotopeName.push_back(IsotopeName); };

    private:   

      // Basic data members (define an Element)

      std::string _elmName;
      std::string _elmSymbol;
      int _elmZeff;
      double _elmAeff;
      int  _nbrIsotope;
      std::vector<double> _weight;   
      std::vector<std::string> _isotopeName;

      friend bool testCdb(const MatElementObj*, const MatElementObj*);
  };
}
#endif





