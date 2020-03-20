//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatIsotopeObj.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatIsotopeObj (transient version) 
//      Header file
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   June 18, 1998 - Talby : created
//-----------------------------------------------------------------------------
#ifndef MATISOTOPEOBJ_HH
#define MATISOTOPEOBJ_HH 1

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include <string>
namespace MatEnv {

  class MatIsotopeObj
  {
    public:   

      // Constructor to create Isotopes
      MatIsotopeObj();
      MatIsotopeObj(const std::string& name, int z, int n, double a);

      // Destructor
      virtual ~MatIsotopeObj();

      // copy constructor
      MatIsotopeObj(const MatIsotopeObj&);

      // print constructor
      void print();

      // Operators (Assignment op)
      MatIsotopeObj& operator= ( const MatIsotopeObj& );   
      int operator==( const MatIsotopeObj& ) const;


      // Selectors
      const std::string& getName() const { return _isoName; };              
      int getZ() const { return _isoZ; };           
      int getN() const { return _isoN; };           
      double getA() const { return _isoA; };             

      void setName(const std::string Name) {_isoName=Name;};
      void setZ(int Z) {_isoZ=Z;};
      void setN(int N) {_isoN=N;};
      void setA(double A) {_isoA=A;};

    private:

      // Basic data members (define an Isotope)

      std::string _isoName;
      int _isoZ;
      int _isoN;
      double _isoA;

      friend bool testCdb(const MatIsotopeObj*, const MatIsotopeObj*);
  };
}
#endif
