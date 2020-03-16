#ifndef DETMTRPROP_HH
#define DETMTRPROP_HH 1

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include <string>
#include <vector>

class MatMaterialObj;
class MtrPropObj;
class ElmPropObj;

#include <sys/types.h>

class DetMtrProp
{

public:   

// Destructor
  virtual ~DetMtrProp(){;};


// Selectors

  virtual const std::string& getName() const = 0;
  virtual double getDensity() const = 0;
  virtual double getCdensity() const = 0;
  virtual double getMdensity() const = 0;
  virtual double getAdensity() const = 0;
  virtual double getX0density() const = 0;
  virtual double getX1density() const = 0;
  virtual double getTaul() const = 0;
  virtual double getDEdxFactor() const = 0;
  virtual double getMeanExciEnergy() const = 0;
  virtual const std::vector< double >& getShellCorrectionVector() const = 0;
  virtual double getEnergyTcut() const = 0;
  virtual const std::vector< double >& getVecTau0() const = 0;
  virtual const std::vector< double >& getVecAlow() const = 0;
  virtual const std::vector< double >& getVecBlow() const = 0;
  virtual const std::vector< double >& getVecClow() const = 0;
  virtual const std::vector< double >& getVecZ() const = 0;

  virtual double getRadLength() const = 0;
  virtual double getIntLength() const = 0;

  virtual std::vector< ElmPropObj* >* getElementVector()  const = 0;
  virtual size_t getNumberOfElements() const = 0;
  virtual const std::vector< double >& getFractionVector()        const = 0;
  virtual const std::vector< int >& getAtomsVector()              const = 0;
  virtual const std::vector< double >& getVecNbOfAtomsPerVolume() const = 0;

  virtual double getZ() const = 0;           
  virtual double getA() const = 0;

  virtual const std::string& getState() const = 0;

};

#endif
