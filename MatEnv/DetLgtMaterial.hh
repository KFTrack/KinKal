// File and Version Information:
//      $Id:  Exp $
//
//  Description:
//  This class defines the standard values and formulas for particles interacting
//  with very light matter as gas.  It is used in the DetectorModel description heirarchy.
//  
//
// Copyright Information:
//
//
//  Authors:
//------------------------------------------------------------------------------
#ifndef DETECTORLIGHTMATERIAL_HH
#define DETECTORLIGHTMATERIAL_HH

//
//  Babar includes
//
#include "MatEnv/DetMaterial.hh"

class DetLgtMaterial : public DetMaterial {
  private:
    double _cutOffEnergybyMean;// variable cut as cutOffEnergybyMean x MeandEdX on step

  public:
    //  Constructor
    // new style
    DetLgtMaterial(const char* detName, const DetMtrProp* detMtrProp);

    virtual ~DetLgtMaterial();

    virtual double energyLoss(double mom,double pathlen,double mass) const;

    virtual double energyGain(double mom,double pathlen, double mass) const;

    virtual double energyDeposit(double mom, double pathlen,double mass) const;

    // raw dE/dx function, used by the above
    virtual double dEdx(double mom,DetMaterial::dedxtype type,
        TrkParticle const& tpart) const;
    double dEdx(double mom,DetMaterial::dedxtype type,double mass,double cutE) const;

    virtual void printAll(std::ostream& os) const;

    double cutOffEnergybyMean() const { return _cutOffEnergybyMean;}
    void setCutOffEnergybyMean(double cutOffEnergybyMean) {_cutOffEnergybyMean = cutOffEnergybyMean;}

};
#endif

