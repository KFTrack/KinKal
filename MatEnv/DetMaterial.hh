// File and Version Information:
//      $Id: DetMaterial.hh,v 1.35 2008/07/21 23:33:11 brownd Exp $
//
//  Description:
//  This class defines the standard values and formulas for particles interacting
//  with matter.  It is used in the DetectorModel description heirarchy.
//
//
// Copyright Information:
//      Copyright (C) 1996      Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 11/21/96
//           Orion Ning, 01/12/21
//------------------------------------------------------------------------------
#ifndef DETECTORMATERIAL_HH
#define DETECTORMATERIAL_HH
//
//  Babar includes
//
#include "KinKal/MatEnv/MtrPropObj.hh"
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <memory>

namespace MatEnv {
  class DetMaterial{
    public:
 //Energy Loss model: choose 'mpv' for the Most Probable Energy Loss, or 'moyalmean' for the mean calculated via the Moyal Distribution approximation, see end of file for more information, as well as discussion about radiative losses
      enum energylossmode {mpv=0, moyalmean};
      //
      //  Constructor
      // new style
      DetMaterial(const char* detName, const MtrPropObj* detMtrProp, energylossmode);

      ~DetMaterial();
      //
      //  Access
      //
      std::string const& name() const {return _name;}

      //
      //  comparison function
      //
      bool operator == (const DetMaterial& other) const {
        return _name == other._name; }

      // total energy loss (ionization and radiation) and variance
      double energyLoss(double mom,double pathlen,double mass) const;
      double energyLossVar(double mom,double pathlen,double mass) const;
      double energyLossRMS(double mom,double pathlen,double mass) const {
        return sqrt(energyLossVar(mom,pathlen,mass));
      }

      // ionization energy loss functions using closed-form Moyal calculations, see end of DetMaterial.cc for more information
      double ionizationEnergyLoss(double mom,double pathlen,double mass) const;
      // most probable value of energy loss
      double ionizationEnergyLossMPV(double mom,double pathlen,double mass) const;
      double ionizationEnergyLossRMS(double mom,double pathlen,double mass) const;
      double ionizationEnergyLossVar(double mom,double pathlen,double mass) const {
        double elrms = ionizationEnergyLossRMS(mom,pathlen,mass);
        return elrms*elrms;
      }
      // radiation (brehmsstrahlung) energy loss calculation. This is relevant only for electrons
      //
      double radiationEnergyLoss(double mom,double pathlen, double mass) const;
      double radiationEnergyLossVar(double mom,double pathlen, double mass) const;
      // Single Gaussian approximation
      double scatterAngleVar(double mom,double pathlen,double mass) const;
      double scatterAngleRMS(double mom,double pathlen,double mass) const {
        return sqrt(scatterAngleVar(mom,pathlen,mass));
      }

      double highlandSigma(double mom,double pathlen, double mass) const;

      static double particleEnergy(double mom,double mass) {
        return sqrt(pow(mom,2)+pow(mass,2)); }

      static double particleMomentum(double energy,double mass) {
        return sqrt(particleMom2(energy,mass)); }

      static double particleMom2(double energy,double mass)  {
        return std::max(pow(energy,2)-pow(mass,2),0.0); }

      static double particleBeta(double mom,double mass) {
        return mom/particleEnergy(mom,mass); }

      static double particleGamma(double mom,double mass)  {
        return particleEnergy(mom,mass)/mass; }

      static double particleKinEnergy(double E, double mass)  {
        return E - mass; }

      static double particleBetaGamma(double mom,double mass) {
        return mom/mass; }

      //
      //  functions used to compute energy loss
      //
      double eloss_xi(double beta,double pathlen) const;
      double densityCorrection(double bg2) const;
      double shellCorrection(double bg2, double tau) const;
      double moyalMean(double deltap, double xi) const;
    protected:
      //
      //  Constants used in material calculations
      //
      static double _dgev; // energy characterizing energy loss
      static const double _alpha; // fine structure constant
      energylossmode _elossmode;
      double _scatterfrac; // fraction of scattering distribution to include in RMS
      double _e_lpm; // brehmsstrahlung energy lower limit
      double _ymax; // bremsstrahlong fractional energy upper limit
      //
      //  Specific data for this material
      //
      std::string _name;
      double _za; // ratio atomic number to atomic weight
      double _zeff; // effective Z of our material
      double _aeff; // effective Z of our material
      double _radthick; // radiation thickness in g/cm**2
      double _intLength; // ineraction length from MatMtrObj in g/cm**2
      double _eexc; // mean ionization energy loss for new e_loss routine
      double _x0; /*  The following specify parameters for energy loss. see
                      Sternheimer etal,'Atomic Data and
                      Nuclear Data Tables', 1984 (40) 267 */
      double _x1;
      double _delta0;
      double _afactor;
      double _mpower;
      double _bigc;
      double _density;

      int _noem;
      std::vector< double >* _shellCorrectionVector;
      double _taul;

      // cached values to speed calculations
      double _invx0;
      double _chic2;
      double _chia2_1;
      double _chia2_2;
      // helper function
      double radiationNormalization(double mom) const;

    public:
      // baseic accessors
      double ZA()const {return _za;}
      double zeff() const { return _zeff;}
      double aeff() const { return _aeff;}
      double radiationLength()const {return _radthick;}
      double intLength()const {return _intLength;}
      double eexc() const { return _eexc; }
      double X0()const {return _x0;}
      double X1()const {return _x1;}
      double delta0()const {return _delta0;}
      double aFactor()const {return _afactor;}
      double mPower()const {return _mpower;}
      double bigC()const {return _bigc;}
      double density()const {return _density;}
      double inverseX0() const { return _invx0; }
      // returns fraction of radiation lengths traversed for a given
      // physical distance through this material
      double radiationFraction(double pathlen) const {
        return _density*pathlen/_radthick; }
      void print(std::ostream& os) const;
      void printAll(std::ostream& os ) const;

      // parameters used in ionization energy loss randomization
      // scattering parameter
      double scatterFraction() const { return _scatterfrac;}
      static constexpr double e_mass_ = 5.10998910E-01; // electron mass in MeVC^2
  };
}
#endif

