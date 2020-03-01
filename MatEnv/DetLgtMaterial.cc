// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: Exp $
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
//
//------------------------------------------------------------------------------
#include "MatEnv/DetLgtMaterial.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "MatEnv/TrkParticle.hh"
#include <iostream>
#include <cfloat>
#include <string>
#include <vector>
using std::endl;
using std::ostream;
using namespace CLHEP;

const double bg2lim = 0.0169;
const double taulim = 8.4146e-3 ;
const double twoln10 = 2.0*log(10.);
//const double betapower = 1.667; // most recent PDG gives beta^-5/3 as dE/dx
const int maxnstep = 10; // maximum number of steps through a single material
//const double cutOffEnergy = 5.;//maximum allowed energy loss
                               //(if particle will loose more - than it will be not reconstructable)
//const double cutOffEnergybyMean=5.;// variable cut as cutOffEnergybyMean x MeandEdX on step

//
//  Constructor
//

DetLgtMaterial::DetLgtMaterial(const char* detMatName, const DetMtrProp* detMtrProp):
                DetMaterial(detMatName, detMtrProp),
                _cutOffEnergybyMean(5.)
{
   _msmom = ((13.6-0.5)*MeV);
  _scatterfrac = 0.99999999;
  _cutOffEnergy = 5.;
}

DetLgtMaterial::~DetLgtMaterial()
{}

double DetLgtMaterial::dEdx(double mom,DetMaterial::dedxtype type,
        TrkParticle const& tpart) const {
      return dEdx(mom,type,tpart.mass(),_cutOffEnergy);
}

double
DetLgtMaterial::dEdx(double mom,DetMaterial::dedxtype type,double mass,double cutE) const {
  if(mom>0.0){
    double Eexc2 = _eexc*_eexc ;
  
  // New energy loss implementation
  
    double Tmax,gamma2,beta2,bg2,rcut,delta,x,sh,dedx ;
    double beta  = particleBeta(mom,mass) ;
    double gamma = particleGamma(mom,mass) ;
    double tau = gamma-1. ;
  
  // high energy part , Bethe-Bloch formula 
  
    beta2 = beta*beta ;
    gamma2 = gamma*gamma ;
    bg2 = beta2*gamma2 ;
 
    static TrkParticle epart(TrkParticle::e_minus);
    double RateMass = epart.mass() / mass;
  
    Tmax = 2.*electron_mass_c2*bg2
      /(1.+2.*gamma*RateMass+RateMass*RateMass) ;

    if(cutE<1e-6) cutE=1e-6;//protect from zero path by 1eV
    rcut =  ( cutE< Tmax) ? cutE/Tmax : 1;

    dedx = log(2.*electron_mass_c2*bg2*Tmax/Eexc2);
    /*if(type == loss)
      dedx -= 2.*beta2;
    else*/
    dedx += log(rcut)-(1.+rcut)*beta2;
  
// density correction 
    x = log(bg2)/twoln10 ;
    if ( x < _x0 ) {
      if(_delta0 > 0) {
	delta = _delta0*pow(10.0,2*(x-_x0));
      }
      else {
	delta = 0.;
      }
    } else {
      delta = twoln10*x - _bigc;
      if ( x < _x1 )
	delta += _afactor * pow((_x1 - x), _mpower);
    } 
  
  // shell correction          
    if ( bg2 > bg2lim ) {
      sh = 0. ;      
      x = 1. ;
      for (int k=0; k<=2; k++) {
	x *= bg2 ;
	sh += (*_shellCorrectionVector)[k]/x;
      }
    }
    else {
      sh = 0. ;      
      x = 1. ;
      for (int k=0; k<2; k++) {
	x *= bg2lim ;
	sh += (*_shellCorrectionVector)[k]/x;
      }
      sh *= log(tau/_taul)/log(taulim/_taul);
    }
    dedx -= delta + sh ;
    dedx *= -_dgev*_density*_za / beta2 ;
    return dedx;
  } else
    return 0.0;
}


double 
DetLgtMaterial::energyLoss(double mom, double pathlen,double mass) const {
// make sure we take positive lengths!
  pathlen = fabs(pathlen);
  double dedx = dEdx(mom,loss,mass,_cutOffEnergy);
  double cutTail=std::min(fabs(dedx)*_cutOffEnergybyMean*pathlen,_cutOffEnergy);
  std::cout<<"dedx "<<dedx<<" cutOff0 "<<fabs(dedx)*_cutOffEnergybyMean*pathlen<<" cutTail "<<cutTail<<std::endl;
  //  double beta  = particleBeta(mom,mass) ;
  // double cutTail=_dgev*_density*_za / beta/beta*fabs(pathlen)*exp(0.2);
  dedx=dEdx(mom,loss,mass,cutTail);
  std::cout<<"dedx_1 "<<dedx<<std::endl;
// see how far I can step, within tolerance, given this energy loss
  double maxstep = maxStepdEdx(mom,mass,dedx);
// if this is larger than my path, I'm done
  if(maxstep>pathlen){
    return dedx*pathlen;
  } else {
// subdivide the material
    unsigned nstep = std::min(int(pathlen/maxstep) + 1,maxnstep);
    double step = pathlen/nstep;
    double energy = particleEnergy(mom,mass);
    double deltae = step*dedx;
    double newenergy(energy+deltae);
    double eloss(deltae);
    for(unsigned istep=0;istep<nstep-1;istep++){
      if(newenergy>mass){
// compute the new dedx given the new momentum
        double newmom = particleMomentum(newenergy,mass);
        deltae = step*dEdx(newmom,loss,mass,cutTail);
// compute the loss in this step
        eloss += deltae;
        newenergy += deltae;
      } else {
// lost all kinetic energy; stop
        eloss = mass-energy;
        break;
      }
    }
    return eloss;
  }
}

double 
DetLgtMaterial::energyGain(double mom, double pathlen, double mass) const {
  // make sure we take positive lengths!
  pathlen = fabs(pathlen);
  double dedx = dEdx(mom,loss,mass,_cutOffEnergy);
  double cutTail=std::min(fabs(dedx)*_cutOffEnergybyMean*pathlen,_cutOffEnergy);
  //double beta  = particleBeta(mom,mass) ;
  //double cutTail=_dgev*_density*_za / beta/beta*fabs(pathlen)*exp(0.2);
  dedx=dEdx(mom,loss,mass,cutTail);
// see how far I can step, within tolerance, given this energy loss
  double maxstep = maxStepdEdx(mom,mass,dedx);
// if this is larger than my path, I'm done
  if(maxstep>pathlen){
    return -dedx*pathlen;
  } else {
// subdivide the material
    unsigned nstep = std::min(int(pathlen/maxstep) + 1,maxnstep);
    double step = pathlen/nstep;
    double energy = particleEnergy(mom,mass);
    double deltae = -step*dedx;
// move to the middle of the slice of material
    double newenergy(energy+deltae);
    double egain(deltae);
    for(unsigned istep=0;istep<nstep-1;istep++){
// compute the new dedx given the new momentum
      double newmom = particleMomentum(newenergy,mass);
      double deltae = -step*dEdx(newmom,loss,mass,cutTail);
      egain += deltae;
      newenergy += deltae;
    }
    return egain;
  }
}  
//
// calculate the energy deposited in an absorber. That's similiar to 
// energyLoss, but the delta electron correction in the Bethe Bloch is
// switched on, and there is no Bremsstrahlung 
//
double 
DetLgtMaterial::energyDeposit(double mom, double pathlen, double mass) const {
  double dedx = dEdx(mom,deposit,mass,_cutOffEnergy);
  return dedx*fabs(pathlen);
}

void
DetLgtMaterial::printAll(ostream& os) const {
  os << "Material " << _name << " has properties : " << endl
  << "  Effective Z = " << _zeff << endl
  << "  Effective A = " << _aeff << endl
  << "  Density (g/cm^3) = " << _density*cm*cm*cm << endl
  << "  Radiation Length (g/cm^2) = " << _radthick*cm*cm<< endl
  << "  Mean Ionization energy (MeV) = " << _eexc << endl
     << "  density X0 = " << _x0 << " ,X1 = " <<_x1
     << "  ,delta0 = " << _delta0 <<" ,A = "<<_afactor
     << "  ,Mpower = "<<_mpower << " ,C = "<<_bigc<<endl
     << "  shell taul = " << _taul 
     << " nshells= "<<_shellCorrectionVector->size()<<" (";
  for(unsigned int i=0;i<_shellCorrectionVector->size();i++) 
    os<<(*_shellCorrectionVector)[i]<<(i+1<_shellCorrectionVector->size()?",":")");
  os<<endl;
}
