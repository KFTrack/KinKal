// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetMaterial.cc,v 1.35 2008/09/08 22:04:29 brownd Exp $
//
//  Description:
//  This class defines the standard values and formulas for particles interacting
//  with matter.  It is used in the DetectorModel description heirarchy.
//
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 11/21/96
//           Matthias Steinke 04/09/99
//------------------------------------------------------------------------------
#include "MatEnv/DetMaterial.hh"
#include <iostream>
#include <cfloat>
#include <string>
#include <vector>
using std::endl;
using std::ostream;
//
// setup the static members
//
//  The multiple scattering constant has an effective thickness
//  for the logrithmic term in the scattering formula corresponding to
//  the average material interior to the tracking volume at dip=~30 degrees
//
namespace MatEnv {

  //double DetMaterial::_msmom = 15.0*MeV;
  double DetMaterial::_dgev = 0.153536e2;
  double DetMaterial::_minkappa(1.0e-3);
  //double DetMaterial::_scatterfrac(0.9999); // integrate 99.99% percent of the tail by default, this should be larger
  // if the materials are very thin.
  const double bg2lim = 0.0169;
  const double taulim = 8.4146e-3 ;
  const double twoln10 = 2.0*log(10.);
  const double betapower = 1.667; // most recent PDG gives beta^-5/3 as dE/dx
  const int maxnstep = 10; // maximum number of steps through a single material
  // should be from a physics class
  const double DetMaterial::_alpha(1.0/137.036);

  //
  //  Constructor
  //

  double cm(10.0); // temporary hack
  DetMaterial::DetMaterial(const char* detMatName, const MtrPropObj* detMtrProp):
    _msmom(15.0),
    _scatterfrac(0.9999),
    _cutOffEnergy(1000.),
    _elossType(loss),
    _name(detMatName),
    _za(detMtrProp->getZ()/detMtrProp->getA()),
    _zeff(detMtrProp->getZ()),
    _aeff(detMtrProp->getA()),
    _radthick(detMtrProp->getRadLength()/cm/cm),
    _intLength(detMtrProp->getIntLength()/detMtrProp->getDensity()),
    _meanion(2.*log(detMtrProp->getMeanExciEnergy()*1.0e6)),
    _eexc(detMtrProp->getMeanExciEnergy()),
    _x0(detMtrProp->getX0density()),
    _x1(detMtrProp->getX1density()),
    _delta0(detMtrProp->getDEdxFactor()),
    _afactor(detMtrProp->getAdensity()),
    _mpower(detMtrProp->getMdensity()),
    _bigc(detMtrProp->getCdensity()),
    _density(detMtrProp->getDensity()/cm/cm/cm),
    _noem(detMtrProp->getNumberOfElements()),
    _taul(detMtrProp->getTaul())
  {
    _shellCorrectionVector = 
      new std::vector< double >(detMtrProp->getShellCorrectionVector());
    _vecNbOfAtomsPerVolume = 
      new std::vector< double >(detMtrProp->getVecNbOfAtomsPerVolume());
    _vecTau0 = 
      new std::vector< double >(detMtrProp->getVecTau0());
    _vecAlow = new std::vector< double >(detMtrProp->getVecAlow());
    _vecBlow = new std::vector< double >(detMtrProp->getVecBlow());
    _vecClow = new std::vector< double >(detMtrProp->getVecClow());
    _vecZ = new std::vector< double >(detMtrProp->getVecZ());
    // compute cached values; these are used in detailed scattering models
    _invx0 = _density/_radthick;
    _nbar = _invx0*1.587e7*pow(_zeff,1.0/3.0)/((_zeff+1)*log(287/sqrt(_zeff)));
    _chic2 = 1.57e1*_zeff*(_zeff+1)/_aeff;  
    _chia2_1 = 2.007e-5*pow(_zeff,2.0/3.0);
    _chia2_2 = 3.34*pow(_zeff*_alpha,2);

    if (detMtrProp->getEnergyTcut()>0.0) {
      _cutOffEnergy = detMtrProp->getEnergyTcut();
      _elossType = deposit;
    }
    if (detMtrProp->getState() == "gas" && detMtrProp->getDensity()<0.01) {
      _scatterfrac = 0.999999;
    }
  }

  DetMaterial::~DetMaterial()
  {
    delete _shellCorrectionVector;
    delete _vecNbOfAtomsPerVolume;
    delete _vecTau0;
    delete _vecAlow;
    delete _vecBlow;
    delete _vecClow;
    delete _vecZ;
  }

  //
  //  Multiple scattering function
  //
  double
    DetMaterial::scatterAngleRMS(double mom, double pathlen,double mass) const {
      if(mom>0.0){
	double beta = particleBeta(mom,mass);
	// pdg formulat
	//    double radfrac = fabs(pathlen*_invx0);
	//    double sigpdg = 0.0136*sqrt(radfrac)*(1.0+0.088*log10(radfrac))/(beta*mom);
	// old Kalman formula
	//    double oldsig = 0.011463*sqrt(radfrac)/(mom*particleBeta(mom,mass));
	// DNB 20/1/2011  Updated to use Dahl-Lynch formula from  NIMB58 (1991)
	double invmom2 = 1.0/pow(mom,2);
	double invb2 = 1.0/pow(beta,2);
	// convert to path in gm/cm^2!!!
	double path = fabs(pathlen)*_density;
	double chic2 = _chic2*path*invb2*invmom2;
	double chia2 = _chia2_1*(1.0 + _chia2_2*invb2)*invmom2;
	double omega = chic2/chia2;
	static double vfactor = 0.5/(1-_scatterfrac);
	double v = vfactor*omega;
	static double sig2factor = 1.0/(1+_scatterfrac*_scatterfrac);
	double sig2 = sig2factor*chic2*( (1+v)*log(1+v)/v - 1);
	// protect against underflow
	double sigdl = sqrt(std::max(0.0,sig2));
	// check
	return sigdl;
      } else
	return 1.0; // 'infinite' scattering
    }

  double
    DetMaterial::dEdx(double mom,dedxtype type,double mass) const {
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


	double RateMass = e_mass_/ mass;

	Tmax = 2.*e_mass_*bg2
	  /(1.+2.*gamma*RateMass+RateMass*RateMass) ;

	dedx = log(2.*e_mass_*bg2*Tmax/Eexc2);
	if(type == loss)
	  dedx -= 2.*beta2;
	else {
	  rcut =  ( _cutOffEnergy< Tmax) ? _cutOffEnergy/Tmax : 1;
	  dedx += log(rcut)-(1.+rcut)*beta2;
	}

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
    DetMaterial::energyLoss(double mom, double pathlen,double mass) const {
      // make sure we take positive lengths!
      pathlen = fabs(pathlen);
      double dedx = dEdx(mom,_elossType,mass);
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
	    deltae = step*dEdx(newmom,_elossType,mass);
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
    DetMaterial::energyGain(double mom, double pathlen, double mass) const {
      // make sure we take positive lengths!
      pathlen = fabs(pathlen);
      double dedx = dEdx(mom,_elossType,mass);
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
	  double deltae = -step*dEdx(newmom,_elossType,mass);
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
    DetMaterial::energyDeposit(double mom, double pathlen, double mass) const {
      double dedx = dEdx(mom,deposit,mass);
      return dedx*fabs(pathlen);
    }  

  /********************** end of New Routines **************************/ 

  //
  //  RMS of energy loss.  This is a gaussian approximation, stolen from
  //  Geant3 (see Phys332)
  //
  double
    DetMaterial::energyLossRMS(double mom,double pathlen,double mass) const {
//      double beta = particleBeta(mom,mass);
//      double emax = eloss_emax(mom,mass);
//      double xi = eloss_xi(beta,fabs(pathlen));
//      double kappa = xi/emax;
//      double gam = sqrt(1.0-0.5*pow(beta,2));
//      // formula comes from GFLUCT.F in gphys dnb Jun 4 2004
//      //
//      // this formula seriously overestimates the rms when kappa<0.001
//      // This only really affects electrons
//      // as for heavier particles resolution effects already dominate when we get to
//      // this range.  I'll truncate
//      if(kappa < _minkappa)kappa = _minkappa;
//      double elossrms = xi*sqrt(gam/kappa);
//      //  cout << "beta = " << beta
//      //       << " emax = " << emax
//      //       << " xi = " << xi
//      //       << " kappa = " << kappa
//      //       << " elossrms = " << elossrms << endl;
//      // this value is way too big: scale it down for now but this function needs a rewrite FIXME!
//      return elossrms;
      return 0.5*energyLoss(mom,pathlen,mass);
    }

  //
  //  Functions needed for energy loss calculation, see reference above
  //
  double
    DetMaterial::eloss_emax(double mom,double mass){
      double beta = particleBeta(mom,mass);
      double gamma = particleGamma(mom,mass);
      double mratio = e_mass_/mass;
      double emax = 2*e_mass_*pow(beta*gamma,2)/
	(1+2*gamma*mratio + pow(mratio,2));
      if(mass <= e_mass_)
	emax *= 0.5;
      return emax;
    }

  double
    DetMaterial::eloss_xi(double beta,double pathlen) const{
      return _dgev*_za*_density*fabs(pathlen)/pow(beta,2);
    }

  void
    DetMaterial::print(ostream& os) const {
      os << "Material " << _name << endl;
    }

  void
    DetMaterial::printAll(ostream& os) const {
      os << "Material " << _name << " has properties : " << endl
	<< "  Effective Z = " << _zeff << endl
	<< "  Effective A = " << _aeff << endl
	<< "  Density (g/cm^3) = " << _density*cm*cm*cm  << endl
	<< "  Radiation Length (g/cm^2) = " << _radthick*cm*cm << endl
	<< "  Interaction Length (g/cm^2) = " << _intLength << endl
	//   << "  Mean Ionization energy (MeV) = " << _meanion << endl
	<< "  Mean Ionization energy (MeV) = " << _eexc << endl;
    }

  double
    DetMaterial::maxStepdEdx(double mom,double mass,double dEdx,double tol) {
      // compute betagamma at entrance
      double betagamma = particleBetaGamma(mom,mass);
      double energy = particleEnergy(mom,mass);
      // basic calculation, based on constant dE/dx
      if(dEdx<0.0){
	double maxstep = -tol*energy/dEdx;
	// Modify for steep rise at low momentum
	if(betagamma<2.0)
	  maxstep *= pow(betagamma,2)/betapower;
	return maxstep;
      }
      else
	return 1.0e6; // large step 
    }


  double
    DetMaterial::nSingleScatter(double mom,double pathlen, double mass) const {
      double beta = particleBeta(mom,mass);
      return pathlen*_nbar/pow(beta,2);
    }


  // note the scaterring momentum here is hard-coded to the simplified Highland formula,
  // this should only be used as a reference, not as a real estimate of the scattering sigma!!!!!
  double
    DetMaterial::highlandSigma(double mom,double pathlen, double mass) const {
      if(mom>0.0){
	double radfrac = _invx0*fabs(pathlen);
	return _msmom*sqrt(radfrac)/(mom*particleBeta(mom,mass));
      } else
	return 1.0;
    }
}
