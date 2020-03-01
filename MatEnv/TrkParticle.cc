// particle description used for tracking
// the implementation depends on the mu2e conditions service, but
// not the interface
#include "MatEnv/TrkParticle.hh"
#include <cmath>


TrkParticle::TrkParticle(TrkParticle::Type ptype) : type_(ptype), mass_(mass(ptype))
{ m2_ = mass_*mass_;}

double TrkParticle::mass(TrkParticle::Type type) {
  // avoid calling the Particle Data table on each call by buffering results in statics
  double retval(-1.0);
  switch (type) {
    case e_minus: case e_plus: return e_mass_; 
    case mu_minus: case mu_plus: return mu_mass_;
    case pi_minus: case pi_plus: return pi_mass_;
    case K_minus: case K_plus: return K_mass_;
    case anti_p_minus: case p_plus: return p_mass_;
    default: {}
  }
  return retval;
}

int TrkParticle::charge() const {
  // avoid calling the Particle Data table on each call by buffering results in statics
  int retval(0);
  switch (type_) {
    case e_minus: case mu_minus: case pi_minus: case K_minus: case anti_p_minus: return -1;
    case e_plus: case mu_plus: case pi_plus: case K_plus: case p_plus: return 1;
    default: {}
  }
  return retval;
}


std::string const&
TrkParticle::name() const {
  // I can't use HepPDT for the names as those include Latex characters and so conflict
  // with many standard software lexicons (like root).
  switch (type_) {
    case e_minus: {
		    static std::string e_minus_name = "eMinus";
		    return e_minus_name;
		  }
    case e_plus: {
		   static std::string e_plus_name = "ePlus";
		   return e_plus_name;
		 }
    case mu_minus: {
		     static std::string mu_minus_name = "muMinus";
		     return mu_minus_name;
		   }
    case mu_plus: {
		    static std::string mu_plus_name = "muPlus";
		    return mu_plus_name;
		  }
    case pi_minus: {
		     static std::string pi_minus_name = "piMinus";
		     return pi_minus_name;
		   }
    case pi_plus: {
		    static std::string pi_plus_name = "piPlus";
		    return pi_plus_name;
		  }
    case K_minus: {
		    static std::string K_minus_name = "KMinus";
		    return K_minus_name;
		  }
    case K_plus: {
		   static std::string K_plus_name = "KPlus";
		   return K_plus_name;
		 }
    case anti_p_minus: {
			 static std::string anti_p_minus_name = "pMinus";
			 return anti_p_minus_name;
		       }
    case p_plus: {
		   static std::string p_plus_name = "pPlus";
		   return p_plus_name;
		 }
    default: {
	       static std::string unknown_name("unknown");
	       return unknown_name;
	     }
  }
}

double
TrkParticle::energy(double momentum) const {
  return sqrt(momentum*momentum + m2_);
}

double
TrkParticle::beta(double momentum) const {
  return fabs(momentum)/energy(momentum);
}

double
TrkParticle::betagamma(double momentum) const {
  return fabs(momentum)/mass_;
}

double
TrkParticle::gamma(double momentum) const {
  return energy(momentum)/mass_;
}


