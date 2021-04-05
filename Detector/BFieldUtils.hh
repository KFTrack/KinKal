#ifndef KinKal_BFieldMapUtils_hh
#define KinKal_BFieldMapUtils_hh
//
//  Utility functions for working with BFieldMap and KTRAJ objects.  Used in the
//  Kinematic Kalman filter fit
//
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Detector/BFieldMap.hh"
#include <algorithm>
#include <cmath>
#include <iostream>
namespace KinKal {
  namespace BFieldUtils {
    // speed of light in units to convert Tesla to mm (bending radius)
    static double constexpr cbar() { return CLHEP::c_light/1000.0; }
    
    // integrate the residual magentic force over the given KTRAJ and range, NOT described by the intrinsic bending, due to the DIFFERENCE
    // between the magnetic field and the nominal field used by the KTRAJ.  Returns the change in momentum
    // = integral of the 'external' force needed to keep the particle onto this trajectory over the specified range;
    template<class KTRAJ> VEC3 integrate(BFieldMap const& bfield, KTRAJ const& ktraj, TimeRange const& trange) {
      // take a fixed number of steps.  This may fail for long ranges FIXME!
      unsigned nsteps(10);
      double dt = trange.range()/nsteps;
      // now integrate
      VEC3 dmom;
      for(unsigned istep=0; istep< nsteps; istep++){
	double tstep = trange.begin() + (0.5+istep)*dt;
	VEC3 vel = ktraj.velocity(tstep);
	VEC3 db = bfield.fieldVect(ktraj.position3(tstep)) - ktraj.bnom(tstep);
	dmom += cbar()*ktraj.charge()*dt*vel.Cross(db);
      }
      return dmom;
    }

    // estimate how long in time from the given start time the trajectory position will stay within the given tolerance
    // compared to the true particle motion, given the true magnetic field.  This measures the impact of the KTRAJ nominal field being
    // different from the true field
    template<class KTRAJ> double rangeInTolerance(double tstart, BFieldMap const& bfield, KTRAJ const& ktraj, double tol,TimeDir tdir = TimeDir::forwards) {
      // compute scaling factor
      double spd = ktraj.speed(tstart);
      double sfac = fabs(cbar()*ktraj.charge()*spd*spd/ktraj.momentum(tstart));
      // estimate step size from initial BFieldMap difference
      VEC3 tpos = ktraj.position3(tstart);
      VEC3 bvec = bfield.fieldVect(tpos);
      auto db = (bvec - ktraj.bnom(tstart)).R();
      // estimate the step size for testing the position deviation.  This comes from 2 components:
      // the (static) difference in field, and the change in field along the trajectory
      double tstep(0.1); // nominal step
      // step increment from static difference from nominal field.  0.2 comes from sagitta geometry
      // protect against nominal field = exact field
      if(db > 1e-4) tstep = std::min(tstep,0.2*sqrt(tol/(sfac*db))); 
      VEC3 dBdt = bfield.fieldDeriv(tpos,ktraj.velocity(tstart));
      // the deviation goes as the cube root of the BFieldMap change.  0.5 comes from cosine expansion
      if(fabs(dBdt.R())>1e-6) tstep = std::min(tstep, 0.5*std::cbrt(tol/(sfac*dBdt.R()))); //
      // loop over the trajectory in fixed steps to compute integrals and domains.
      // step size is defined by momentum direction tolerance.
      double tend = tstart;
      double dx(0.0);
      // advance till spatial distortion exceeds position tolerance or we reach the range limit
      do{
	// increment the range
	tend += (tdir == TimeDir::forwards) ? tstep : -tstep;
	tpos = ktraj.position3(tend);
	bvec = bfield.fieldVect(tpos);
	// BFieldMap diff with nominal
	auto db = (bvec - ktraj.bnom(tend)).R();
	// spatial distortion accumulation; this goes as the square of the time times the field difference
	dx += sfac*(tend-tstart)*tstep*db;
      } while(fabs(dx) < tol && ktraj.range().inRange(tend));
//      std::cout << "tstep " << tstep << " tstart " << tstart << " tend " << tend  << std::endl;
      return tend;
    }
  }

}
#endif
