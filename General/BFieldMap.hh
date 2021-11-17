#ifndef KinKal_BFieldMap_hh
#define KinKal_BFieldMap_hh
// class defining a BFieldMap Map interface for use in KinKal.
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/PhysicalConstants.h"
#include "KinKal/Fit/Config.hh"
#include "Math/SMatrix.h"
#include <vector>
#include <limits>
#include <algorithm>
#include <cstdarg>
#include <cmath>
#include <ostream>

namespace KinKal {
  class BFieldMap {
    public:
      using Grad = ROOT::Math::SMatrix<double,3>; // field gradient: ie dBi/d(x,y,z)
      // return value of the field at a point
      virtual VEC3 fieldVect(VEC3 const& position) const = 0;
      // return BFieldMap gradient = dB_i/dx_j, at a given point
      virtual Grad fieldGrad(VEC3 const& position) const = 0;
      // return the BFieldMap derivative at a given point along a given velocity, WRT time
      virtual VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const = 0;
      virtual ~BFieldMap(){}
      virtual void print(std::ostream& os ) const = 0;
      virtual double zMin() const = 0;
      virtual double zMax() const = 0;
      BFieldMap(){}
      // disallow copy and equivalence
      BFieldMap(BFieldMap const& ) = delete;
      BFieldMap& operator =(BFieldMap const& ) = delete;
      // speed of light in units to convert Tesla to mm (bending radius)
      static double constexpr cbar() { return CLHEP::c_light/1000.0; }
      // templated interface for interacting with kinematic trajectory classes
      // how far can you go along the given kinematic trajectory in the given direction from the start time till BField inhomogeneity effects are out-of tolerance
      template<class KTRAJ> double rangeInTolerance(KTRAJ const& ktraj, double tstart, double tol, bool local=true) const;
      // divide a kinematic trajectory range into magnetic 'domains' within which the BField inhomogeneity effects are within tolerance
      template<class KTRAJ> void setDomains(KTRAJ const& ktraj, TimeRange const& range, Config const& config, std::vector<TimeRange>& ranges) const;
      // integrate the residual magentic force over the given kinematic trajectory and range due to the difference between the true field and the nominal field in the
      template<class KTRAJ> VEC3 integrate(KTRAJ const& ktraj, TimeRange const& trange) const;
  };

  template<class KTRAJ> VEC3 BFieldMap::integrate(KTRAJ const& ktraj, TimeRange const& trange) const {
    // take a fixed number of steps.  This may fail for long ranges FIXME!
    unsigned nsteps(10);
    double dt = trange.range()/nsteps;
    // now integrate
    VEC3 dmom;
    for(unsigned istep=0; istep< nsteps; istep++){
      double tstep = trange.begin() + (0.5+istep)*dt;
      VEC3 vel = ktraj.velocity(tstep);
      VEC3 db = fieldVect(ktraj.position3(tstep)) - ktraj.bnom(tstep);
      dmom += cbar()*ktraj.charge()*dt*vel.Cross(db);
    }
    return dmom;
  }

  // estimate how long in time from the given start time the trajectory position will stay within the given tolerance
  // compared to the true particle motion, given the true magnetic field.  This measures the impact of the KTRAJ nominal field being
  // both fixed and different from the true field
  template<class KTRAJ> double BFieldMap::rangeInTolerance(KTRAJ const& ktraj, double tstart, double tol, bool local) const {
    VEC3 tpos = ktraj.position3(tstart);
    VEC3 dBdt = fieldDeriv(tpos,ktraj.velocity(tstart));
    double dbdt = dBdt.R();
    // if there's no gradient, don't extend
    if(dbdt < 1e-20) return ktraj.range().end()+1.0e-10;
    // compute scaling factor
    double spd = ktraj.speed(tstart);
    double sfac = fabs(cbar()*ktraj.charge()*spd*spd/ktraj.momentum(tstart));
    // estimate step size from initial BFieldMap difference
    VEC3 bref;
    if(local)
      bref = fieldVect(tpos);
    else// fixed reference comes from the traj
      bref = ktraj.bnom(tstart);
    // estimate the step size for testing the position deviation.  This comes from 2 components:
    // the (static) difference in field, and the change in field along the trajectory
    double tstep(0.1); // maximum step.  I shouldn't need this FIXME
    // step increment from static difference from nominal field.  0.2 comes from sagitta geometry
    // protect against nominal field = exact field
    //    auto db = (bvec - ktraj.bnom(tstart)).R();
    //    if(db > 1e-4) tstep = std::min(tstep,0.2*sqrt(tol/(sfac*db)));
    // use the magnitude as a worst-case
    // the deviation goes as the cube root of the BFieldMap change.  0.5 comes from cosine expansion
    // this calculation needs testing/fixing: the response to tolerance is not as expected FIXME!
    if(dbdt >1e-6) tstep = std::min(tstep, 0.5*std::cbrt(tol/(sfac*dbdt)));
    // loop over the trajectory in fixed steps to compute integrals and domains.
    // step size is defined by momentum direction tolerance.
    double tend = tstart;
    double dx(0.0);
    // advance till spatial distortion exceeds position tolerance or we reach the range limit
    VEC3 bvec;
    do{
      // increment the range
      tend += tstep;
      tpos = ktraj.position3(tend);
      bvec = fieldVect(tpos);
      // BFieldMap diff with reference
      auto db = (bvec - bref).R();
      // spatial distortion accumulation; this goes as the square of the time times the field difference
      dx += sfac*(tend-tstart)*tstep*db;
    } while(fabs(dx) < tol && ktraj.range().inRange(tend));
    return tend;
  }

  // divide a trajectory into magnetic 'domains' within which BField changes are within tolerance
  template<class KTRAJ> void BFieldMap::setDomains(KTRAJ const& ktraj, TimeRange const& range, Config const& config, std::vector<TimeRange>& ranges) const {
    double tstart, tend;
    tstart = range.begin();
    do {
      // see how far we can go on the current traj before the BField change causes it to go out of tolerance
      // that defines the end of this domain
      tend = rangeInTolerance(ktraj,tstart,config.tol_,config.localBFieldCorr());
      ranges.emplace_back(tstart,tend);
      // start the next domain and the end of this one
      tstart = tend;
    } while(tstart < range.end());
  }

  // trivial instance of the above, used for testing
  class UniformBFieldMap : public BFieldMap {
    public:
      VEC3 fieldVect(VEC3 const& position) const override { return fvec_; }
      Grad fieldGrad(VEC3 const& position) const override { return Grad(); }
      VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override { return VEC3(); }
      double zMin() const override { return -std::numeric_limits<float>::max(); }
      double zMax() const override { return std::numeric_limits<float>::max(); }
      void print(std::ostream& os =std::cout) const override { os << "Uniform BField, B = " << fvec_ << std::endl; }
      UniformBFieldMap(VEC3 const& bnom) : fvec_(bnom) {}
      UniformBFieldMap(double BZ) : UniformBFieldMap(VEC3(0.0,0.0,BZ)) {}
      virtual ~UniformBFieldMap(){}
      // disallow copy and equivalence
      UniformBFieldMap(UniformBFieldMap const& ) = delete;
      UniformBFieldMap& operator =(UniformBFieldMap const& ) = delete;
    private:
      VEC3 fvec_; // constant field
  };

  // use superposition to create a composite field
  class CompositeBFieldMap : public BFieldMap {
    public:
      using FCOL = std::vector<const BFieldMap*>;
      VEC3 fieldVect(VEC3 const& position) const override;
      Grad fieldGrad(VEC3 const& position) const override;
      VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override;
      double zMin() const override;
      double zMax() const override;
      void print(std::ostream& os =std::cout) const override;
      CompositeBFieldMap () {}
      CompositeBFieldMap(FCOL const& fields) : fields_(fields) {}
      void addField(BFieldMap const& field) { fields_.push_back(&field); }
      virtual ~CompositeBFieldMap() {}
      // disallow copy and equivalence
      CompositeBFieldMap(CompositeBFieldMap const& ) = delete;
      CompositeBFieldMap& operator =(CompositeBFieldMap const& ) = delete;

    private:
      FCOL fields_; // fields
  };

  // simple Z gradient field, used to test Field corrections
  class GradientBFieldMap : public BFieldMap {
    public:
      GradientBFieldMap(double b0, double b1, double zg0, double zg1);
      VEC3 fieldVect(VEC3 const& position) const override;
      Grad fieldGrad(VEC3 const& position) const override { return fgrad_; }
      VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override;
      double zMin() const override { return -std::numeric_limits<float>::max(); }
      double zMax() const override { return std::numeric_limits<float>::max(); }
      void print(std::ostream& os =std::cout) const override { os << "BField with  constant gradient of " << grad_ << " Tesla/mm" << std::endl; }
      double gradient() const { return grad_; }
      virtual ~GradientBFieldMap(){}
      // disallow copy and equivalence
      GradientBFieldMap(GradientBFieldMap const& ) = delete;
      GradientBFieldMap& operator =(GradientBFieldMap const& ) = delete;
    private:
      double b0_, b1_;
      double z0_;
      double grad_; // gradient in tesla/mm, computed from the fvec values
      Grad fgrad_;
  };

}
#endif
