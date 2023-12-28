#ifndef KinKal_BFieldMap_hh
#define KinKal_BFieldMap_hh
// class defining a BFieldMap Map interface for use in KinKal.
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/PhysicalConstants.h"
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
      // is the point inside the range of this map?
      virtual bool inRange(VEC3 const& position) const = 0;
      virtual ~BFieldMap(){}
      virtual void print(std::ostream& os ) const = 0;
      BFieldMap(){}
      // disallow copy and equivalence
      BFieldMap(BFieldMap const& ) = delete;
      BFieldMap& operator =(BFieldMap const& ) = delete;
      // speed of light in units to convert Tesla to mm (bending radius)
      static double constexpr cbar() { return CLHEP::c_light/1000.0; }
      // templated interface for interacting with kinematic trajectory classes
      // how far can you go along the given kinematic trajectory till BField inhomogeneity makes the momentum accuracy out of (fractional) tolerance
      template<class KTRAJ> double rangeInTolerance(KTRAJ const& ktraj, double tstart, double tol) const;
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

  // estimate how long the momentum vector from the given trajectory will stay within the given (fractional) tolerance given the field spatial variation
  // ie mag(P_true(tstart+dt) - P_traj(tstart+dt)) < tol.  This is a 2nd order local calculation
  template<class KTRAJ> double BFieldMap::rangeInTolerance(KTRAJ const& ktraj, double tstart, double tol) const {
    auto tpos = ktraj.position3(tstart); // starting position
    double dp = ktraj.momentum(tstart)*tol; // fractional tolerance on momentum
    auto vel = ktraj.velocity(tstart); // starting velocity
    auto dBdt = fieldDeriv(tpos,vel); // change in field WRT time along this velocity
    double d2pdt2 = (dBdt.Cross(vel)).R()*cbar()*fabs(ktraj.charge()); // 2nd derivative of momentum due to B change along the path
    if(d2pdt2 > 1e-10)
      return sqrt(dp/d2pdt2);
    else
      return ktraj.range().range()
  }

  // trivial instance of the above, used for testing
  class UniformBFieldMap : public BFieldMap {
    public:
      VEC3 fieldVect(VEC3 const& position) const override { return fvec_; }
      Grad fieldGrad(VEC3 const& position) const override { return Grad(); }
      VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override { return VEC3(); }
      void print(std::ostream& os =std::cout) const override { os << "Uniform BField, B = " << fvec_ << std::endl; }
      UniformBFieldMap(VEC3 const& bnom) : fvec_(bnom) {}
      UniformBFieldMap(double BZ) : UniformBFieldMap(VEC3(0.0,0.0,BZ)) {}
      bool inRange(VEC3 const& position) const override { return true; };
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
      void print(std::ostream& os =std::cout) const override;
      bool inRange(VEC3 const& position) const override;
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
      bool inRange(VEC3 const& position) const override { return true; }
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
