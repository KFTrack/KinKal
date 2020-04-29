#ifndef KinKal_BField_hh
#define KinKal_BField_hh
// class defining a BField Map interface for use in KinKal.
#include "KinKal/Vectors.hh"
#include "KinKal/TRange.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "Math/SMatrix.h"
#include <vector>
#include <cstdarg>

namespace KinKal {
  class BField {
    public:
      typedef ROOT::Math::SMatrix<float,3> Grad; // field gradient: ie dBi/d(x,y,z)
      virtual void fieldVect(Vec3 const& position, Vec3& field) const = 0; 
      // return BField gradient = dB_i/dx_j, at a given point
      virtual void fieldGrad(Vec3 const& position, Grad& fgrad) const = 0;
      // return the BField derivative at a given point along a given velocity, WRT time
      virtual void fieldDeriv(Vec3 const& position, Vec3 const& velocity, Vec3& dBdt) const = 0;
      // integrate the magentic force over the given trajectory and range, due to the DIFFERENCE
      // between this field and the field vector referenced by the traj.  Returns the change in momentum
      // ( = force needed to force the particle onto this trajectory);
      template <class KTRAJ> void integrate(KTRAJ const& ktraj, TRange const& range, Vec3& dp) const;
      // speed of light in units to convert Tesla to mm
      static float cbar() { static float retval = CLHEP::c_light/1000.0; return retval; }

      virtual ~BField(){}
  };

  // trivial instance of the above, used for testing
  class UniformBField : public BField {
    public:
      virtual void fieldVect(Vec3 const&position, Vec3& fvec) const override { fvec = fvec_; }
      virtual void fieldGrad(Vec3 const&position, Grad& fgrad) const override { fgrad = Grad(); }
      virtual void fieldDeriv(Vec3 const& position, Vec3 const& velocity, Vec3& dBdt) const override { dBdt = Vec3(); }
      UniformBField(Vec3 const& bnom) : fvec_(bnom) {}
      UniformBField(float BZ) : UniformBField(Vec3(0.0,0.0,BZ)) {}
      virtual ~UniformBField(){}
    private:
      Vec3 fvec_; // constant field
  };

// use superposition to create a composite field
  class CompositeBField : public BField {
    public:
      virtual void fieldVect(Vec3 const&position, Vec3& fvec) const override;
      virtual void fieldGrad(Vec3 const&position, Grad& fgrad) const override;
      virtual void fieldDeriv(Vec3 const& position, Vec3 const& velocity, Vec3& dBdt) const override;
      CompositeBField () {}
      CompositeBField(int fcount, ...);
      void addField(BField const& field) { fields_.push_back(&field); }
      virtual ~CompositeBField() {}
    private:
      std::vector<const BField*> fields_; // fields
  };

  // simple Z gradient field, used to test Field corrections
  class GradBField : public BField {
    public:
      GradBField(float b0, float b1, float zg0, float zg1);
      virtual void fieldVect(Vec3 const&position, Vec3& fvec) const override;
      virtual void fieldGrad(Vec3 const&position, Grad& fgrad) const override;
      virtual void fieldDeriv(Vec3 const& position, Vec3 const& velocity, Vec3& dBdt) const override;
      virtual ~GradBField(){}
    private:
      float b0_, b1_;
      float z0_; 
      float grad_; // gradient in tesla/mm, computed from the fvec values
      Grad fgrad_;
  };

  template<class KTRAJ> void BField::integrate(KTRAJ const& ktraj, TRange const& trange, Vec3& dp) const {
    dp = Vec3();
    // compare this field with the noiminal at the range mid: only bother to integrate if this is outside
    // tolerance.  Smarter would be to test at the begining, middle and end FIXME!
    Vec3 midpos, midfield;
    ktraj.position(trange.mid(),midpos);
    fieldVect(midpos,midfield);
    auto db = midfield - ktraj.bnom(trange.mid());
    if(db.R() > 1e-6){ // tolerance should be a parameter FIXME!
      unsigned nsteps(10); // this should be computed from field rate-of-change over this range FIXME!
      float tstep = trange.range()/float(nsteps);
      for(unsigned istep=0; istep< nsteps; istep++){
	float tsamp = trange.low() + istep*tstep;
	Vec3 tpos, vel, bvec;
	ktraj.position(tsamp,tpos);
	ktraj.velocity(tsamp,vel);
	fieldVect(tpos,bvec);
	Vec3 db = bvec - ktraj.bnom(tsamp);
	dp += cbar()*ktraj.charge()*tstep*vel.Cross(db); // check sign FIXME!
      }
    }
  }

}
#endif
