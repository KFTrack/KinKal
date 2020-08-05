#ifndef KinKal_BField_hh
#define KinKal_BField_hh
// class defining a BField Map interface for use in KinKal.
#include "KinKal/Vectors.hh"
#include "KinKal/TRange.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "Math/SMatrix.h"
#include <vector>
#include <algorithm>
#include <cstdarg>
#include <cmath>

namespace KinKal {
  class BField {
    public:
      typedef ROOT::Math::SMatrix<double,3> Grad; // field gradient: ie dBi/d(x,y,z)
      // return value of the field at a point
      virtual Vec3 fieldVect(Vec3 const& position) const = 0; 
      // return BField gradient = dB_i/dx_j, at a given point
      virtual Grad fieldGrad(Vec3 const& position) const = 0;
      // return the BField derivative at a given point along a given velocity, WRT time
      virtual Vec3 fieldDeriv(Vec3 const& position, Vec3 const& velocity) const = 0;
      virtual ~BField(){}
  };

  // trivial instance of the above, used for testing
  class UniformBField : public BField {
    public:
      virtual Vec3 fieldVect(Vec3 const& position) const override { return fvec_; }
      virtual Grad fieldGrad(Vec3 const& position) const override { return Grad(); }
      virtual Vec3 fieldDeriv(Vec3 const& position, Vec3 const& velocity) const override { return Vec3(); }
      UniformBField(Vec3 const& bnom) : fvec_(bnom) {}
      UniformBField(double BZ) : UniformBField(Vec3(0.0,0.0,BZ)) {}
      virtual ~UniformBField(){}
    private:
      Vec3 fvec_; // constant field
  };

// use superposition to create a composite field
  class CompositeBField : public BField {
    public:
      virtual Vec3 fieldVect(Vec3 const& position) const override;
      virtual Grad fieldGrad(Vec3 const& position) const override;
      virtual Vec3 fieldDeriv(Vec3 const& position, Vec3 const& velocity) const override;
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
      GradBField(double b0, double b1, double zg0, double zg1);
      virtual Vec3 fieldVect(Vec3 const& position) const override;
      virtual Grad fieldGrad(Vec3 const& position) const override { return fgrad_; }
      virtual Vec3 fieldDeriv(Vec3 const& position, Vec3 const& velocity) const override;
      virtual ~GradBField(){}
    private:
      double b0_, b1_;
      double z0_; 
      double grad_; // gradient in tesla/mm, computed from the fvec values
      Grad fgrad_;
  };

}
#endif
