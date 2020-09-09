#ifndef KinKal_BFieldMap_hh
#define KinKal_BFieldMap_hh
// class defining a BFieldMap Map interface for use in KinKal.
#include "KinKal/Vectors.hh"
#include "KinKal/TimeRange.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "Math/SMatrix.h"
#include <vector>
#include <algorithm>
#include <cstdarg>
#include <cmath>

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
  };

  // trivial instance of the above, used for testing
  class UniformBFieldMap : public BFieldMap {
    public:
      virtual VEC3 fieldVect(VEC3 const& position) const override { return fvec_; }
      virtual Grad fieldGrad(VEC3 const& position) const override { return Grad(); }
      virtual VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override { return VEC3(); }
      UniformBFieldMap(VEC3 const& bnom) : fvec_(bnom) {}
      UniformBFieldMap(double BZ) : UniformBFieldMap(VEC3(0.0,0.0,BZ)) {}
      virtual ~UniformBFieldMap(){}
    private:
      VEC3 fvec_; // constant field
  };

// use superposition to create a composite field
  class CompositeBFieldMap : public BFieldMap {
    public:
      using FCOL = std::vector<const BFieldMap*>;
      virtual VEC3 fieldVect(VEC3 const& position) const override;
      virtual Grad fieldGrad(VEC3 const& position) const override;
      virtual VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override;
      CompositeBFieldMap () {}
      CompositeBFieldMap(FCOL const& fields) : fields_(fields) {}
      void addField(BFieldMap const& field) { fields_.push_back(&field); }
      virtual ~CompositeBFieldMap() {}
    private:
      FCOL fields_; // fields
  };

  // simple Z gradient field, used to test Field corrections
  class GradBFieldMap : public BFieldMap {
    public:
      GradBFieldMap(double b0, double b1, double zg0, double zg1);
      virtual VEC3 fieldVect(VEC3 const& position) const override;
      virtual Grad fieldGrad(VEC3 const& position) const override { return fgrad_; }
      virtual VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override;
      virtual ~GradBFieldMap(){}
    private:
      double b0_, b1_;
      double z0_; 
      double grad_; // gradient in tesla/mm, computed from the fvec values
      Grad fgrad_;
  };

}
#endif
