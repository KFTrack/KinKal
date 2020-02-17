#ifndef KinKal_BField_hh
#define KinKal_BField_hh
// class defining a BField Map interface for use in KinKal.
#include "BTrk/KinKal/Types.hh"
namespace KinKal {
  class BField {
    public:
      virtual void fieldVect(Vec3& field, Vec3 const& position=Vec3()) const = 0; // nominal field defined at the origin
  };

  // trivial instance of the above, used for testing
  class UniformBField : public BField {
    public:
      virtual void fieldVect(Vec3& field, Vec3 const& position=Vec3()) const override { field = field_; }
      UniformBField(Vec3 const& bnom) : field_(bnom) {}
      UniformBField(double BZ) : UniformBField(Vec3(0.0,0.0,BZ)) {}
    private:
      Vec3 field_; // constant field
  };
}
#endif
