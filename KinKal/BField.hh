#ifndef KinKal_BField_hh
#define KinKal_BField_hh
// class defining a BField Map interface for use in KinKal.
#include "KinKal/Vectors.hh"
#include <vector>

namespace KinKal {
  class BField {
    public:
      virtual void fieldVect(Vec3& field, Vec3 const& position=Vec3()) const = 0; // nominal field defined at the origin
      virtual ~BField(){}
      // add interface for path integration FIXME!
  };

  // trivial instance of the above, used for testing
  class UniformBField : public BField {
    public:
      virtual void fieldVect(Vec3& fvec, Vec3 const& position=Vec3()) const override { fvec = fvec_; }
      UniformBField(Vec3 const& bnom) : fvec_(bnom) {}
      UniformBField(double BZ) : UniformBField(Vec3(0.0,0.0,BZ)) {}
      virtual ~UniformBField(){}
    private:
      Vec3 fvec_; // constant field
  };

// use superposition to create a composite field
  class CompositeBField : public BField {
    public:
      virtual void fieldVect(Vec3& fvec, Vec3 const& position=Vec3()) const override {
	fvec = Vec3();
	for(auto const field : fields_ ){
	  Vec3 temp;
	  field->fieldVect(temp,position);
	  fvec += temp;
	}
      }
      CompositeBField () {}
      void addField(BField const& field) { fields_.push_back(&field); }
    private:
      std::vector<const BField*> fields_; // fields
  };



  // simple Z gradient field, used to test Field corrections
  class GradBField : public BField {
    public:
      GradBField(double b0, double b1, double zg0, double zg1) :
	f0_(0.0,0.0,b0), f1_(0.0,0.0,b1), z0_(zg0), z1_(zg1), grad_((b1 - b0)/(zg1-zg0)) {}
      virtual void fieldVect(Vec3& fvec, Vec3 const& position=Vec3()) const override {
	if(position.z() < z0_)
	  fvec = f0_;
	else if(position.z() > z1_)
	  fvec = f1_; 
	else {
	  double bgrad = grad_*(position.z()-z0_);
	  // work in cylindrical coordinates
	  double bz = f0_.z()+bgrad;
	  double bx = -0.5*grad_*position.x();
	  double by = -0.5*grad_*position.y();
	  fvec = Vec3(bx,by,bz);
	}
      }
      virtual ~GradBField(){}
    private:
      Vec3 f0_, f1_;
      double z0_, z1_;
      double grad_; // gradient in tesla/mm, computed from the fvec values
  };
}
#endif
