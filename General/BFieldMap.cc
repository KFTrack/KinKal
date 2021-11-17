#include "KinKal/General/BFieldMap.hh"

namespace KinKal {

   VEC3 CompositeBFieldMap::fieldVect(VEC3 const&position) const  {
     VEC3 fvec;
     for(auto const field : fields_ ){
       fvec += field->fieldVect(position);
     }
     return fvec;
   }

   BFieldMap::Grad CompositeBFieldMap::fieldGrad(VEC3 const&position) const  {
     Grad fgrad;
     for(auto const field : fields_ ){
       fgrad += field->fieldGrad(position);
     }
     return fgrad;
   }

   VEC3 CompositeBFieldMap::fieldDeriv(VEC3 const& position, VEC3 const& velocity) const  {
     VEC3 dBdt;
     for(auto const field : fields_ ){
       dBdt += field->fieldDeriv(position,velocity);
     }
     return dBdt;
   }

   void CompositeBFieldMap::print(std::ostream& os ) const {
     os << "Composite BField with constituents as follows:" << std::endl;
     for(auto const& field : fields_) {
       field->print(os);
     }
   }

   double CompositeBFieldMap::zMin() const {
     double retval = -std::numeric_limits<float>::max();
     for(auto const& field : fields_) {
       retval = std::max(retval, field->zMin());
     }
     return retval;
   }

   double CompositeBFieldMap::zMax() const {
     double retval = std::numeric_limits<float>::max();
     for(auto const& field : fields_) {
       retval = std::min(retval, field->zMax());
     }
     return retval;
   }

   GradientBFieldMap::GradientBFieldMap(double b0, double b1, double zg0, double zg1) :
     b0_(b0), b1_(b1), z0_(zg0), grad_((b1_ - b0_)/(zg1-zg0)) {
       fgrad_[0][0] = -0.5*grad_;
       fgrad_[1][1] = -0.5*grad_;
       fgrad_[2][2] = -grad_;
     }

   VEC3 GradientBFieldMap::fieldVect(VEC3 const&position) const  {
     return VEC3(-0.5*grad_*position.X(), -0.5*grad_*position.Y(), b0_ + grad_*(position.Z()-z0_));
   }

   VEC3 GradientBFieldMap::fieldDeriv(VEC3 const& position, VEC3 const& velocity) const  {
     return VEC3(-0.5*grad_*velocity.X(),-0.5*grad_*velocity.Y(),grad_*velocity.Z());
   }

}
