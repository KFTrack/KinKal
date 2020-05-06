#include "KinKal/BField.hh"

namespace KinKal {
   CompositeBField::CompositeBField(int fcount, ...) {
     std::va_list args;
     va_start(args,fcount);
     for (int ifield=0;ifield<fcount;ifield++) {
       fields_.push_back(va_arg(args,const BField*));
     }
   }

   Vec3 CompositeBField::fieldVect(Vec3 const&position) const  {
     Vec3 fvec;
     for(auto const field : fields_ ){
       fvec += field->fieldVect(position);
     }
     return fvec;
   }

   BField::Grad CompositeBField::fieldGrad(Vec3 const&position) const  {
     Grad fgrad;
     for(auto const field : fields_ ){
       fgrad += field->fieldGrad(position);
     }
     return fgrad;
   }

   Vec3 CompositeBField::fieldDeriv(Vec3 const& position, Vec3 const& velocity) const  {
     Vec3 dBdt;
     for(auto const field : fields_ ){
       dBdt += field->fieldDeriv(position,velocity);
     }
     return dBdt;
   }

   GradBField::GradBField(double b0, double b1, double zg0, double zg1) :
     b0_(b0), b1_(b1), z0_(zg0), grad_((b1_ - b0_)/(zg1-zg0)) {
       std::cout << "BGrad = " << grad_ << std::endl;
       fgrad_[0][0] = -0.5*grad_;
       fgrad_[1][1] = -0.5*grad_;
       fgrad_[2][2] = -grad_;
     }

   Vec3 GradBField::fieldVect(Vec3 const&position) const  {
     return Vec3(-0.5*grad_*position.X(), -0.5*grad_*position.Y(), b0_ + grad_*(position.Z()-z0_));
   }

   Vec3 GradBField::fieldDeriv(Vec3 const& position, Vec3 const& velocity) const  {
     return Vec3(-0.5*grad_*velocity.X(),-0.5*grad_*velocity.Y(),grad_*velocity.Z());
   }

}
