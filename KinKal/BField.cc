#include "KinKal/BField.hh"

namespace KinKal {
   CompositeBField::CompositeBField(int fcount, ...) {
     std::va_list args;
     va_start(args,fcount);
     for (int ifield=0;ifield<fcount;ifield++) {
       fields_.push_back(va_arg(args,const BField*));
     }
   }
   void CompositeBField::fieldVect(Vec3 const&position, Vec3& fvec) const  {
     fvec = Vec3();
     for(auto const field : fields_ ){
       Vec3 temp;
       field->fieldVect(position,temp);
       fvec += temp;
     }
   }
   void CompositeBField::fieldGrad(Vec3 const&position, Grad& fgrad) const  {
     fgrad = Grad();
     for(auto const field : fields_ ){
       Grad temp;
       field->fieldGrad(position,temp);
       fgrad += temp;
     }
   }
   void CompositeBField::fieldDeriv(Vec3 const& position, Vec3 const& velocity, Vec3& dBdt) const  {
     dBdt = Vec3();
     for(auto const field : fields_ ){
       Vec3 dbdt;
       field->fieldDeriv(position,velocity,dbdt);
       dBdt += dbdt;
     }
   }

   GradBField::GradBField(float b0, float b1, float zg0, float zg1) :
     b0_(b0), b1_(b1), z0_(zg0), grad_((b1_ - b0_)/(zg1-zg0)) {
       std::cout << "BGrad = " << grad_ << std::endl;
       fgrad_[0][0] = -0.5*grad_;
       fgrad_[1][1] = -0.5*grad_;
       fgrad_[2][2] = -grad_;
     }
   void GradBField::fieldVect(Vec3 const&position, Vec3& fvec) const  {
     float bz = b0_ + grad_*(position.z()-z0_);
     float bx = -0.5*grad_*position.x();
     float by = -0.5*grad_*position.y();
     fvec = Vec3(bx,by,bz);
   }
   void GradBField::fieldGrad(Vec3 const&position, Grad& fgrad) const  {
      fgrad = fgrad_;
   }
   void GradBField::fieldDeriv(Vec3 const& position, Vec3 const& velocity, Vec3& dBdt) const  {
     dBdt = Vec3(-0.5*grad_*velocity.X(),-0.5*grad_*velocity.Y(),grad_*velocity.Z());
   }

}
