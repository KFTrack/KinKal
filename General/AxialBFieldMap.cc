//
#include "KinKal/General/AxialBFieldMap.hh"
#include <cmath>
namespace KinKal {

  AxialBFieldMap::AxialBFieldMap(double zmin, double zmax, std::vector<double> const& field) : zmin_(zmin), zmax_(zmax), axial_(field) {
    zstep_=(zmax_-zmin_)/double(field.size()-1);
  }
  VEC3 AxialBFieldMap::fieldVect(VEC3 const& position) const {
    size_t ilow = lowBound(position.Z());
    size_t ihigh = ilow+1;
    double db = axial_[ihigh]-axial_[ilow];
    double grad = db/zstep_;
    return VEC3(-0.5*grad*position.X(), -0.5*grad*position.Y(), axial_[ilow] + grad*(position.Z()-zval(ilow)));
  }

  AxialBFieldMap::Grad AxialBFieldMap::fieldGrad(VEC3 const& position) const {
    size_t ilow = lowBound(position.Z());
    size_t ihigh = ilow+1;
    double db = axial_[ihigh]-axial_[ilow];
    double grad = db/zstep_;
    Grad retval;
    retval[0][0] = retval[1][1] = -0.5*grad;
    retval[2][2] = -grad;
    return retval;
  }

  VEC3 AxialBFieldMap::fieldDeriv(VEC3 const& position, VEC3 const& velocity) const {
    size_t ilow = lowBound(position.Z());
    size_t ihigh = ilow+1;
    double db = axial_[ihigh]-axial_[ilow];
    double grad = db/zstep_;
    return VEC3(-0.5*grad*velocity.X(),-0.5*grad*velocity.Y(),grad*velocity.Z());
  }

  size_t AxialBFieldMap::lowBound(double zval) const {
    int low = (int)std::floor((zval-zmin_)/zstep_);
    if(low < 0 )
      return (size_t)0;
    else if(low >= axial_.size())
      return (size_t)axial_.size()-2;
    else
      return (size_t)low;
  }

}
