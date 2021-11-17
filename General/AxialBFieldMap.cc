//
#include "KinKal/General/AxialBFieldMap.hh"
#include <cmath>
#include <fstream>
namespace KinKal {

  AxialBFieldMap::AxialBFieldMap(double zmin, double zmax, std::vector<double> const& field) : zmin_(zmin), zmax_(zmax), axial_(field) {
    zstep_=(zmax_-zmin_)/double(field.size()-1);
  }

  AxialBFieldMap::AxialBFieldMap(std::string const& file) {
    // read the file
    std::ifstream spectrum_stream(file);
    if ( (spectrum_stream.rdstate() & std::ifstream::failbit ) != 0 ){
      std::string errmsg = std::string("can't open Axial field file ") + file;
      throw std::invalid_argument(errmsg.c_str());
    }
    std::string line;
    bool first(true);
    static std::string comment("#");
    double bz, zold(zmin_);
    axial_.reserve(1000);
    while (std::getline(spectrum_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        std::istringstream iss(line);
        if (iss >> zmax_ >> bz ) {
          if(first){
            first = false;
            zmin_ = zmax_;
          }
          if(axial_.size()>0){
            // check uniformity
            double dz = zmax_ - (zmin_ + (zmax_-zold)*axial_.size());
            if(fabs(dz) > 1e-5) throw std::invalid_argument("field spacing isn't uniform!");
          }
          zold = zmax_;
          axial_.push_back(bz);
        }
      }
    }
    zstep_=(zmax_-zmin_)/double(axial_.size()-1);
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
    double dz = zstep_;
    // patch for null spots
    while(fabs(db)<1e-8 && ihigh < axial_.size()-1){
      ihigh++;
      dz += zstep_;
      db = axial_[ihigh]-axial_[ilow];
    }
    double grad = db/dz;
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
    else if(low >= (int)axial_.size())
      return (size_t)axial_.size()-2;
    else
      return (size_t)low;
  }

}
