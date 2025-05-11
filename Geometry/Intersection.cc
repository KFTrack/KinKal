#include "KinKal/Geometry/Intersection.hh"
#include <stdexcept>
#include <iostream>
#include <cstdio>
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, KinKal::Intersection const& inter) {
    if(inter.onsurface_){
      ost << "Intersection on surface ";
      if(inter.inbounds_) ost << " in bounds ";
      if(inter.inrange_) ost << " in time range ";
    } else {
      ost << "No Intersection";
    }
    if(inter.good()){
      ost << " at time  " << inter.time_ << " position " << inter.pos_ << " surface normal " << inter.norm_;
    }
    return ost;
  }
}
