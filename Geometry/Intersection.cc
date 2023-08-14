#include "KinKal/Geometry/Intersection.hh"
#include <stdexcept>
#include <iostream>
#include <cstdio>

std::ostream& operator <<(std::ostream& ost, KinKal::Intersection const& inter) {
  if(inter.onsurface_){
    ost << "Intersection at time  " << inter.time_ << " position " << inter.pos_ << " surface normal " << inter.norm_;
    if(inter.inbounds_){
      ost << " in bounds ";
    }
  }else {
    ost << "No Intersection";
  }
  return ost;
}
