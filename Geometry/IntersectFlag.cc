#include "KinKal/Geometry/IntersectFlag.hh"
#include <stdexcept>
#include <iostream>
#include <cstdio>

std::ostream& operator <<(std::ostream& ost, KinKal::IntersectFlag const& iflag) {
  if(iflag.onsurface_){
    ost << "On surface ";
    if(iflag.inbounds_)ost << " in bounds ";
    if(iflag.inrange_)ost << " in range ";
  }else {
    ost << "No Intersection";
  }
  return ost;
}
