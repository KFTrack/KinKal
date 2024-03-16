#include "KinKal/Geometry/LineSegment.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

namespace KinKal {

  void LineSegment::print(std::ostream& ost, int detail) const {
    ost << *this << std::endl;
  }

  std::ostream& operator <<(std::ostream& ost, LineSegment const& lineseg) {
    ost << " LineSegment starting " << lineseg.start()
    << " direction " << lineseg.direction()
    << " length " << lineseg.length();
    return ost;
  }
}
