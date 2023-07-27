//
//  Description of an annular planar section
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Disk_hh
#define KinKal_Disk_hh
#include "KinKal/Geometry/Annulus.hh"
namespace KinKal {
  class Disk : public Annulus {
    public:
      // construct from necessary parameters
      Disk(VEC3 const& norm,VEC3 const& udir, VEC3 const& center, double radius) : Annulus(norm,udir,center,0.0,radius) {}
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Disk const& disk);
#endif
