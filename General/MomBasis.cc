#include "KinKal/General/MomBasis.hh"
#include <stdexcept>
namespace KinKal {
  VEC3 MomBasis::direction(Direction tdir, VEC3 momdir) {
    VEC3 retval;
    if(tdir == momdir_){
      retval = momdir.Unit();
    } else {
      // check the direction; the results are undefined for momdir
      double perp2 = momdir.Perp2();
      if(perp2 <= 0.0)throw std::invalid_argument("Perpendicular directions undefined when momentum direction is along Z");
      if( tdir == phidir_){
        VEC3 phidir(momdir.Y(),-momdir.X(),0.0);
        retval = phidir.Unit();
      } else if( tdir == perpdir_){
        if(momdir.Z() == 0.0)throw std::invalid_argument("Perp directions undefined when momentum is perpendicular to Z");
        VEC3 perpdir = momdir;
        perpdir.SetZ(-perp2/momdir.Z());
        retval = perpdir.Unit();
      } else {
        throw std::invalid_argument("Unknown Direction" );
      }
    }
    return retval;
  }
}
