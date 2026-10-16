#include "KinKal/Fit/MetaIterConfig.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, MetaIterConfig const& miconfig ) {
      ost << "Meta-Iteration temp " << miconfig.temperature();
      ost << " with " << miconfig.nUpdaters() << " Dedicated Updaters";
      ost << ", reference piece hits " << miconfig.hitRefPiece() << " xings " << miconfig.xingRefPiece();
      return ost;
  }
}
