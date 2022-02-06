#include "KinKal/Fit/MetaIterConfig.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, MetaIterConfig const& miconfig ) {
      ost << "Meta-Iteration temp " << miconfig.temp_;
      ost << " with " << miconfig.updaters_.size() << " Dedicated Updaters";
      return ost;
  }
}
