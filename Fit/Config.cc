#include "KinKal/Fit/Config.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, MetaIterConfig const& miconfig ) {
      ost << "Meta-Iteration " << miconfig.miter_ << " temp " << miconfig.temp_;
      ost << " converge, diverge delta-chisq," << miconfig.convdchisq_ << " "<< miconfig.divdchisq_ << " ";
      ost << miconfig.updaters_.size() << " Dedicated Updaters";
      return ost;
  }

  std::ostream& operator <<(std::ostream& ost, Config const& kkconfig ) {
    ost << "Config maxniter " << kkconfig.maxniter_
      << " dweight " << kkconfig.dwt_
      << " dpar chisq " << kkconfig.pdchi2_
      << " time buffer " << kkconfig.tbuff_
      << " fractional momentum tolerance " << kkconfig.tol_
      << " min NDOF " << kkconfig.minndof_
      << " BField correction " << kkconfig.bfcorr_
      << " with " << kkconfig.schedule().size()
      << " Meta-iterations:" << std::endl;
    for(auto const& miconfig : kkconfig.schedule() ) {
      ost << miconfig << std::endl;
    }
    return ost;
  }
}
