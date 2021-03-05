#include "KinKal/Fit/Config.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, MetaIterConfig const& miconfig ) {
      ost << "Meta-Iteration " << miconfig.miter_ << " temp " << miconfig.temp_;
      ost << " time precision " << miconfig.tprec_;
      ost << " converge, diverge delta-chisq," << miconfig.convdchisq_ << " "<< miconfig.divdchisq_ << " ";
      ost << miconfig.updaters_.size() << " Dedicated Updaters" << std::endl;
      return ost;
  }

  std::ostream& operator <<(std::ostream& ost, Config const& kkconfig ) {
    ost << "Config maxniter " << kkconfig.maxniter_ << " dweight " << kkconfig.dwt_
      << " min NDOF " << kkconfig.minndof_ << " BField correction " << kkconfig.bfcorr_
      << " with " << kkconfig.schedule().size() << " Meta-iterations:" << std::endl;
    for(auto const& miconfig : kkconfig.schedule() ) {
      ost << miconfig << std::endl;
    }
    return ost;
  }
}
