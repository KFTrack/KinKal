#include "KinKal/Fit/Config.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, MetaIterConfig miconfig ) {
      ost << "Meta-Iteration " << miconfig.miter_ << " temp " << miconfig.temp_;
      ost << miconfig.updaters_.size() << " Dedicated Updaters" << std::endl;
      ost << " converge, diverge, oscillating dchisq " << miconfig.convdchisq_ 
      << " "<< miconfig.divdchisq_  
      << " "<< miconfig.oscdchisq_ << " time precision " << miconfig.tprec_;
    return ost;
  }

  std::ostream& operator <<(std::ostream& ost, Config kkconfig ) {
    ost << "Config maxniter " << kkconfig.maxniter_ << " dweight " << kkconfig.dwt_
      << " min NDOF " << kkconfig.minndof_ << " BField correction " << kkconfig.bfcorr_
      << " with " << kkconfig.schedule().size() << " Meta-iterations:" << std::endl;
    for(auto const& miconfig : kkconfig.schedule() ) {
      ost << miconfig << std::endl;
    }
    return ost;
  }
}
