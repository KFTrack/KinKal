#include "KinKal/Fit/Config.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, Config const& kkconfig ) {
    ost << "Config maxniter " << kkconfig.maxniter_
      << " dweight " << kkconfig.dwt_
      << " converge dchisq/dof " << kkconfig.convdchisq_
      << " diverge dchisq/dof " << kkconfig.divdchisq_
      << " diverge dpar chisq " << kkconfig.pdchisq_
      << " diverge traj gap (mm) " << kkconfig.divgap_
      << " fractional momentum tolerance " << kkconfig.tol_
      << " min domain step (ns) " << kkconfig.mindtstep_
      << " min field (T) " << kkconfig.minfield_
      << " domain margin (ns) " << kkconfig.domainmargin_
      << " max domains " << kkconfig.maxdomains_
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
