#include "KinKal/KKConfig.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, KKConfig kkconfig ) {
    ost << "KKConfig maxniter " << kkconfig.maxniter_ << " dweight " << kkconfig.dwt_
      << " conv. dchisq " << kkconfig.convdchisq_ << " min NDOF " << kkconfig.minndof_ 
      << " with " << kkconfig.schedule().size() << " Meta-iterations:" << std::endl;
    unsigned imeta(0);
    for(auto const& mconfig : kkconfig.schedule() ) {
      ost << "Meta-iteration " << imeta++ << " temp " << mconfig.temp_;
      if(mconfig.processmat_)
	ost << " Material Processing";
      else
	ost << " No Material Processing";
      if(mconfig.processbfield_)
	ost << " BField Processing";
      else
	ost << " No BField Processing";
      if(mconfig.updatehits_)
	ost << " Hit Updating with ";
      else
	ost << " No Hit Updating with ";
      ost << mconfig.hitupdateparams_.size() << " Hit update parameters" << std::endl;
    }
    return ost;
  }
}
