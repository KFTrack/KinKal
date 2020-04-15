#include "KinKal/KKConfig.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, KKConfig kkconfig ) {
    ost << "KKConfig maxniter " << kkconfig.maxniter_ << " dweight " << kkconfig.dwt_
      << " conv. dchisq " << kkconfig.convdchisq_ << " min NDOF " << kkconfig.minndof_ 
      << " with " << kkconfig.schedule().size() << " Meta-iterations:" << std::endl;
    unsigned imeta(0);
    for(auto const& mconfig : kkconfig.schedule() ) {
      ost << "Meta-iteration " << imeta++ << " temp " << mconfig.temp_;
      if(mconfig.updatemat_)
	ost << " Update Material Xings";
      if(mconfig.updatebfcorr_)
	ost << " Update BField Correction";
      if(mconfig.updatehits_){
	ost << " Update Hit Internals with ";
	ost << mconfig.hitupdateparams_.size() << " Hit update parameters" << std::endl;
      }
    }
    return ost;
  }
}
