#include "KinKal/KKConfig.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, MConfig mconfig ) {
      ost << "Meta-Iteration " << mconfig.miter_ << " temp " << mconfig.temp_;
      if(mconfig.updatemat_)
	ost << " Update Material Xings";
      if(mconfig.updatebfcorr_)
	ost << " Update BField Correction";
      if(mconfig.updatehits_){
	ost << " Update Hit Internals with ";
	ost << mconfig.hitupdateparams_.size() << " Hit update parameters" << std::endl;
      }
      ost << " converge, diverge, oscillating dchisq " << mconfig.convdchisq_ 
      << " "<< mconfig.divdchisq_  
      << " "<< mconfig.oscdchisq_;
    return ost;
  }

  std::ostream& operator <<(std::ostream& ost, KKConfig kkconfig ) {
    ost << "KKConfig maxniter " << kkconfig.maxniter_ << " dweight " << kkconfig.dwt_
      << " min NDOF " << kkconfig.minndof_ 
      << " with " << kkconfig.schedule().size() << " Meta-iterations:" << std::endl;
    for(auto const& mconfig : kkconfig.schedule() ) {
      ost << mconfig << std::endl;
    }
    return ost;
  }
}
