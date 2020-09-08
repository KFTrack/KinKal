#include "KinKal/KKConfig.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, MIConfig miconfig ) {
      ost << "Meta-Iteration " << miconfig.miter_ << " temp " << miconfig.temp_;
      if(miconfig.updatemat_)
	ost << " Update Material Xings";
      if(miconfig.updatebfcorr_)
	ost << " Update BField Correction";
      if(miconfig.updatehits_){
	ost << " Update Hit Internals with ";
	ost << miconfig.hitupdaters_.size() << " Hit updaters" << std::endl;
      }
      ost << " converge, diverge, oscillating dchisq " << miconfig.convdchisq_ 
      << " "<< miconfig.divdchisq_  
      << " "<< miconfig.oscdchisq_;
    return ost;
  }

  std::ostream& operator <<(std::ostream& ost, KKConfig kkconfig ) {
    ost << "KKConfig maxniter " << kkconfig.maxniter_ << " dweight " << kkconfig.dwt_
      << " min NDOF " << kkconfig.minndof_ 
      << " with " << kkconfig.schedule().size() << " Meta-iterations:" << std::endl;
    for(auto const& miconfig : kkconfig.schedule() ) {
      ost << miconfig << std::endl;
    }
    return ost;
  }
}
