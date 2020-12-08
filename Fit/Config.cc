#include "KinKal/Fit/Config.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, MetaIterConfig miconfig ) {
      ost << "Meta-Iteration " << miconfig.miter_ << " temp " << miconfig.temp_;
      if(miconfig.updatemat_)
	ost << " Update Material Xings";
      if(miconfig.updatebfcorr_)
	ost << " Update BFieldMap Correction";
      if(miconfig.updatehits_){
	ost << " Update Hit Internals rith ";
	ost << miconfig.hitupdaters_.size() << " Hit updaters" << std::endl;
      }
      ost << " converge, diverge, oscillating dchisq " << miconfig.convdchisq_ 
      << " "<< miconfig.divdchisq_  
      << " "<< miconfig.oscdchisq_;
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
