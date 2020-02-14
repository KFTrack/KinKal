#ifndef KinKal_KalSite_hh
#define KinKal_KalSite_hh
//
// Primary class of the Kinematic Kalman fit.  This class owns the state describing
// the fit (hits, material interactions, BField corrections) and coordinates the
// fit processing.  It uses
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the kinematic trajectory used in this fit
//
#include "BTrk/KinKal/PTTraj.hh"
#include "BTrk/KinKal/KKSite.hh"
#include <list>
#include <memory>

namespace KinKal {
  template<class KTraj> class KKTrk {
    public:
      typedef KKSite<KTraj> KKSITE;
      enum FitStatus {current=0,needsfit,unconverged,failed}; // fit status

    private:
      FitStatus status_; // current fit status
      PTTraj<KTraj> reftraj_; // reference against which the derivatives were evaluated
      PTTraj<KTraj> fittraj_; // result of the current fit
      std::list<std::unique_ptr<KKSITE> > sites_; // sites used in this fit.


  };
}
#endif
