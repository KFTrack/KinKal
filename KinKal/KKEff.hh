#ifndef KinKal_KKEff_hh
#define KinKal_KKEff_hh
//
// Class representing a discrete effect along the kinematic Kalman filter track fit
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the trajectory class representing the particle in this fit
//
#include "KinKal/PKTraj.hh"
#include "KinKal/KKData.hh"
#include "KinKal/KKEffBase.hh"
#include <array>
#include <memory>
#include <ostream>

namespace KinKal {

  template<class KTRAJ> class KKEff : public KKEffBase {
    public:
      // type of the data payload used for processing the fit
      typedef KKData<KTRAJ::PDATA::PDim()> KKDATA;
      typedef WData<KTRAJ::PDATA::PDim()> WDATA;
      typedef typename KTRAJ::PDATA PDATA;
      typedef PKTraj<KTRAJ> PKTRAJ;
      // Add this effect to the ongoing fit in a give direction.  Return value indicates success
      virtual bool process(KKDATA& kkdata,TDir tdir) = 0;
      virtual double chisq(PDATA const& pars) const = 0; // compute chisquared WRT some parameters
      // update this effect for a new refernce trajectory.  This must be overriden, but the base class implementation is still useful
      virtual bool update(PKTRAJ const& ref) = 0;
      // append this effects trajectory change (if appropriate)
      virtual bool append(PKTRAJ& fit) = 0;
      virtual ~KKEff(){} 
    protected:
      KKEff() {}
  };
  
  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, KKEff<KTRAJ> const& eff) {
    ost << (eff.isActive() ? "Active " : "Inactive ") << "time " << eff.time() << " status " <<
    TDir::forwards << " " << KKEffBase::statusName(eff.status(TDir::forwards))  << " : " <<
    TDir::backwards << " " << KKEffBase::statusName(eff.status(TDir::backwards));
    return ost;
  }

}

#endif
