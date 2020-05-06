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
#include "KinKal/KKConfig.hh"
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
      virtual double time() const = 0; // time of this effect
      virtual unsigned nDOF() const {return 0; }; // how/if this effect contributes to the measurement NDOF
      virtual bool isActive() const = 0; // whether this effect is/was used in the fit
       // Add this effect to the ongoing fit in a give direction.
      virtual void process(KKDATA& kkdata,TDir tdir) = 0;
      virtual double fitChi() const { return 0.0;} // unbiased chi contribution of this effect after fitting
      virtual double chisq(PDATA const& pdata) const { return 0.0;} // chisq contribution WRT parameters 
      // update this effect for a new refernce trajectory
      virtual void update(PKTRAJ const& ref) = 0;
      // update this effect for a new configuration and reference trajectory
      virtual void update(PKTRAJ const& ref, MConfig const& mconfig) = 0;
      // append this effects trajectory change (if appropriate)
      virtual void append(PKTRAJ& fit) {};
      virtual void print(std::ostream& ost=std::cout,int detail=0) const =0;
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
