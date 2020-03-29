#ifndef KinKal_KKEff_hh
#define KinKal_KKEff_hh
//
// Class representing a discrete effect along the kinematic Kalman filter track fit
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the trajectory class representing the particle in this fit
//
#include "KinKal/TDir.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/KKData.hh"
#include <array>
#include <memory>

namespace KinKal {
// base class for KKEff, for untemplated functions and content
  class KKEffBase {
    public:
      enum Status{unprocessed=0,processed,updated,failed};
      // time of this effect 
      virtual double time() const = 0;
      virtual unsigned nDOF() const = 0;
      virtual bool isActive() const = 0;
      void updateStatus() { status_[0] = status_[1] = updated; }
      Status status(TDir tdir) const { return status_[static_cast<std::underlying_type<TDir>::type>(tdir)]; }
      KKEffBase() : status_{{unprocessed,unprocessed}} {}
      virtual ~KKEffBase(){}
    protected:
      void setStatus(TDir tdir, Status status) { status_[static_cast<std::underlying_type<TDir>::type>(tdir)] = status; }
      std::array<Status,2> status_; // status of processing in each direction
  };

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

}

#endif
