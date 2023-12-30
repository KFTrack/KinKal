#ifndef KinKal_DomainWall_hh
#define KinKal_DomainWall_hh
//
// Effect describing the change in fit parameters for the change in BField crossing between 2 domains
// This effect adds no information content or noise (presently), just transports the parameters
//
#include "KinKal/Fit/Effect.hh"
#include "KinKal/Fit/Domain.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/Fit/Config.hh"
#include <iostream>
#include <stdexcept>
#include <array>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class DomainWall : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      double time() const override { return domain_.begin(); } // set by convention that this bounds the upper domain
      bool active() const override { return true; } // always active
      void process(FitState& kkdata,TimeDir tdir) override;
      void updateState(MetaIterConfig const& miconfig,bool first) override {}; // nothing to do here
      void updateConfig(Config const& config) override {}
      void updateReference(PTRAJ const& ptraj) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      void append(PTRAJ& fit,TimeDir tdir) override;
      Chisq chisq(Parameters const& pdata) const override { return Chisq();} // no information added
      auto const& parameterChange() const { return dpfwd_; }
      virtual ~DomainWall(){}
      // disallow copy and equivalence
      DomainWall(DomainWall const& ) = delete;
      DomainWall& operator =(DomainWall const& ) = delete;
      // specific DomainWall interface
      // create from the domain and BField
      DomainWall(BFieldMap const& bfield,Domain const& domain,PTRAJ const& ptraj);
      auto const& domain() const { return domain_; }
      // time at the middle of the PREVIOUS domain (approximate)
      double prevTime() const { return domain_.begin() - 0.5*domain_.range(); }
       // time at the middle of the NEXT domain
      double nextTime() const { return domain_.mid(); }

    private:
      BFieldMap const& bfield_; // bfield
      Domain domain_;  // the upper domain bounded by this wall
      DVEC dpfwd_; // parameter change across this domain wall in the forwards time direction
  };

  template<class KTRAJ> DomainWall<KTRAJ>::DomainWall(BFieldMap const& bfield,Domain const& domain,PTRAJ const& ptraj) :
    bfield_(bfield), domain_(domain) {
      updateReference(ptraj);
    }

  template<class KTRAJ> void DomainWall<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    kkdata.append(dpfwd_,tdir);
    // rotate the covariance matrix for the change in reference BField.  TODO
  }

  template<class KTRAJ> void DomainWall<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    // sample BField in the domains bounded by this domain wall
    auto bnext = bfield_.fieldVect(ptraj.position3(nextTime()));
    auto const& refpiece = ptraj.nearestPiece(prevTime()); // by convention, use previous domains parameters to define the derivative
    dpfwd_ = refpiece.dPardB(this->time(),bnext);
  }

  template<class KTRAJ> void DomainWall<KTRAJ>::append(PTRAJ& ptraj,TimeDir tdir) {
    double etime = time();
    // make sure the piece is appendable
    if((tdir == TimeDir::forwards && ptraj.back().range().begin() > etime) ||
        (tdir == TimeDir::backwards && ptraj.front().range().end() < etime) )
      throw std::invalid_argument("DomainWall: Can't append piece");
    // The prepend direction is awkward, as it must assume the previous domain has the same range as this
    TimeRange newrange = (tdir == TimeDir::forwards) ? TimeRange(domain_.begin(),std::max(ptraj.range().end(),domain_.end())) :
      TimeRange(std::min(ptraj.range().begin(),domain_.begin()-domain_.range()),domain_.begin());
    auto const& oldpiece = (tdir == TimeDir::forwards) ? ptraj.back() : ptraj.front();
    KTRAJ newpiece(oldpiece);
    newpiece.range() = newrange;
    if( tdir == TimeDir::forwards){
      auto bnext = bfield_.fieldVect(oldpiece.position3(nextTime()));
      newpiece.setBNom(etime,bnext);
      ptraj.append(newpiece);
      // update the parameters for the next iteration
      dpfwd_ = newpiece.params().parameters() - oldpiece.params().parameters();
    } else {
      auto bprev = bfield_.fieldVect(oldpiece.position3(prevTime()));
      newpiece.setBNom(etime,bprev);
      ptraj.prepend(newpiece);
      // update the parameters for the next iteration
      dpfwd_ = oldpiece.params().parameters() - newpiece.params().parameters();
    }
  }

  template<class KTRAJ> void DomainWall<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "DomainWall " << static_cast<Effect<KTRAJ>const&>(*this);
    ost << " effect " << dpfwd_ << " domain " << domain_ << std::endl;
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, DomainWall<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }

}
#endif
