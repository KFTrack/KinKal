#ifndef KinKal_BField_hh
#define KinKal_BField_hh
//
// Effect to correct the fit parameters for the change in BField along a small piece of the trajectory.
// This effect adds no information content or noise (presently), just transports the parameters
//
#include "KinKal/Fit/Effect.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/Fit/Config.hh"
#include <iostream>
#include <stdexcept>
#include <array>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class BField : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      double time() const override { return drange_.mid(); } // apply the correction at the middle of the range
      bool active() const override { return bfcorr_; }
      void process(FitState& kkdata,TimeDir tdir) override;
      void updateState(MetaIterConfig const& miconfig,bool first) override {}
      void updateConfig(Config const& config) override { bfcorr_ = config.bfcorr_; }
      void updateReference(KTRAJPTR const& ltrajptr) override {} // nothing explicit here
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      void append(PTRAJ& fit,TimeDir tdir) override;
      Chisq chisq(Parameters const& pdata) const override { return Chisq();}
      auto const& parameterChange() const { return dpfwd_; }
      virtual ~BField(){}
      // disallow copy and equivalence
      BField(BField const& ) = delete;
      BField& operator =(BField const& ) = delete;
      // create from the domain range, the effect, and the
      BField(Config const& config, BFieldMap const& bfield,TimeRange const& drange) :
        bfield_(bfield), drange_(drange), bfcorr_(config.bfcorr_) {}
      TimeRange const& range() const { return drange_; }

    private:
      BFieldMap const& bfield_; // bfield
      TimeRange drange_; // extent of this effect.  The middle is at the transition point between 2 bfield domains (domain transition)
      DVEC dpfwd_; // aggregate effect in parameter space of BFieldMap change over this domain in the forwards time direction
      bool bfcorr_; // apply correction or not
  };

  template<class KTRAJ> void BField<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    if(bfcorr_){
      kkdata.append(dpfwd_,tdir);
      // rotate the covariance matrix for the change in BField.  This requires 2nd derivatives TODO
    }
  }

  template<class KTRAJ> void BField<KTRAJ>::append(PTRAJ& ptraj,TimeDir tdir) {
    if(bfcorr_){
      double etime = time();
      // make sure the piece is appendable
      if((tdir == TimeDir::forwards && ptraj.back().range().begin() > etime) ||
          (tdir == TimeDir::backwards && ptraj.front().range().end() < etime) )
        throw std::invalid_argument("BField: Can't append piece");
      // assume the next domain has ~about the same range
      TimeRange newrange = (tdir == TimeDir::forwards) ? TimeRange(etime,std::max(ptraj.range().end(),drange_.end())) :
        TimeRange(std::min(ptraj.range().begin(),drange_.begin()),etime);
      // update the parameters according to the change in bnom across this domain
      // This corresponds to keeping the physical position and momentum constant, but referring to the BField
      // at the end vs the begining of the domain
      // Use the 1st order approximation: the exact method tried below doesn't do any better (slightly worse)
      VEC3 bend = (tdir == TimeDir::forwards) ? bfield_.fieldVect(ptraj.position3(drange_.end())) : bfield_.fieldVect(ptraj.position3(drange_.begin()));
      // update the parameter change due to the BField change.  Note this assumes the traj piece
      // at the begining of the domain has the same bnom as the BField at that point in space
      KTRAJ newpiece = (tdir == TimeDir::forwards) ? ptraj.back() : ptraj.front();
      newpiece.setBNom(etime,bend);
      newpiece.range() = newrange;
      // extract the parameter change for the next processing BEFORE appending
      // This should really be done in updateState, but it's easier here and has no knock-on effects
      dpfwd_ = (tdir == TimeDir::forwards) ? newpiece.params().parameters() - ptraj.back().params().parameters() : ptraj.back().params().parameters() - newpiece.params().parameters();
      if( tdir == TimeDir::forwards){
        ptraj.append(newpiece);
      } else {
        ptraj.prepend(newpiece);
      }
      // exact calculation (for reference)
      // extract the particle state at this transition
      // auto pstate = ptraj.back().stateEstimate(etime);
      // re-compute the trajectory at the domain end using this state
      // KTRAJ newpiece(pstate,bend,newrange);
      // set the parameter change for the next processing BEFORE appending
      // dpfwd_ = newpiece.params().parameters()-ptraj.back().params().parameters();
      // ptraj.append(newpiece);
    }
  }

  template<class KTRAJ> void BField<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "BField " << static_cast<Effect<KTRAJ>const&>(*this);
    ost << " effect " << dpfwd_ << " domain range " << drange_ << std::endl;
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, BField<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }

}
#endif
