#ifndef KinKal_TrackEnd_hh
#define KinKal_TrackEnd_hh
//
// End of the track effect.  This initiates the fit and trajectory building
//
#include "KinKal/Fit/Effect.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/General/BFieldMap.hh"
#include <stdexcept>
#include <limits>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class TrackEnd : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      // provide interface
      void update(PKTRAJ const& ref) override;
      void update(PKTRAJ const& ref, MetaIterConfig const& miconfig) override {
        vscale_ = miconfig.varianceScale(); // annealing scale for covariance deweighting, to avoid numerical effects
        return update(ref); }
      double time() const override { return (tdir_ == TimeDir::forwards) ? -std::numeric_limits<double>::max() : std::numeric_limits<double>::max(); } // make sure this is always at the end
      bool active() const override { return true; }
      void process(FitState& kkdata,TimeDir tdir) override;
      void append(PKTRAJ& fit) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~TrackEnd(){}
      // accessors
      TimeDir const& tDir() const { return tdir_; }
      double deWeighting() const { return config_.dwt_; }
      KTRAJ const& endTraj() const { return endtraj_; }
      Weights const& endEffect() const { return endeff_; }

      // construct from trajectory and direction.  Deweighting must be tuned to balance stability vs bias
      TrackEnd(Config const& config, BFieldMap const& bfield, PKTRAJ const& pktraj,TimeDir tdir);
      // disallow
      TrackEnd() = delete;
      TrackEnd(TrackEnd const& other) = delete;
      TrackEnd& operator =(TrackEnd const& other) = delete;
    private:
      Config const& config_; // cache configuration
      BFieldMap const& bfield_; // BField; needed to define reference
      TimeDir tdir_; // direction for this effect; note the early end points forwards, the late backwards
      VEC3 bnom_; // nominal BField
      double vscale_; // variance scale (from annealing)
      Weights endeff_; // wdata representation of this effect's constraint/measurement
      KTRAJ endtraj_; // cache of parameters at the end of processing this direction, used in traj creation
      static double tbuff_; // buffer to range end
  };

  template <class KTRAJ> double TrackEnd<KTRAJ>::tbuff_ = 1.0; // this should come from the config FIXME!

  template <class KTRAJ> TrackEnd<KTRAJ>::TrackEnd(Config const& config, BFieldMap const& bfield, PKTRAJ const& pktraj, TimeDir tdir) :
    config_(config), bfield_(bfield), tdir_(tdir) , vscale_(1.0),
    endtraj_(tdir == TimeDir::forwards ? pktraj.front() : pktraj.back()){
      update(pktraj);
    }

  template<class KTRAJ> void TrackEnd<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    if(tdir == tdir_)
      // start the fit with the de-weighted info cached from the previous iteration or seed
      kkdata.append(endeff_);
    else
      // at the opposite end, cache the final parameters
      endtraj_.params() = kkdata.pData();
    KKEFF::setState(tdir,KKEFF::processed);
  }

  template<class KTRAJ> void TrackEnd<KTRAJ>::update(PKTRAJ const& ref) {
    auto refend = (tdir_ == TimeDir::forwards) ? ref.front().params() : ref.back().params();
    refend.covariance() *= (config_.dwt_/vscale_);
    // convert this to a weight (inversion)
    endeff_ = Weights(refend);
    // set the range; this should buffer the original traj
    if(tdir_ == TimeDir::forwards){
      endtraj_.setRange(TimeRange(ref.range().begin()-tbuff_,ref.range().end()));
    } else {
      endtraj_.setRange(TimeRange(ref.range().begin(),ref.range().end()+tbuff_));
    }
    // update BField reference
    double endtime = (tdir_ == TimeDir::forwards) ? ref.range().begin() : ref.range().end();
    bnom_ = bfield_.fieldVect(ref.position3(endtime));
    KKEFF::updateState();
  }

  template<class KTRAJ> void TrackEnd<KTRAJ>::append(PKTRAJ& fit) {
    // Test the fit is empty and we're going in the right direction
    if(tdir_ == TimeDir::forwards) {
      if(fit.pieces().size() == 0){
        // take the end cache and  seed the fit with it
        // if we're using local BField, update accordingly
        if(config_.bfcorr_ == Config::variable || config_.bfcorr_ == Config::both){
          endtraj_.setBNom(endtraj_.range().begin(),bnom_);
        }
        // append this to the (empty) fit
        fit.append(endtraj_);
      } else
        throw std::invalid_argument("Input ParticleTrajectory isn't empty");
    }
  }

  template<class KTRAJ> void TrackEnd<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "TrackEnd " << static_cast<Effect<KTRAJ>const&>(*this) << " direction " << tDir() << " deweight " << deWeighting() << std::endl;
    ost << "EndTraj ";
    endTraj().print(ost,detail);
    if(detail > 0){
      ost << "EndWeight " << endeff_ << std::endl;
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, TrackEnd<KTRAJ> const& kkend) {
    kkend.print(ost,0);
    return ost;
  }


}
#endif
