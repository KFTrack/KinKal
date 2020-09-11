#ifndef KinKal_TrackEnd_hh
#define KinKal_TrackEnd_hh
//
// End of the KK fit effect set, used transiently to start the fit and trajectory building
//
#include "KinKal/Effect.hh"
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
      bool isActive() const override { return true; }
      void process(FitState& kkdata,TimeDir tdir) override;
      void append(PKTRAJ& fit) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~TrackEnd(){}
      // accessors
      TimeDir const& tDir() const { return tdir_; }
      double deWeighting() const { return dwt_; }
      KTRAJ const& endTraj() const { return endtraj_; }
      Weights const& endEffect() const { return endeff_; }

      // construct from trajectory and direction.  Deweighting must be tuned to balance stability vs bias
      TrackEnd(PKTRAJ const& pktraj,TimeDir tdir, double dweight=1e6); 
    private:
      TimeDir tdir_; // direction for this effect; note the early end points forwards, the late backwards
      double dwt_; // deweighting factor
      double vscale_; // variance scale (from annealing)
      Weights endeff_; // wdata representation of this effect's constraint/measurement
      KTRAJ endtraj_; // cache of parameters at the end of processing this direction, used in traj creation
 };

  template <class KTRAJ> TrackEnd<KTRAJ>::TrackEnd(PKTRAJ const& pktraj, TimeDir tdir, double dweight) :
    tdir_(tdir) , dwt_(dweight), vscale_(1.0), endtraj_(tdir == TimeDir::forwards ? pktraj.front() : pktraj.back()){
      update(pktraj);
    }


  template<class KTRAJ> void TrackEnd<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    if(tdir == tdir_) 
      // start the fit with the de-weighted info cached from the previous iteration or seed
      kkdata.append(endeff_);
    else
    // at the opposite end, cache the final parameters
      endtraj_.params() = kkdata.pData();
    KKEFF::setStatus(tdir,KKEFF::processed);
  }

  template<class KTRAJ> void TrackEnd<KTRAJ>::update(PKTRAJ const& ref) {
    auto refend = ref.nearestPiece(time()).params();
    refend.covariance() *= (dwt_/vscale_);
    // convert this to a weight (inversion)
    endeff_ = Weights(refend);
    KKEFF::updateStatus();
  }

  template<class KTRAJ> void TrackEnd<KTRAJ>::append(PKTRAJ& fit) {
    // if the fit is empty and we're going in the right direction, take the end cache and
    // seed the fit with it
    if(tdir_ == TimeDir::forwards) {
      if(fit.pieces().size() == 0){
	// start with a very large range 
	endtraj_.range() = TimeRange(-std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
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
