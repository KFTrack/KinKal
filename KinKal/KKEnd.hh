#ifndef KinKal_KKEnd_hh
#define KinKal_KKEnd_hh
//
// End of the KK fit effect set, used transiently to start the fit and trajectory building
//
#include "KinKal/KKEff.hh"
#include <stdexcept>
#include <limits>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class KKEnd : public KKEff<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef typename KTRAJ::PDATA PDATA; // forward derivative type
      typedef typename KKEFF::WDATA WDATA; // forward the typedef
      typedef typename KKEFF::KKDATA KKDATA;
      // provide interface
      virtual void update(PKTRAJ const& ref) override;
      virtual void update(PKTRAJ const& ref, MConfig const& mconfig) override { 
	vscale_ = mconfig.varianceScale(); // annealing scale for covariance deweighting, to avoid numerical effects
	return update(ref); }
      virtual double time() const override { return (tdir_ == TDir::forwards) ? -std::numeric_limits<double>::max() : std::numeric_limits<double>::max(); } // make sure this is always at the end
      virtual bool isActive() const override { return true; }
      virtual void process(KKDATA& kkdata,TDir tdir) override;
      virtual void append(PKTRAJ& fit) override;
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~KKEnd(){}
      // accessors
      TDir const& tDir() const { return tdir_; }
      double deWeighting() const { return dwt_; }
      KTRAJ const& endTraj() const { return endtraj_; }
      WDATA const& endEffect() const { return endeff_; }

      // construct from trajectory and direction.  Deweighting must be tuned to balance stability vs bias
      KKEnd(PKTRAJ const& pktraj,TDir tdir, double dweight=1e6); 
    private:
      TDir tdir_; // direction for this effect; note the early end points forwards, the late backwards
      double dwt_; // deweighting factor
      double vscale_; // variance scale (from annealing)
      WDATA endeff_; // wdata representation of this effect's constraint/measurement
      KTRAJ endtraj_; // cache of parameters at the end of processing this direction, used in traj creation
 };

  template <class KTRAJ> KKEnd<KTRAJ>::KKEnd(PKTRAJ const& pktraj, TDir tdir, double dweight) :
    tdir_(tdir) , dwt_(dweight), vscale_(1.0), endtraj_(tdir == TDir::forwards ? pktraj.front() : pktraj.back()){
      update(pktraj);
    }


  template<class KTRAJ> void KKEnd<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    if(tdir == tdir_) 
      // start the fit with the de-weighted info cached from the previous iteration or seed
      kkdata.append(endeff_);
    else
    // at the opposite end, cache the final parameters
      endtraj_.params() = kkdata.pData();
    KKEffBase::setStatus(tdir,KKEffBase::processed);
  }

  template<class KTRAJ> void KKEnd<KTRAJ>::update(PKTRAJ const& ref) {
    auto refend = ref.nearestPiece(time()).params();
    refend.covariance() *= (dwt_/vscale_);
    // convert this to a weight (inversion)
    endeff_ = WData<PKTRAJ::NParams()>(refend);
    KKEffBase::updateStatus();
  }

  template<class KTRAJ> void KKEnd<KTRAJ>::append(PKTRAJ& fit) {
    // if the fit is empty and we're going in the right direction, take the end cache and
    // seed the fit with it
    if(tdir_ == TDir::forwards) {
      if(fit.pieces().size() == 0){
	// start with a very large range 
	endtraj_.range() = TRange(-std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
	// append this to the (empty) fit
	fit.append(endtraj_);
      } else
	throw std::invalid_argument("Input PKTraj isn't empty");
    }
  }

  template<class KTRAJ> void KKEnd<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "KKEnd " << static_cast<KKEff<KTRAJ>const&>(*this) << " direction " << tDir() << " deweight " << deWeighting() << std::endl;
    ost << "EndTraj ";
    endTraj().print(ost,detail);
    if(detail > 0){
      ost << "EndWeight " << endeff_ << std::endl;
    }
  }
  
  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, KKEnd<KTRAJ> const& kkend) {
    kkend.print(ost,0);
    return ost;
  }


}
#endif
