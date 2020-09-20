#ifndef KinKal_DetectorHit_hh
#define KinKal_DetectorHit_hh
//
//  Base class to describe a measurement that constraints some parameters of the fit
//  DetectorHits must have a measurement value (WRT a reference trajectory) and covariance, but can internally
//  be of any dimension and constrain any physical aspect of the fit (time, position, time+position, momentum, ...)
//  The measurement may be associated with a piece of detector material as well
//  Used as part of the kinematic Kalman fit
//
#include "General/Weights.hh"
#include "General/Parameters.hh"
#include "Detector/DetectorXing.hh"
#include "Trajectory/ParticleTrajectory.hh"
#include "Fit/Config.hh"
#include <memory>
#include <ostream>

namespace KinKal {
  template <class KTRAJ> class DetectorHit {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using DXING = DetectorXing<KTRAJ>;
      using DXINGPTR = std::shared_ptr<DXING>;
     // default
      DetectorHit(){}
      virtual ~DetectorHit(){}
      // disallow copy and equivalence
      DetectorHit(DetectorHit const& ) = delete; 
      DetectorHit& operator =(DetectorHit const& ) = delete;
      // compute the constraint this hit implies WRT the current reference, expressed as a weight
      virtual Weights weight() const = 0;
      // number of degrees of freedom constrained by this measurement (typically 1)
      virtual unsigned nDOF() const = 0;
      // compute the distance between this measurement and some reference parameters, scaled by errors
      virtual double chi(Parameters const& pdata) const = 0;
      // time of this hit: this is WRT the reference trajectory
      virtual double time() const = 0;
      // update the internals of the hit, specific to this meta-iteraion
      virtual void update(PKTRAJ const& pktraj, MetaIterConfig const& config) = 0;
      // update to a new reference, without changing any conditions
      virtual void update(PKTRAJ const& pktraj) = 0;
      // hits may inactive
      virtual bool isActive() const =0;
      // associated material information; null means no material
      virtual DXINGPTR const& detXingPtr() const = 0;
      bool hasMaterial() const { return (bool)detXingPtr(); }
      virtual void print(std::ostream& ost=std::cout,int detail=0) const = 0;
  };

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, DetectorHit<KTRAJ> const& thit) {
    thit.print(ost,0);
    return ost;
 }

}
#endif

