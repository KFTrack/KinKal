#ifndef KinKal_PKTraj_hh
#define KinKal_PKTraj_hh
//
//  Piecewise trajectory.  Templated on a KTraj.
//  used as part of the kinematic kalman fit
//
#include "BTrk/KinKal/PTTraj.hh"
#include "BTrk/KinKal/KTraj.hh"
#include <deque>
namespace KinKal {

  template <class KT> class PKTraj : public PTTraj<KT>, public KTraj {
    public:
      // base class implementation
// construct from an initial piece
      PKTraj(KT const& piece);
// base class overrides
      virtual void dirVector(trajdir dir,double time,Vec3& unit) const override;
      virtual void momentum(double t,Mom4& mom) const override; 
// specialize the piece accessor
      KT const& nearestPiece(double time) const;
  };
}
#endif

