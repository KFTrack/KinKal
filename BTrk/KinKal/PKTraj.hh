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
// construct from an initial piece, which also provides kinematic information
      PKTraj(KT const& piece) : PTTraj<KT>(piece), KTraj(piece.mass(),piece.charge()) {}
// base class overrides; these just rely on the PTTraj to find the appropriate piece
      virtual void dirVector(trajdir dir,double time,Vec3& unit) const override {
	PTTraj<KT>::nearestPiece(time).dirVector(dir,time,unit);
      }
      virtual void momentum(double time,Mom4& mom) const override {
	PTTraj<KT>::nearestPiece(time).momentum(time,mom);
      }

  };
}
#endif

