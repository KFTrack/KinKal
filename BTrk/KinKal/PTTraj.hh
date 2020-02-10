#ifndef KinKal_PTTraj_hh
#define KinKal_PTTraj_hh
//
//  TTraj subclass describing a piecewise trajectory.  Templated on a
//  simple TTraj.
//  used as part of the kinematic kalman fit
//
#include "BTrk/KinKal/TTraj.hh"
#include "BTrk/KinKal/TrkDirection.hh"
#include <deque>
namespace KinKal {

  template <class TT> class PTTraj : public TTraj{
    public:
      // base class implementation
      virtual void position(Vec4& pos) const override;
      virtual void position(double time, Vec3& pos) const override;
      virtual void velocity(double time, Vec3& vel) const override;
      virtual void direction(double time, Vec3& dir) const override;
// construct from an initial piece
      PTTraj(TT const& piece);
// append or prepend a piece
      addPiece(TT const& newpiece, trkDirection tdir=KinKal::tpos);
// Find the piece associated with a particular time
      virtual TT const& nearestPiece(double time) const;

    private:
      deque<TT> pieces_; // constituent pieces
  };
}
#endif

