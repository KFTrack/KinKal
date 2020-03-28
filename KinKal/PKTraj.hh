#ifndef KinKal_PKTraj_hh
#define KinKal_PKTraj_hh
//
//  Piecewise trajectory with kinematic information.  Templated on a KTraj.
//  used as part of the kinematic kalman fit
//
#include "KinKal/PTTraj.hh"
#include "KinKal/KInter.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class PKTraj : public PTTraj<KTRAJ>, public KInter {
    public:
      typedef typename KTRAJ::PDER PDER; // forward derivative type from the 0th traj parameters
      constexpr static size_t NParams() { return KTRAJ::NParams(); } // forward the parameter space dimension
      // base class implementation
      // construct from an initial piece, which also provides kinematic information
      PKTraj(KTRAJ const& piece) : PTTraj<KTRAJ>(piece), KInter(piece.mass(),piece.charge()) {}
      PKTraj(TRange const& range, double mass, int charge) : PTTraj<KTRAJ>(range), KInter(mass,charge) {}
      virtual ~PKTraj(){}
      // override append and prepend to check mass and charge consistency
      virtual bool append(KTRAJ const& newpiece, bool allowremove=false) override {
	if(fabs(newpiece.mass()-mass())>1e-6 || newpiece.charge() != charge()) throw std::invalid_argument("Invalid particle parameters");
	return PTTraj<KTRAJ>::append(newpiece,allowremove);
      }
      virtual bool prepend(KTRAJ const& newpiece, bool allowremove=false) override {
	if(newpiece.mass() != mass() || newpiece.charge() != charge()) throw std::invalid_argument("Invalid particle parameters");
	return PTTraj<KTRAJ>::prepend(newpiece,allowremove);
      }
      // base class overrides; these just rely on the PKTraj to find the appropriate piece
      virtual void dirVector(MDir mdir,double time,Vec3& unit) const override {
	PTTraj<KTRAJ>::nearestPiece(time).dirVector(mdir,time,unit);
      }
      virtual void momentum(double time,Mom4& mom) const override {
	PTTraj<KTRAJ>::nearestPiece(time).momentum(time,mom);
      }
      double momentum(double time) const override { return PTTraj<KTRAJ>::nearestPiece(time).momentum(time); }
      double energy(double time) const override { return PTTraj<KTRAJ>::nearestPiece(time).energy(time); }
      void rangeInTolerance(TRange& range, BField const& bfield, double tol) const override {
	PTTraj<KTRAJ>::nearestPiece(range.low()).rangeInTolerance(range,bfield,tol); }
  };
}
#endif

