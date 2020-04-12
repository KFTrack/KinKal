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
      typedef PTTraj<KTRAJ> PTTRAJ;
      typedef typename KTRAJ::PDER PDER; // forward derivative type from the 0th traj parameters
      constexpr static size_t NParams() { return KTRAJ::NParams(); } // forward the parameter space dimension
      // base class implementation
      // construct from an initial piece, which also provides kinematic information
      PKTraj(KTRAJ const& piece) : PTTRAJ(piece), KInter(piece.mass(),piece.charge()) {}
      PKTraj(TRange const& range, double mass, int charge) : PTTRAJ(range), KInter(mass,charge) {}
      //  append and prepend to check mass and charge consistency
      void append(KTRAJ const& newpiece, bool allowremove=false)  {
	if(fabs(newpiece.mass()-mass())>1e-6 || newpiece.charge() != charge()) throw std::invalid_argument("Invalid particle parameters");
	PTTRAJ::append(newpiece,allowremove);
      }
      void prepend(KTRAJ const& newpiece, bool allowremove=false)  {
	if(newpiece.mass() != mass() || newpiece.charge() != charge()) throw std::invalid_argument("Invalid particle parameters");
	PTTRAJ::prepend(newpiece,allowremove);
      }
      // base class s; these just rely on the PKTraj to find the appropriate piece
      void dirVector(MDir mdir,float time,Vec3& unit) const  {
	PTTRAJ::nearestPiece(time).dirVector(mdir,time,unit);
      }
      void momentum(float time,Mom4& mom) const  {
	PTTRAJ::nearestPiece(time).momentum(time,mom);
      }
      double momentum(float time) const  { return PTTRAJ::nearestPiece(time).momentum(time); }
      double momentumVar(float time) const  { return PTTRAJ::nearestPiece(time).momentumVar(time); }
      double energy(float time) const  { return PTTRAJ::nearestPiece(time).energy(time); }
      void rangeInTolerance(TRange& range, BField const& bfield, double tol) const  {
	PTTRAJ::nearestPiece(range.low()).rangeInTolerance(range,bfield,tol); }
  };
}
#endif

