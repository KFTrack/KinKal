#ifndef KinKal_PKTraj_hh
#define KinKal_PKTraj_hh
//
//  Piecewise trajectory with kinematic information.  Templated on a KTraj.
//  used as part of the kinematic kalman fit
//
#include "KinKal/PTTraj.hh"
#include "KinKal/BField.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class PKTraj : public PTTraj<KTRAJ> {
    public:
      typedef PTTraj<KTRAJ> PTTRAJ;
      typedef typename KTRAJ::DVEC DVEC; // forward derivative type from the ktraj parameters
      constexpr static size_t NParams() { return KTRAJ::NParams(); } // forward the parameter space dimension
      // base class implementation
      // construct from an initial piece, which also provides kinematic information
      PKTraj(KTRAJ const& piece) : PTTRAJ(piece) {}
      PKTraj() : PTTRAJ() {}
      //  append and prepend to check mass and charge consistency
      void append(KTRAJ const& newpiece, bool allowremove=false)  {
	if(PTTRAJ::pieces().size() > 0){
	  if(fabs(newpiece.mass()-mass())>1e-6 || newpiece.charge() != charge()) throw std::invalid_argument("Invalid particle parameters");
	}
	PTTRAJ::append(newpiece,allowremove);
      }
      void prepend(KTRAJ const& newpiece, bool allowremove=false)  {
	if(fabs(newpiece.mass()-mass())>1e-6 || newpiece.charge() != charge()) throw std::invalid_argument("Invalid particle parameters");
	PTTRAJ::prepend(newpiece,allowremove);
      }
      Mom4 momentum(double time) const  { return PTTRAJ::nearestPiece(time).momentum(time); }
      double momentumMag(double time) const  { return PTTRAJ::nearestPiece(time).momentumMag(time); }
      double momentumVar(double time) const  { return PTTRAJ::nearestPiece(time).momentumVar(time); }
      double energy(double time) const  { return PTTRAJ::nearestPiece(time).energy(time); }
      double mass() const { return PTTRAJ::front().mass(); } // this will throw for empty
      double charge() const { return PTTRAJ::front().charge(); } // this will throw for empty 
      void rangeInTolerance(TRange& range, BField const& bfield, double tol) const  {
      // this could have a smarter implementation FIXME!
	PTTRAJ::nearestPiece(range.low()).rangeInTolerance(range,bfield,tol); }
      Vec3 const& bnom(double time) const { return PTTRAJ::nearestPiece(time).bnom(); }
  };
}
#endif

