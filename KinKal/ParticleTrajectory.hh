#ifndef KinKal_ParticleTrajectory_hh
#define KinKal_ParticleTrajectory_hh
//
//  Particle trajectory, based on a piecewise trajectory with kinematic information, templated on a simple kinetic trajectory (KTRAJ)
//  used as part of the kinematic kalman fit
//
#include "KinKal/PieceTrajectory.hh"
#include "KinKal/ParticleState.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class ParticleTrajectory : public PieceTrajectory<KTRAJ> {
    public:
      using PTTRAJ = PieceTrajectory<KTRAJ>;
      
      // base class implementation
      // construct from an initial piece, which also provides kinematic information
      ParticleTrajectory(KTRAJ const& piece) : PTTRAJ(piece) {}
      ParticleTrajectory() : PTTRAJ() {}
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
      MOM4 momentum(double time) const  { return PTTRAJ::nearestPiece(time).momentum(time); }
      double momentumMag(double time) const  { return PTTRAJ::nearestPiece(time).momentumMag(time); }
      double momentumVar(double time) const  { return PTTRAJ::nearestPiece(time).momentumVar(time); }
      double energy(double time) const  { return PTTRAJ::nearestPiece(time).energy(time); }
      double mass() const { return PTTRAJ::front().mass(); } // this will throw for empty
      double charge() const { return PTTRAJ::front().charge(); } // this will throw for empty 
      VEC3 const& bnom(double time) const { return PTTRAJ::nearestPiece(time).bnom(); }
      ParticleState state(double time) const { return PTTRAJ::nearestPiece(time).state(time); }
      ParticleStateMeasurement measurementState(double time) const { return PTTRAJ::nearestPiece(time).measurementState(time); }
  };
}
#endif

