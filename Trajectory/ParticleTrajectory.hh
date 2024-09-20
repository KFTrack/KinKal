#ifndef KinKal_ParticleTrajectory_hh
#define KinKal_ParticleTrajectory_hh
//
//  Particle trajectory, based on a piecewise trajectory with kinematic information, templated on a simple kinetic trajectory (KTRAJ)
//  used as part of the kinematic kalman fit
//
#include "KinKal/Trajectory/PiecewiseTrajectory.hh"
#include "KinKal/General/ParticleStateEstimate.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class ParticleTrajectory : public PiecewiseTrajectory<KTRAJ> {
    public:
      using PTTRAJ = PiecewiseTrajectory<KTRAJ>;

      // base class implementation
      // construct from an initial piece, which also provides kinematic information
      ParticleTrajectory(KTRAJ const& piece) : PTTRAJ(piece) {}
      ParticleTrajectory() : PTTRAJ() {}
      //  append and prepend to check mass and charge consistency
      void append(KTRAJ const& newpiece, bool allowremove=false)  {
        if(PTTRAJ::pieces().size() > 0){
          if(fabs(newpiece.mass()-mass())>1e-6) throw std::invalid_argument("Invalid particle parameters");
        }
        PTTRAJ::append(newpiece,allowremove);
      }
      void prepend(KTRAJ const& newpiece, bool allowremove=false)  {
        if(fabs(newpiece.mass()-mass())>1e-6) throw std::invalid_argument("Invalid particle parameters");
        PTTRAJ::prepend(newpiece,allowremove);
      }
      // kinematic interface
      VEC4 position4(double time) const  { return PTTRAJ::nearestPiece(time).position4(time); }
      VEC3 momentum3(double time) const  { return PTTRAJ::nearestPiece(time).momentum3(time); }
      MOM4 momentum4(double time) const  { return PTTRAJ::nearestPiece(time).momentum4(time); }
      double momentum(double time) const  { return PTTRAJ::nearestPiece(time).momentum(time); }
      double beta(double time) const  { return PTTRAJ::nearestPiece(time).beta(); }
      double momentumVariance(double time) const  { return PTTRAJ::nearestPiece(time).momentumVariance(time); }
      double energy(double time) const  { return PTTRAJ::nearestPiece(time).energy(time); }
      double mass() const { return PTTRAJ::front().mass(); } // this will throw for empty
      double charge() const { return PTTRAJ::front().charge(); } // this will throw for empty
      VEC3 const& bnom(double time) const { return PTTRAJ::nearestPiece(time).bnom(); }
      ParticleState state(double time) const { return PTTRAJ::nearestPiece(time).state(time); }
      ParticleStateEstimate stateEstimate(double time) const { return PTTRAJ::nearestPiece(time).stateEstimate(time); }
  };
}
#endif

