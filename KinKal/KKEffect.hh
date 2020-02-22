#ifndef KinKal_KKEffect_hh
#define KinKal_KKEffect_hh
//
// Class representing a discrete effect along the kinematic Kalman filter track fit
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the trajectory class representing the particle in this fit
//
#include <array>
#include "KinKal/Types.hh"
#include "KinKal/TTraj.hh"
#include "KinKal/WData.hh"
namespace KinKal {

  template<class KTRAJ> class KKEffect {
    public:
      enum Status{unprocessed=0,processed,reset,failed};
      // type of the data payload used in this trajectory
      typedef typename KTRAJ::PDATA PDATA; // forward the parameter type 
      typedef WData<KTRAJ::PDATA::PDim()> WDATA; // declare the associated weights type as well 
      // process this site given the adjacent site
      virtual bool process(KKEffect const& other,TDir tdir) = 0;
      // update this site for a new refernce trajectory.  This must be overriden, but the base class implementation is still useful
      virtual void update(KTRAJ const& ref) = 0;
      // time of this site 
      virtual double time() const = 0;
      // Access the processed data from this site in a particular direction; uses lazy conversion between parameter and weights space
      PDATA const& params(TDir tdir) const {
	if(!params_[tdir].matrixOK() && weights_[tdir].matrixOK()) params_[tdir].invert(weights_[tdir]);
	return params_[tdir];
      }
      WDATA const& weights(TDir tdir) const {
	if(!weights_[tdir].matrixOK() && params_[tdir].matrixOK()) weights_[tdir].invert(params_[tdir]);
	return weights_[tdir];
      }
      Status status(TDir tdir) const { return status_[tdir]; }
      KTRAJ const& referenceTraj() const { return *reftraj_; }
      bool isActive() const { return active_; }
      virtual ~KKEffect(){} 
    protected:
      KKEffect(KTRAJ const& reftraj) : reftraj_(&reftraj), status_{{unprocessed,unprocessed}}, active_(true) {}
      KTRAJ const* reftraj_; // reference trajectory for this site
      std::array<Status,2> status_; // status of processing in each direction
      bool active_; // activity flag
      mutable std::array<PDATA,2> params_; // params after processing
      mutable std::array<WDATA,2> weights_; // weights after processing
  };

}

#endif
