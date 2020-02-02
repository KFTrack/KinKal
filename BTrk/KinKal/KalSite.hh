#ifndef KinKal_KalSite_hh
#define KinKal_KalSite_hh
//
// Class representing a processing step of the kinematic Kalman filter
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the kinematic trajectory used in this fit
//
#include <array>
#include "BTrk/KinKal/Types.hh"
namespace KinKal {

  template<class KTRAJ> class KalSite {
    public:
      enum Status{unprocessed=0,processed,failed};
      // type of the data payload used in this trajectory
      typedef typename KTRAJ::TDATA TDATA; // forward the type definition
      // actions
      virtual bool process(KalSite const& other,TDir tdir);
      //accessors
      TDATA const& params(TDir dir) const { return params_[dir]; }
      TDATA const& weights(TDir dir) const { return weights_[dir]; }
      // non-const accessors perform inversion if necessary 
      TDATA const& params(TDir dir) {
	if (params_[dir].isInvalid() && weights_[dir].isValid()) params_[dir].inverse(weights_[dir]); 
	return params_[dir];
      }
      TDATA const& weights(TDir dir) {
	if (weights_[dir].isInvalid() && params_[dir].isValid()) weights_[dir].inverse(params_[dir]); 
	return weights_[dir];
      }
      Status status(TDir dir) const { return status_[dir]; }
      KTRAJ const& referenceTraj() const { return reftraj_; }
      bool isActive() const { return active_; }
    protected:
      KalSite(double time) : status_{unprocessed}, time_(time), active_(true) {}
      std::array<TDATA,2> params_; // cache of the parameter space representation of the fit in each direction
      std::array<TDATA,2> weights_; // cache of the weight space representation of the fit in each direction
      std::array<Status,2> status_; // status of processing in each direction
      KTRAJ reftraj_; // reference trajectory for this site
      double time_; // time of this site
      bool active_; //
  };

}

#endif
