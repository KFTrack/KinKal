#ifndef KinKal_TData_hh
#define KinKal_TData_hh
//
//  Data object describing fit parameters or weights
//  templated on the parameter vector dimension
//  used as part of the kinematic kalman fit
//
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include <stdexcept>

namespace KinKal {
  template <size_t DDIM> class TData {
    public:
      enum Status {unknown=-1, valid=0, invalid}; // not clear if this is the right place for this FIXME!
      // define the parameter types
      typedef ROOT::Math::SVector<double,DDIM> DVEC; // data vector
      typedef ROOT::Math::SMatrix<double,DDIM,DDIM,ROOT::Math::MatRepSym<double,DDIM> > DMAT;  // associated matrix
      Status status() const { return status_; }
      bool matrixOK() const { return status_ == valid; }
      // set status
      void setStatus(Status status) { status_ = status; }
      // construct from vector and matrix
      TData(DVEC const& vec, DMAT const& mat = ROOT::Math::SVector<double,DDIM>::SMatrixIdentity()) : vec_(vec), mat_(mat), status_(valid) {}// assume valid?  FIXME!
      TData() : status_(invalid) {}
      // copy with optional inversion
      TData(TData const& tdata, bool inv) : TData(tdata) { if (inv) invert(); }
      // accessors
      DVEC const& vec() const { return vec_; }
      DMAT const& mat() const { return mat_; }
      DVEC& vec() { return vec_; }
      DMAT& mat() { return mat_; }
      // scale the matrix
      void scale(double sfac) { mat_ *= sfac; }
      // inversion changes from params <-> weight. 
      // Invert in-place, overriding status
      void invert() {
	// first invert the matrix
	if(mat_.Invert()){
	  vec_ = mat_*vec_;
	  setStatus(valid);
	} else {
	  setStatus(invalid);
	  throw std::runtime_error("Inversion failure");
	}
	// check
	if(isnan(mat_(0,0)))throw std::runtime_error("Inversion failure");
      }
     // append
      TData & operator -= (TData const& other) {
	vec_ -= other.vec();
	mat_ -= other.mat();
	//  assume if one is OK they both are (??) FIXME! 
	if(status() == valid || other.status() == valid)status_=valid;
	return *this;
      }
      TData & operator += (TData const& other) {
	vec_ += other.vec();
	mat_ += other.mat();
	//  assume if one is OK they both are (??) FIXME! 
	if(status() == valid || other.status() == valid)status_=valid;
	return *this;
      }
    private:
      DVEC vec_; // parameters
      DMAT mat_; // parameter covariance
      Status status_; // matrix status
      // allow PData and WData to access
  };
}
#endif
