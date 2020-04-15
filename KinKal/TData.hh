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
      // define the parameter types
      typedef ROOT::Math::SVector<double,DDIM> DVEC; // data vector
      typedef ROOT::Math::SMatrix<double,DDIM,DDIM,ROOT::Math::MatRepSym<double,DDIM> > DMAT;  // associated matrix
      // construct from vector and matrix
      TData(DVEC const& vec, DMAT const& mat) : vec_(vec), mat_(mat) {}
      TData(DVEC const& vec) : vec_(vec)  {}
      TData() {}
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
      // Invert in-place
      void invert() {
	// first invert the matrix
	if(mat_.Invert()){
	  vec_ = mat_*vec_;
	} else {
	  throw std::runtime_error("Inversion failure");
	}
	// check
	if(isnan(mat_(0,0)))throw std::runtime_error("Inversion failure");
      }
     // append
      TData & operator -= (TData const& other) {
	vec_ -= other.vec();
	mat_ -= other.mat();
	return *this;
      }
      TData & operator += (TData const& other) {
	vec_ += other.vec();
	mat_ += other.mat();
	return *this;
      }
    private:
      DVEC vec_; // parameters
      DMAT mat_; // parameter covariance
  };
}
#endif
