#ifndef KinKal_FitData_hh
#define KinKal_FitData_hh
//
//  Data object describing fit parameters or weights
//  dimension is always 6 for kinematic parameterizations
//  used as part of the kinematic kalman fit
//
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include "KinKal/Vectors.hh"
#include <stdexcept>

namespace KinKal {
  class FitData {
    public:
      // construct from vector and matrix
      FitData(DVEC const& vec, DMAT const& mat) : vec_(vec), mat_(mat) {}
      FitData(DVEC const& vec) : vec_(vec)  {}
      FitData() {}
      // copy with optional inversion
      FitData(FitData const& tdata, bool inv=false) : vec_(tdata.vec_), mat_(tdata.mat_) { if (inv) invert(); }
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
   if(std::isnan(mat_(0,0)))throw std::runtime_error("Inversion failure");

      }
     // append
      FitData & operator -= (FitData const& other) {
	vec_ -= other.vec();
	mat_ -= other.mat();
	return *this;
      }
      FitData & operator += (FitData const& other) {
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
