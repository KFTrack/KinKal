#ifndef KinKal_Data_hh
#define KinKal_Data_hh
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
  class Data {
    public:
      // construct from vector and matrix
      Data(DVEC const& vec, DMAT const& mat) : vec_(vec), mat_(mat) {}
      Data(DVEC const& vec) : vec_(vec)  {}
      Data() {}
      // copy with optional inversion
      Data(Data const& tdata, bool inv=false) : vec_(tdata.vec_), mat_(tdata.mat_) { if (inv) invert(); }
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
      Data & operator -= (Data const& other) {
	vec_ -= other.vec();
	mat_ -= other.mat();
	return *this;
      }
      Data & operator += (Data const& other) {
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
