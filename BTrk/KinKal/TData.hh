#ifndef kinkal_TData_hh
#define kinkal_TData_hh
//
//  Data payload to describing either parameters or weights
//  Parameters describe values and covariance of fit results and
//  describe a physical trajectory and the mathematical inverse of
//  weights.  Weights express constraints and are the mathematical
//  inverse of parameters.
//  used as part of the kinematic kalman fit
//
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include <iostream>

namespace KinKal {

// untemplated content
  class TDataBase {
    public:
      enum DataType{param=0,weight};
      enum InversionStatus{success=0,fail};
    private:
      // specify whether the object represents parameters or weights
      DataType dtype_;
      InversionStatus istat_;
    protected:
      // override default constructor
      TDataBase(DataType dtype=param) : dtype_(dtype), istat_(success) {}
      void setFail() { istat_ = fail; }
      void invert() {
	istat_ = success;
	switch (dtype_ ){
	  case param: 
	    dtype_ = weight;
	    break;
	  case weight:
	    dtype_ = param;
	    break;
	  default:
	    std::cout << "invalid data type " << dtype_ << std::endl;
	  // should throw here FIXME!
	}
      }
  };


// main class is templated on the data vector dimension
  template <size_t DDIM> class TData : public TDataBase {
    public:
      // define the parameter types
      typedef ROOT::Math::SVector<double,DDIM> DVec; // data vector
      typedef ROOT::Math::SMatrix<double,DDIM,DDIM,ROOT::Math::MatRepSym<double,DDIM> > DMat;  // associated matrix
      constexpr static size_t PDim() { return DDIM; }

      // construct from vector and matrix
      TData(DVec const& pars, DMat const& pcov) : vec_(pars), mat_(pcov) {}
      TData(DataType dtype=param) : TDataBase(dtype) {}
      // inversion changes from params <-> weight. 
      // Invert in-place, overriding status
      void invert() {
	// first invert the matrix
	if(mat_.invert()){
	  vec_ = mat_*vec_;
	  // update base class members
	  TDataBase::invert();
	} else
	  setFail();
      }
      // invert a different object
      void invert(TData const& other) { 
	*this = other;
	invert();
      }

      // override copy constructor, with optional inversion
      TData(TData const& other,bool doinvert=false) : TDataBase(other), vec_(other.vec()), mat_(other.mat()) {
	if(doinvert)invert();
      }

      // accessors
      DVec const& vec() const { return vec_; }
      DVec const& mat() const { return mat_; }
      DVec& vec() { return vec_; }
      DVec& mat() { return mat_; }
    private:
      DVec vec_; // parameter or weight vector
      DMat mat_; // covariance or weight matrix
  };

}
#endif

