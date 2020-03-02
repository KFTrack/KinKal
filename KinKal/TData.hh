#ifndef KinKal_TData_hh
#define KinKal_TData_hh
//
//  Data object describing fit parameters or weights
//  templated on the parameter vector dimension
//  used as part of the kinematic kalman fit
//
#include "Math/SVector.h"
#include "Math/SMatrix.h"

namespace KinKal {
  template <size_t DDIM> class TData {
    public:
      enum Status {unknown=-1, valid=0, invalid}; // not clear if this is the right place for this FIXME!
      constexpr static size_t PDim() { return DDIM; }
      Status status() const { return status_; }
      bool matrixOK() const { return status_ == valid; }
    protected:
      // define the parameter types
      typedef ROOT::Math::SVector<double,DDIM> DVec; // data vector
      typedef ROOT::Math::SMatrix<double,DDIM,DDIM,ROOT::Math::MatRepSym<double,DDIM> > DMat;  // associated matrix
      // construct from vector and matrix
      TData(DVec const& vec, DMat const& mat = ROOT::Math::SVector<double,DDIM>::SMatrixIdentity()) : vec_(vec), mat_(mat), status_(valid) {}// assume valid?  FIXME!
      TData() : status_(invalid) {}
      // copy with optional inversion
      TData(TData const& tdata, bool inv) : TData(tdata) { if (inv) invert(); }
      // accessors
      DVec const& vec() const { return vec_; }
      DMat const& mat() const { return mat_; }
      DVec& vec() { return vec_; }
      DMat& mat() { return mat_; }
      // scale the matrix
      void scale(double sfac) { mat_ *= sfac; }
      // set status
      void setStatus(Status status) { status_ = status; }

      // inversion changes from params <-> weight. 
      // Invert in-place, overriding status
      void invert() {
	// first invert the matrix
	if(mat_.Invert()){
	  vec_ = mat_*vec_;
	  setStatus(valid);
	} else
	  setStatus(invalid);
      }
      // invert a different object
      void invert(TData const& other) { 
	*this = other;
	invert();
      }
    private:
      DVec vec_; // parameters
      DMat mat_; // parameter covariance
      Status status_; // matrix status
  };
}
#endif
