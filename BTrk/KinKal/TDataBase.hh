#ifndef KinKal_TDataBase_hh
#define KinKal_TDataBase_hh
//
// untemplated base for TData class
//
#include <stdexcept>
namespace KinKal {
  class TDataBase {
    public:
      enum DataType{param=0,weight};
      enum Status{valid=0,invalid};
      bool isValid() const { return status_ == valid; }
      bool isInvalid() const { return status_ != valid; }
    private:
      // specify whether the object represents parameters or weights
      DataType dtype_;
      Status status_;
    protected:
      // override default constructor
      TDataBase(DataType dtype=param) : dtype_(dtype), status_(invalid) {}
      void setStatus(Status status) { status_ = status; }
      void invert() {
	switch (dtype_ ){
	  case param: 
	    dtype_ = weight;
	    break;
	  case weight:
	    dtype_ = param;
	    break;
	  default:
	    throw std::invalid_argument("invalid data type ");
	    setStatus(invalid);
	}
      }
  };
}
#endif
