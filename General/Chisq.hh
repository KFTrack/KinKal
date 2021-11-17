#ifndef KinKal_Chisq_hh
#define KinKal_Chisq_hh
//
// Simple struct representing chisquared from the least-squares fit
//
#include "TMath.h"
#include <ostream>
namespace KinKal {
  class Chisq {
    public:
      // default constructor
      Chisq() : chisq_(0.0), ndof_(0) {}
      Chisq(double chisq, unsigned ndof) : chisq_(chisq), ndof_(ndof) {}
      double chisq() const { return chisq_; }
      int nDOF() const { return ndof_; }
      double chisqPerNDOF() const {
        if(ndof_ > 0)
          return chisq_/ndof_;
        else
          return -1.0;
      }
      // chisquared is an additive quantity
      Chisq& operator += (const Chisq& other) {
        chisq_ += other.chisq();
        ndof_ += other.nDOF();
        return *this;
      }
      // chisquared probability (= consistency)
      double probability() const {
        if(ndof_ > 0 && chisq_ > 0.0)
          return TMath::Prob(chisq_,ndof_);
        else
          return -1.0;}
    private:
      double chisq_; // value of chisquared
      int ndof_; // associated number of degrees of freedom
  };
  std::ostream& operator <<(std::ostream& os, Chisq const& chisq );
}
#endif
