#ifndef KinKal_PTTraj_hh
#define KinKal_PTTraj_hh
//
//  TTraj subclass describing a piecewise trajectory.  Templated on a
//  simple TTraj.
//  used as part of the kinematic kalman fit
//
#include "BTrk/KinKal/TTraj.hh"
#include "BTrk/KinKal/TrkDir.hh"
#include <deque>
#include <stdexcept>
namespace KinKal {

  template <class TT> class PTTraj : public TTraj{
    public:
      typedef typename std::deque<TT> DTT;
      // base class implementation
      virtual void position(Vec4& pos) const override;
      virtual void position(double time, Vec3& pos) const override;
      virtual void velocity(double time, Vec3& vel) const override;
      virtual void direction(double time, Vec3& dir) const override;
// must have at least 1 initial piece
      PTTraj(TT const& piece);
// append or prepend a piece, at the time of the corresponding end of the new trajectory.  The last 
// piece will be shortened or extended as necessary to keep time contiguous.
// Optionally allow truncate existing pieces to accomodate this piece.
// If appending requires truncation and allowremove=false, the return is false
      bool append(TT const& newpiece, bool allowremove=false);
      bool prepend(TT const& newpiece, bool allowremove=false);
      bool add(TT const& newpiece, TrkDir tdir=KinKal::ttpos, bool allowremove=false);
// Find the piece associated with a particular time
      TT const& nearestPiece(double time) const { return pieces_[nearestIndex(time)]; }
      size_t nearestIndex(double time) const;
      DTT const& pieces() const { return pieces_; }
      // test for spatial gaps
      void gaps(double& largest, size_t& ilargest, double& average) const;
    private:
      PTTraj() = delete; // no default/null constructor
      DTT pieces_; // constituent pieces
  };

  // implementation: just return the values from the piece
  template <class TT> void PTTraj<TT>::position(Vec4& pos) const {
    return nearestPiece(pos.T()).position(pos);
  }
  template <class TT> void PTTraj<TT>::position(double time, Vec3& pos) const {
    return nearestPiece(time).position(time,pos);
  }
  template <class TT> void PTTraj<TT>::velocity(double time, Vec3& vel) const {
    return nearestPiece(time).position(time,vel);
  }
  template <class TT> void PTTraj<TT>::direction(double time, Vec3& dir) const {
    return nearestPiece(time).position(time,dir);
  }

  template <class TT> PTTraj<TT>::PTTraj(TT const& piece) : TTraj(piece.range()),pieces_(1,piece)
  {}

  template <class TT> bool PTTraj<TT>::add(TT const& newpiece, TrkDir tdir, bool allowremove){
    bool retval(false);
    switch (tdir) {
      case KinKal::ttpos:
	retval = append(newpiece,allowremove);
	break;
      case KinKal::ttneg:
	retval = prepend(newpiece,allowremove);
	break;
      default:
	throw std::invalid_argument("Invalid direction");
    }
    return  retval;
  }

  template <class TT> bool PTTraj<TT>::prepend(TT const& newpiece, bool allowremove) {
    if(pieces_.empty())throw std::length_error("Empty PTTraj!"); 
    bool retval(allowremove);// if we allow removal this function will always succeed
    // if the new piece completely contains the existing pieces, overwrite or fail
    if(newpiece.range().contains(range())){
      if(allowremove)
	*this = PTTraj(newpiece);
    } else {
      // find the piece that needs to be modified
      size_t ipiece = nearestPiece(newpiece.range.high());
      // see if truncation is needed
      if( allowremove){
	while(ipiece >0 ) 
	  pieces_.pop_front();
	ipiece--;
      }
      // if we're at the start, prepend
      if(ipiece == 0){
	// update ranges and add the piece
	pieces_.front().range().low() = newpiece.range().high();
	range().low() = newpiece.range().low();
	pieces_.push_front(newpiece);
	// subtract a small buffer to prevent overlaps
	pieces_.front().range().high() -= TRange::tbuff_;
	retval = true;
      }
    }
    return retval;
  }

  template <class TT> bool PTTraj<TT>::append(TT const& newpiece, bool allowremove) {
    if(pieces_.empty())throw std::length_error("Empty PTTraj!"); 
    bool retval(allowremove);// if we allow removal this function will always succeed
    // if the new piece completely contains the existing pieces, overwrite or fail
    if(newpiece.range().contains(range())){
      if(allowremove)
	*this = PTTraj(newpiece);
    } else {
      // find the piece that needs to be modified
      size_t ipiece = nearestIndex(newpiece.range().low());
      // see if truncation is needed
      if( allowremove){
	while(ipiece < pieces_.size()-1) {
	  pieces_.pop_back();
	}
      }
      // if we're at the end, append
      if(ipiece == pieces_.size()-1){
	// update ranges and add the piece.  Leave a buffer on the upper range to prevent overlap
	pieces_.back().range().high() = newpiece.range().low()-TRange::tbuff_;
	range().high() = newpiece.range().high();
	pieces_.push_back(newpiece);
	retval = true;
      }
    }
    return retval;
  }

  template <class TT> size_t PTTraj<TT>::nearestIndex(double time) const {
    size_t retval;
    if(pieces_.empty())throw std::length_error("Empty PTTraj!");
    if(time <= range().low()){
      retval = 0;
    } else if(time >= range().high()){
      retval = pieces_.size()-1;
    } else {
      // scan
      retval = 0;
      while(retval < pieces_.size() && !pieces_[retval].range().inRange(time)){
	retval++;
      }
    }
    return retval;
  }

  template <class TT> void PTTraj<TT>::gaps(double& largest,  size_t& ilargest, double& average) const {
    largest = average = 0.0;
    ilargest =0;
    // loop over adjacent pairs
    for(size_t ipair=1; ipair<pieces_.size();++ipair){
      double jtime = pieces_[ipair].range().low(); // time of the junction of this piece with its preceeding piece
      Vec3 p0,p1;
      pieces_[ipair].position(jtime,p0);
      pieces_[ipair-1].position(jtime,p1);
      double diff = sqrt((p1 - p0).Mag2());
      average += diff;
      if(diff > largest){
	largest = diff;
	ilargest = ipair;
      }
    }
    if(pieces_.size() > 1)
      average /= (pieces_.size()-1);
  }
}

#endif

