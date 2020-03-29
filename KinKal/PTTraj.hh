#ifndef KinKal_PTTraj_hh
#define KinKal_PTTraj_hh
//
//  TTraj subclass describing a piecewise trajectory.  Templated on a
//  simple TTraj.
//  used as part of the kinematic kalman fit
//
#include "KinKal/TTraj.hh"
#include "KinKal/TDir.hh"
#include <deque>
#include <ostream>
#include <stdexcept>
#include <typeinfo>

namespace KinKal {
  template <class TTRAJ> class PTTraj : public TTraj {
    public:
      constexpr static size_t NParams() { return TTRAJ::NParams(); }
      typedef typename std::deque<TTRAJ> DTTRAJ;
      // base class implementation
      virtual void position(Vec4& pos) const override;
      virtual void position(double time, Vec3& pos) const override;
      virtual void velocity(double time, Vec3& vel) const override;
      virtual double speed(double time) const override;
      virtual void direction(double time, Vec3& dir) const override;
      virtual void setRange(TRange const& trange) override;
// construct without any pieces, but specify the range
      PTTraj(TRange const& range) : TTraj(range) {}
// one initial piece
      PTTraj(TTRAJ const& piece);
      virtual ~PTTraj(){} 
// append or prepend a piece, at the time of the corresponding end of the new trajectory.  The last 
// piece will be shortened or extended as necessary to keep time contiguous.
// Optionally allow truncate existing pieces to accomodate this piece.
// If appending requires truncation and allowremove=false, the piece is not appended and the return code is false
      virtual bool append(TTRAJ const& newpiece, bool allowremove=false);
      virtual bool prepend(TTRAJ const& newpiece, bool allowremove=false);
      bool add(TTRAJ const& newpiece, TDir tdir=TDir::forwards, bool allowremove=false);
// Find the piece associated with a particular time
      TTRAJ const& nearestPiece(double time) const { return pieces_[nearestIndex(time)]; }
      TTRAJ const& front() const { return pieces_.front(); }
      TTRAJ const& back() const { return pieces_.back(); }
      size_t nearestIndex(double time) const;
      DTTRAJ const& pieces() const { return pieces_; }
      // test for spatial gaps
      void gaps(double& largest, size_t& ilargest, double& average) const;
    private:
      PTTraj() = delete; // no default/null constructor
      DTTRAJ pieces_; // constituent pieces
  };

  template <class TTRAJ> std::ostream& operator <<(std::ostream& os, PTTraj<TTRAJ> const& pttraj) {
    os << "Piecewise trajectory with " << pttraj.pieces().size() << " pieces of " << typeid(TTRAJ).name() << pttraj.range() << std::endl;
    for(auto const& piece : pttraj.pieces()) os << piece << std::endl;
//    << "Front " << pttraj.front() << std::endl
//    << "Back " << pttraj.back();
    return os;
  }

  // implementation: just return the values from the piece
  template <class TTRAJ> void PTTraj<TTRAJ>::position(Vec4& pos) const {
    nearestPiece(pos.T()).position(pos);
  }
  template <class TTRAJ> void PTTraj<TTRAJ>::position(double time, Vec3& pos) const {
    nearestPiece(time).position(time,pos);
  }
  template <class TTRAJ> void PTTraj<TTRAJ>::velocity(double time, Vec3& vel) const {
    nearestPiece(time).velocity(time,vel);
  }
  template <class TTRAJ> double PTTraj<TTRAJ>::speed(double time) const {
    return nearestPiece(time).speed(time);
  }
  template <class TTRAJ> void PTTraj<TTRAJ>::direction(double time, Vec3& dir) const {
    nearestPiece(time).direction(time,dir);
  }
  template <class TTRAJ> void PTTraj<TTRAJ>::setRange(TRange const& trange) {
    TTraj::setRange(trange);
// trim pieces as necessary
    while(pieces_.size() > 0 && trange.low() > pieces_.front().range().high() ) pieces_.pop_front();
    while(pieces_.size() > 0 && trange.high() < pieces_.back().range().low() ) pieces_.pop_back();
// update piece range
    pieces_.front().setRange(TRange(trange.low(),pieces_.front().range().high()));
    pieces_.back().setRange(TRange(pieces_.front().range().low(),trange.high()));
  }

  template <class TTRAJ> PTTraj<TTRAJ>::PTTraj(TTRAJ const& piece) : TTraj(piece.range()),pieces_(1,piece)
  {}

  template <class TTRAJ> bool PTTraj<TTRAJ>::add(TTRAJ const& newpiece, TDir tdir, bool allowremove){
    bool retval(false);
    switch (tdir) {
      case TDir::forwards:
	retval = append(newpiece,allowremove);
	break;
      case TDir::backwards:
	retval = prepend(newpiece,allowremove);
	break;
      default:
	throw std::invalid_argument("Invalid direction");
    }
    return  retval;
  }

  template <class TTRAJ> bool PTTraj<TTRAJ>::prepend(TTRAJ const& newpiece, bool allowremove) {
  // new piece can't have infinite range
    if(newpiece.range().infinite())throw std::invalid_argument("Can't prepend infinite range traj");

    bool retval(false);
    if(pieces_.empty()){
      pieces_.push_back(newpiece);
      // override the range
      pieces_.front().setRange(range());
      retval = true;
    } else {
      retval = allowremove;// if we allow removal this function will always succeed
      // if the new piece completely contains the existing pieces, overwrite or fail
      if(newpiece.range().contains(range())){
	if(allowremove)
	  *this = PTTraj(newpiece);
      } else {
	// find the piece that needs to be modified
	size_t ipiece = nearestIndex(newpiece.range().high());
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
	  range().low() = std::min(range().low(),newpiece.range().low());
	  pieces_.push_front(newpiece);
	  // subtract a small buffer to prevent overlaps
	  pieces_.front().range().high() -= TRange::tbuff_;
	  retval = true;
	}
      }
    }
    return retval;
  }

  template <class TTRAJ> bool PTTraj<TTRAJ>::append(TTRAJ const& newpiece, bool allowremove) {
  // new piece can't have infinite range
    if(newpiece.range().infinite())throw std::invalid_argument("Can't append infinite range traj");
    bool retval(false);
    if(pieces_.empty()){
      pieces_.push_back(newpiece);
      // override the range
//      pieces_.front().setRange(range());
      range() = newpiece.range();
      retval = true;
    } else {
     retval = allowremove;// if we allow removal this function will always succeed
      // if the new piece completely contains the existing pieces, overwrite or fail
      if(newpiece.range().low() < range().low()){
	if(allowremove)
	  *this = PTTraj(newpiece);
	else
	  throw std::invalid_argument("range overlap");
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
	// check ranges
	  // update ranges and add the piece.  Leave a buffer on the upper range to prevent overlap
	  pieces_.back().range().high() = newpiece.range().low()-TRange::tbuff_;
	  range().high() = std::max(range().high(),newpiece.range().high());
	  pieces_.push_back(newpiece);
	  retval = true;
	} else {
	  throw std::invalid_argument("range error");
	}
      }
    }
    return retval;
  }

  template <class TTRAJ> size_t PTTraj<TTRAJ>::nearestIndex(double time) const {
    size_t retval;
    if(pieces_.empty())throw std::length_error("Empty PTTraj!");
    if(time <= range().low()){
      retval = 0;
    } else if(time >= range().high()){
      retval = pieces_.size()-1;
    } else {
      // scan
      retval = 0;
      while(retval < pieces_.size() && !pieces_[retval].range().inRange(time) && time > pieces_[retval].range().high()){
	retval++;
      }
      if(retval == pieces_.size())throw std::range_error("Failed PTraj range search");
    }
    return retval;
  }

  template <class TTRAJ> void PTTraj<TTRAJ>::gaps(double& largest,  size_t& ilargest, double& average) const {
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

