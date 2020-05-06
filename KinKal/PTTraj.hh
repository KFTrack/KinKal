#ifndef KinKal_PTTraj_hh
#define KinKal_PTTraj_hh
//
//  TTraj subclass describing a piecewise trajectory.  Templated on a
//  simple TTraj.
//  used as part of the kinematic kalman fit
//
#include "KinKal/TDir.hh"
#include "KinKal/Vectors.hh"
#include "KinKal/LocalBasis.hh"
#include "KinKal/TRange.hh"
#include <deque>
#include <ostream>
#include <stdexcept>
#include <typeinfo>

namespace KinKal {
  template <class TTRAJ> class PTTraj {
    public:
      constexpr static size_t NParams() { return TTRAJ::NParams(); }
      typedef typename std::deque<TTRAJ> DTTRAJ;
      // forward calls to the pieces 
      void position(Vec4& pos) const {nearestPiece(pos.T()).position(pos); }
      Vec3 position(double time) const { return nearestPiece(time).position(time); }
      Vec3 velocity(double time) const { return nearestPiece(time).velocity(time); }
      double speed(double time) const { return nearestPiece(time).speed(time); }
      Vec3 direction(double time, LocalBasis::LocDir mdir=LocalBasis::momdir) const { return nearestPiece(time).direction(time,mdir); }
      TRange range() const { return TRange(pieces_.front().range().low(),pieces_.back().range().high()); }
      void setRange(TRange const& trange, bool trim=false);
// construct without any content.  Any functions except append or prepend will throw in this state
      PTTraj() {}
// construct from an initial piece
      PTTraj(TTRAJ const& piece);
// append or prepend a piece, at the time of the corresponding end of the new trajectory.  The last 
// piece will be shortened or extended as necessary to keep time contiguous.
// Optionally allow truncate existing pieces to accomodate this piece.
// If appending requires truncation and allowremove=false, the piece is not appended and the return code is false
      void append(TTRAJ const& newpiece, bool allowremove=false);
      void prepend(TTRAJ const& newpiece, bool allowremove=false);
      void add(TTRAJ const& newpiece, TDir tdir=TDir::forwards, bool allowremove=false);
// Find the piece associated with a particular time
      TTRAJ const& nearestPiece(double time) const { return pieces_[nearestIndex(time)]; }
      TTRAJ const& front() const { return pieces_.front(); }
      TTRAJ const& back() const { return pieces_.back(); }
      TTRAJ& front() { return pieces_.front(); }
      TTRAJ& back() { return pieces_.back(); }
      size_t nearestIndex(double time) const;
      DTTRAJ const& pieces() const { return pieces_; }
      // test for spatial gaps
      double gap(size_t ihigh) const;
      void gaps(double& largest, size_t& ilargest, double& average) const;
      void print(std::ostream& ost, int detail) const ;
    private:
      DTTRAJ pieces_; // constituent pieces
  };

  template <class TTRAJ> void PTTraj<TTRAJ>::setRange(TRange const& trange, bool trim) {
// trim pieces as necessary
    if(trim){
      while(pieces_.size() > 1 && trange.low() > pieces_.front().range().high() ) pieces_.pop_front();
      while(pieces_.size() > 1 && trange.high() < pieces_.back().range().low() ) pieces_.pop_back();
    } else if(trange.low() > pieces_.front().range().high() || trange.high() < pieces_.back().range().low())
      throw std::invalid_argument("Invalid Range");
    // update piece range
    pieces_.front().setRange(TRange(trange.low(),pieces_.front().range().high()));
    pieces_.back().setRange(TRange(pieces_.front().range().low(),trange.high()));
  }

  template <class TTRAJ> PTTraj<TTRAJ>::PTTraj(TTRAJ const& piece) : pieces_(1,piece)
  {}

  template <class TTRAJ> void PTTraj<TTRAJ>::add(TTRAJ const& newpiece, TDir tdir, bool allowremove){
    switch (tdir) {
      case TDir::forwards:
	append(newpiece,allowremove);
	break;
      case TDir::backwards:
	prepend(newpiece,allowremove);
	break;
      default:
	throw std::invalid_argument("Invalid direction");
    }
  }

  template <class TTRAJ> void PTTraj<TTRAJ>::prepend(TTRAJ const& newpiece, bool allowremove) {
  // new piece can't have infinite range
    if(newpiece.range().infinite())throw std::invalid_argument("Can't prepend infinite range traj");
    if(pieces_.empty()){
      pieces_.push_back(newpiece);
      //  the range
      range() = newpiece.range();
    } else {
      // if the new piece completely contains the existing pieces, overwrite or fail
      if(newpiece.range().contains(range())){
	if(allowremove)
	  *this = PTTraj(newpiece);
	else
	  throw std::invalid_argument("range overlap");
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
	  double tmin = std::min(newpiece.range().low(),pieces_.front().range().low());
	  pieces_.front().range().low() = newpiece.range().high() +TRange::tbuff_; 
	  pieces_.push_front(newpiece);
	  pieces_.front().range().low() = tmin;
	} else {
//	  throw std::invalid_argument("range error");
	}
      }
    }
  }

  template <class TTRAJ> void PTTraj<TTRAJ>::append(TTRAJ const& newpiece, bool allowremove) {
  // new piece can't have infinite range
    if(newpiece.range().infinite())throw std::invalid_argument("Can't append infinite range traj");
    if(pieces_.empty()){
      pieces_.push_back(newpiece);
      // override the range
      range() = newpiece.range();
    } else {
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
	  // update ranges and add the piece.
	  // first, make sure we don't loose range
	  double tmax = std::max(newpiece.range().high(),pieces_.back().range().high());
	  // truncate the range of the current back to match with the start of the new piece.  Leave a buffer on the upper range to prevent overlap
	  pieces_.back().range().high() = newpiece.range().low()-TRange::tbuff_;
	  pieces_.push_back(newpiece);
	  pieces_.back().range().high() = tmax;
	} else {
//	  throw std::invalid_argument("range error");
	}
      }
    }
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

  template <class TTRAJ> double PTTraj<TTRAJ>::gap(size_t ihigh) const {
    double retval(0.0);
    if(ihigh>0 && ihigh < pieces_.size()){
      double jtime = pieces_[ihigh].range().low(); // time of the junction of this piece with its preceeding piece
      Vec3 p0,p1;
      p0 = pieces_[ihigh].position(jtime);
      p1 = pieces_[ihigh-1].position(jtime);
      retval = (p1 - p0).R();
    }
    return retval;
  }

  template <class TTRAJ> void PTTraj<TTRAJ>::gaps(double& largest,  size_t& ilargest, double& average) const {
    largest = average = 0.0;
    ilargest =0;
    // loop over adjacent pairs
    for(size_t ipair=1; ipair<pieces_.size();++ipair){
      double gval = gap(ipair);
      average += gval;
      if(gval > largest){
	largest = gval;
	ilargest = ipair;
      }
    }
    if(pieces_.size() > 1)
      average /= (pieces_.size()-1);
  }

  template <class TTRAJ> void PTTraj<TTRAJ>::print(std::ostream& ost, int detail) const {
    double maxgap, avggap;
    size_t igap;
    gaps(maxgap,igap,avggap);
    ost << "PTTraj with " << range()  << " pieces " << pieces().size() << " gaps max "<< maxgap << " avg " << avggap << std::endl;
    if(detail ==1 && pieces().size() > 0){
      ost << "Front ";
      front().print(ost,detail);
      if(pieces().size() > 1){
	ost << "Back ";
	back().print(ost,detail);
      }
    } else if (detail >1){
      unsigned ipiece(0);
      for (auto const& piece : pieces_) {
	ost << "Piece " << ipiece++ << " ";
	piece.print(ost,detail);
      }
    }
  }
  
  template <class TTRAJ> std::ostream& operator <<(std::ostream& ost, PTTraj<TTRAJ> const& pttraj) {
    pttraj.print(ost,0);
    return ost;
  }

}

#endif

