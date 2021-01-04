#ifndef KinKal_PiecewiseTrajectory_hh
#define KinKal_PiecewiseTrajectory_hh
//
//  class describing a piecewise trajectory.  Templated on a simple time-based trajectory
//  used as part of the kinematic kalman fit
//
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/MomBasis.hh"
#include "KinKal/General/TimeRange.hh"
#include <deque>
#include <ostream>
#include <stdexcept>
#include <typeinfo>

namespace KinKal {
  template <class TTRAJ> class PiecewiseTrajectory {
    public:
      using DTTRAJ = std::deque<TTRAJ>;
      // forward calls to the pieces 
      void position3(VEC4& pos) const {nearestPiece(pos.T()).position3(pos); }
      VEC3 position3(double time) const { return nearestPiece(time).position3(time); }
      VEC3 velocity(double time) const { return nearestPiece(time).velocity(time); }
      double speed(double time) const { return nearestPiece(time).speed(time); }
      VEC3 direction(double time, MomBasis::Direction mdir=MomBasis::momdir_) const { return nearestPiece(time).direction(time,mdir); }
      TimeRange range() const { return TimeRange(pieces_.front().range().begin(),pieces_.back().range().end()); }
      void setRange(TimeRange const& trange, bool trim=false);
// construct without any content.  Any functions except append or prepend will throw in this state
      PiecewiseTrajectory() {}
// construct from an initial piece
      PiecewiseTrajectory(TTRAJ const& piece);
// append or prepend a piece, at the time of the corresponding end of the new trajectory.  The last 
// piece will be shortened or extended as necessary to keep time contiguous.
// Optionally allow truncate existing pieces to accomodate this piece.
// If appending requires truncation and allowremove=false, the piece is not appended and the return code is false
      void append(TTRAJ const& newpiece, bool allowremove=false);
      void prepend(TTRAJ const& newpiece, bool allowremove=false);
      void add(TTRAJ const& newpiece, TimeDir tdir=TimeDir::forwards, bool allowremove=false);
// Find the piece associated with a particular time
      TTRAJ const& nearestPiece(double time) const { return pieces_[nearestIndex(time)]; }
      TTRAJ const& piece(size_t index) const { return pieces_[index]; }
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

  template <class TTRAJ> void PiecewiseTrajectory<TTRAJ>::setRange(TimeRange const& trange, bool trim) {
// trim pieces as necessary
    if(trim){
      while(pieces_.size() > 1 && trange.begin() > pieces_.front().range().end() ) pieces_.pop_front();
      while(pieces_.size() > 1 && trange.end() < pieces_.back().range().begin() ) pieces_.pop_back();
    } else if(trange.begin() > pieces_.front().range().end() || trange.end() < pieces_.back().range().begin())
      throw std::invalid_argument("Invalid Range");
    // update piece range
    pieces_.front().setRange(TimeRange(trange.begin(),pieces_.front().range().end()));
    pieces_.back().setRange(TimeRange(pieces_.front().range().begin(),trange.end()));
  }

  template <class TTRAJ> PiecewiseTrajectory<TTRAJ>::PiecewiseTrajectory(TTRAJ const& piece) : pieces_(1,piece)
  {}

  template <class TTRAJ> void PiecewiseTrajectory<TTRAJ>::add(TTRAJ const& newpiece, TimeDir tdir, bool allowremove){
    switch (tdir) {
      case TimeDir::forwards:
	append(newpiece,allowremove);
	break;
      case TimeDir::backwards:
	prepend(newpiece,allowremove);
	break;
      default:
	throw std::invalid_argument("Invalid direction");
    }
  }

  template <class TTRAJ> void PiecewiseTrajectory<TTRAJ>::prepend(TTRAJ const& newpiece, bool allowremove) {
  // new piece can't have infinite range
    if(newpiece.range().infinite())throw std::invalid_argument("Can't prepend infinite range traj");
    if(pieces_.empty()){
      pieces_.push_back(newpiece);
    } else {
      // if the new piece completely contains the existing pieces, overwrite or fail
      if(newpiece.range().contains(range())){
	if(allowremove)
	  *this = PiecewiseTrajectory(newpiece);
	else
	  throw std::invalid_argument("range overlap");
      } else {
	// find the piece that needs to be modified
	size_t ipiece = nearestIndex(newpiece.range().end());
	// see if truncation is needed
	if( allowremove){
	  while(ipiece >0 ) 
	    pieces_.pop_front();
	  ipiece--;
	}
	// if we're at the start, prepend
	if(ipiece == 0){
	  // update ranges and add the piece
	  double tmin = std::min(newpiece.range().begin(),pieces_.front().range().begin());
	  pieces_.front().range().begin() = newpiece.range().end() +TimeRange::tbuff_; 
	  pieces_.push_front(newpiece);
	  pieces_.front().range().begin() = tmin;
	} else {
	  throw std::invalid_argument("range error");
	}
      }
    }
  }

  template <class TTRAJ> void PiecewiseTrajectory<TTRAJ>::append(TTRAJ const& newpiece, bool allowremove) {
  // new piece can't have infinite range
    if(newpiece.range().infinite())throw std::invalid_argument("Can't append infinite range traj");
    if(pieces_.empty()){
      pieces_.push_back(newpiece);
    } else {
      // if the new piece completely contains the existing pieces, overwrite or fail
      if(newpiece.range().begin() < range().begin()){
	if(allowremove)
	  *this = PiecewiseTrajectory(newpiece);
	else
	  throw std::invalid_argument("range overlap");
      } else {
	// find the piece that needs to be modified
	size_t ipiece = nearestIndex(newpiece.range().begin());
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
	  double tmax = std::max(newpiece.range().end(),pieces_.back().range().end());
	  // truncate the range of the current back to match with the start of the new piece.  Leave a buffer on the upper range to prevent overlap
	  pieces_.back().range().end() = newpiece.range().begin()-TimeRange::tbuff_;
	  pieces_.push_back(newpiece);
	  pieces_.back().range().end() = tmax;
	} else {
	  throw std::invalid_argument("range error");
	}
      }
    }
  }

  template <class TTRAJ> size_t PiecewiseTrajectory<TTRAJ>::nearestIndex(double time) const {
    size_t retval;
    if(pieces_.empty())throw std::length_error("Empty PiecewiseTrajectory!");
    if(time <= range().begin()){
      retval = 0;
    } else if(time >= range().end()){
      retval = pieces_.size()-1;
    } else {
      // scan
      retval = 0;
      while(retval < pieces_.size() && (!pieces_[retval].range().inRange(time)) && time > pieces_[retval].range().end()){
	retval++;
      }
      if(retval == pieces_.size())throw std::range_error("Failed PTraj range search");
    }
    return retval;
  }

  template <class TTRAJ> double PiecewiseTrajectory<TTRAJ>::gap(size_t ihigh) const {
    double retval(0.0);
    if(ihigh>0 && ihigh < pieces_.size()){
      double jtime = pieces_[ihigh].range().begin(); // time of the junction of this piece with its preceeding piece
      VEC3 p0,p1;
      p0 = pieces_[ihigh].position3(jtime);
      p1 = pieces_[ihigh-1].position3(jtime);
      retval = (p1 - p0).R();
    }
    return retval;
  }

  template <class TTRAJ> void PiecewiseTrajectory<TTRAJ>::gaps(double& largest,  size_t& ilargest, double& average) const {
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

  template <class TTRAJ> void PiecewiseTrajectory<TTRAJ>::print(std::ostream& ost, int detail) const {
    double maxgap, avggap;
    size_t igap;
    gaps(maxgap,igap,avggap);
    ost << "PiecewiseTrajectory of " << TTRAJ::trajName() << " with " << range()  << " pieces " << pieces().size() << " gaps max "<< maxgap << " avg " << avggap << std::endl;
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
  
  template <class TTRAJ> std::ostream& operator <<(std::ostream& ost, PiecewiseTrajectory<TTRAJ> const& pttraj) {
    pttraj.print(ost,0);
    return ost;
  }

}

#endif

