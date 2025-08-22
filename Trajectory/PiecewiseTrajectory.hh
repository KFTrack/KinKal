#ifndef KinKal_PiecewiseTrajectory_hh
#define KinKal_PiecewiseTrajectory_hh
//
//  class describing a piecewise trajectory.  Templated on a simple time-based trajectory
//  used as part of the kinematic kalman fit
//
#include "KinKal/General/CloneContext.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/MomBasis.hh"
#include "KinKal/General/TimeRange.hh"
#include <deque>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <typeinfo>

namespace KinKal {
  template <class KTRAJ> class PiecewiseTrajectory {
    public:
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using DKTRAJ = std::deque<KTRAJPTR>;
      // forward calls to the pieces
      VEC3 position3(double time) const { return nearestPiece(time).position3(time); }
      VEC3 velocity(double time) const { return nearestPiece(time).velocity(time); }
      double speed(double time) const { return nearestPiece(time).speed(time); }
      VEC3 direction(double time, MomBasis::Direction mdir=MomBasis::momdir_) const { return nearestPiece(time).direction(time,mdir); }
      VEC3 const& bnom(double time) const { return nearestPiece(time).bnom(time); }
      double t0() const;
      TimeRange range() const { if(pieces_.size() > 0) return TimeRange(pieces_.front()->range().begin(),pieces_.back()->range().end()); else return TimeRange(); }
      void setRange(TimeRange const& trange, bool trim=false);
      // construct without any content.  Any functions except append or prepend will throw in this state
      PiecewiseTrajectory() {}
      // construct from an initial piece
      PiecewiseTrajectory(KTRAJ const& piece);
      // clone op for reinstantiation
      PiecewiseTrajectory(PiecewiseTrajectory<KTRAJ> const&);
      std::shared_ptr< PiecewiseTrajectory<KTRAJ> > clone(CloneContext&) const;
      PiecewiseTrajectory<KTRAJ>& operator=(PiecewiseTrajectory<KTRAJ> const&) = default;
      // append or prepend a piece, at the time of the corresponding end of the new trajectory.  The last
      // piece will be shortened or extended as necessary to keep time contiguous.
      // Optionally allow truncate existing pieces to accomodate this piece.
      // If appending requires truncation and allowremove=false, the piece is not appended and the return code is false
      void append(KTRAJ const& newpiece, bool allowremove=false);
      void prepend(KTRAJ const& newpiece, bool allowremove=false);
      void add(KTRAJ const& newpiece, TimeDir tdir=TimeDir::forwards, bool allowremove=false);
      // Find the piece associated with a particular time
      KTRAJPTR const& nearestTraj(double time) const { return pieces_[nearestIndex(time)]; }
      KTRAJPTR const& indexTraj(size_t index) const { return pieces_[index]; }
      KTRAJ const& nearestPiece(double time) const { return *pieces_[nearestIndex(time)]; }
      KTRAJ const& piece(size_t index) const { return *pieces_[index]; }
      KTRAJ const& front() const { return *pieces_.front(); }
      KTRAJ const& back() const { return *pieces_.back(); }
      KTRAJ& front() { return *pieces_.front(); }
      KTRAJ& back() { return *pieces_.back(); }
      KTRAJPTR const& frontPtr() const { return pieces_.front(); }
      KTRAJPTR const& backPtr() const { return pieces_.back(); }
      size_t nearestIndex(double time) const;
      DKTRAJ const& pieces() const { return pieces_; }
      // test for spatial gaps
      double gap(size_t ihigh) const;
      void gaps(double& largest, size_t& ilargest, double& average) const;
      void print(std::ostream& ost, int detail) const ;
    private:
      DKTRAJ pieces_; // constituent pieces
  };

  template <class KTRAJ> void PiecewiseTrajectory<KTRAJ>::setRange(TimeRange const& trange, bool trim) {
    // trim pieces as necessary
    if(trim){
      while(pieces_.size() > 1 && trange.begin() > pieces_.front()->range().end() ) pieces_.pop_front();
      while(pieces_.size() > 1 && trange.end() < pieces_.back()->range().begin() ) pieces_.pop_back();
    } else if(trange.begin() > pieces_.front()->range().end() || trange.end() < pieces_.back()->range().begin())
      throw std::invalid_argument("Invalid Range");
    // update piece range
    pieces_.front()->setRange(TimeRange(trange.begin(),pieces_.front()->range().end()));
    pieces_.back()->setRange(TimeRange(pieces_.back()->range().begin(),trange.end()));
  }

  template <class KTRAJ> PiecewiseTrajectory<KTRAJ>::PiecewiseTrajectory(KTRAJ const& piece) : pieces_(1,std::make_shared<KTRAJ>(piece))
  {}

  template <class KTRAJ> void PiecewiseTrajectory<KTRAJ>::add(KTRAJ const& newpiece, TimeDir tdir, bool allowremove){
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

  template <class KTRAJ> void PiecewiseTrajectory<KTRAJ>::prepend(KTRAJ const& newpiece, bool allowremove) {
    // new piece can't have null range
    if(newpiece.range().null())throw std::invalid_argument("Can't prepend null range traj");
    if(pieces_.empty()){
      pieces_.push_back(std::make_shared<KTRAJ>(newpiece));
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
          double tmin = std::min(newpiece.range().begin(),pieces_.front()->range().begin());
          pieces_.front()->range() = TimeRange(newpiece.range().end(),pieces_.front()->range().end());
          pieces_.push_front(std::make_shared<KTRAJ>(newpiece));
          pieces_.front()->range() = TimeRange(tmin,pieces_.front()->range().end());
        } else {
          throw std::invalid_argument("range error");
        }
      }
    }
  }

  template <class KTRAJ> void PiecewiseTrajectory<KTRAJ>::append(KTRAJ const& newpiece, bool allowremove) {
    // new piece can't have null range
    if(newpiece.range().null())throw std::invalid_argument("Can't append null range traj");
    if(pieces_.empty()){
      pieces_.push_back(std::make_shared<KTRAJ>(newpiece));
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
          double tmax = std::max(newpiece.range().end(),pieces_.back()->range().end());
          // truncate the range of the current back to match with the start of the new piece.
          pieces_.back()->range() = TimeRange(pieces_.back()->range().begin(),newpiece.range().begin());
          pieces_.push_back(std::make_shared<KTRAJ>(newpiece));
          pieces_.back()->range() = TimeRange(pieces_.back()->range().begin(),tmax);
        } else {
          throw std::invalid_argument("range error");
        }
      }
    }
  }

  template <class KTRAJ> size_t PiecewiseTrajectory<KTRAJ>::nearestIndex(double time) const {
    size_t retval;
    if(pieces_.empty())throw std::length_error("Empty PiecewiseTrajectory!");
    if(time <= range().begin()){
      retval = 0;
    } else if(time >= range().end()){
      retval = pieces_.size()-1;
    } else {
      // scan
      retval = 0;
      while(retval < pieces_.size() && (!pieces_[retval]->range().inRange(time)) && time > pieces_[retval]->range().end()){
        retval++;
      }
      if(retval == pieces_.size())throw std::range_error("Failed PTraj range search");
    }
    return retval;
  }

  template <class KTRAJ> double PiecewiseTrajectory<KTRAJ>::gap(size_t ihigh) const {
    double retval(0.0);
    if(ihigh>0 && ihigh < pieces_.size()){
      double jtime = pieces_[ihigh]->range().begin(); // time of the junction of this piece with its preceeding piece
      VEC3 p0,p1;
      p0 = pieces_[ihigh]->position3(jtime);
      p1 = pieces_[ihigh-1]->position3(jtime);
      retval = (p1 - p0).R();
    }
    return retval;
  }

  template <class KTRAJ> void PiecewiseTrajectory<KTRAJ>::gaps(double& largest,  size_t& ilargest, double& average) const {
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

  template <class KTRAJ> double PiecewiseTrajectory<KTRAJ>::t0() const {
  // find a self-consistent t0
    if(pieces_.empty())throw std::length_error("Empty PiecewiseTrajectory!");
    double t0 = pieces_.front()->t0();
    auto index = nearestIndex(t0);
    int oldindex = -1;
    unsigned niter(0);
    while(static_cast<int>(index) != oldindex && niter < 10){
      ++niter;
      t0 = pieces_[index]->t0();
      oldindex = index;
      index = nearestIndex(t0);
    }
    return t0;
  }

  template <class KTRAJ> void PiecewiseTrajectory<KTRAJ>::print(std::ostream& ost, int detail) const {
    double maxgap, avggap;
    size_t igap;
    gaps(maxgap,igap,avggap);
    ost << "PiecewiseTrajectory of " << KTRAJ::trajName() << " with " << range()  << " pieces " << pieces().size() << " gaps max "<< maxgap << " avg " << avggap << std::endl;
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
        piece->print(ost,detail);
      }
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, PiecewiseTrajectory<KTRAJ> const& pttraj) {
    pttraj.print(ost,0);
    return ost;
  }

  // clone op for reinstantiation
  template<class KTRAJ>
  PiecewiseTrajectory<KTRAJ>::PiecewiseTrajectory(PiecewiseTrajectory<KTRAJ> const& rhs){
    for (const auto& ptr: rhs.pieces()){
      auto piece = std::make_shared<KTRAJ>(*ptr);
      pieces_.push_back(piece);
    }
  }

  template<class KTRAJ>
  std::shared_ptr< PiecewiseTrajectory<KTRAJ> > PiecewiseTrajectory<KTRAJ>::clone(CloneContext& context) const {
    auto rv = std::make_shared< PiecewiseTrajectory<KTRAJ> >();
    for (auto const& ptr : pieces_){
      auto piece = context.get(ptr);
      rv->append(*piece);
    }
    return rv;
  }
}

#endif

