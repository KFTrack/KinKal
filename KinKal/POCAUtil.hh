#ifndef KinKal_POCAUtil_hh
#define KinKal_POCAUtil_hh

#include "KinKal/Vectors.hh"
#include "Math/Rotation3D.h"
namespace KinKal {
  class POCAUtil{

    public:
      POCAUtil( Vec3 const& p1, Vec3 const& t1,Vec3 const& p2,Vec3 const& t2,double closeToParallelCut = 1.e-8);
      ~POCAUtil();

      double dca()   const { return _dca; }
      double dca2d() const { return _dca2d; }

      Vec3 const& point1() const { return _pca1; }
      Vec3 const& point2() const { return _pca2; }

      double s1() const { return _s1; }
      double s2() const { return _s2; }

      bool closeToParallel() const { return _closeToParallel; }

    private:

      Vec3 _p1, _t1;
      Vec3 _p2, _t2;
      double _s1, _s2;
      Vec3 _pca1, _pca2;

      double _dca;
      double _dca2d;

      bool _closeToParallel;

      double _cut;
    };
}
#endif
