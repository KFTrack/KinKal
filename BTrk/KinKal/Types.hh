#ifndef KinKal_Types_hh
#define KinKal_Types_hh
// typedefs used throughout the KinKal package
#include "Math/Cartesian3D.h"
#include "Math/Vector4D.h"
//#include "Math/Vector4Dfwd.h"
namespace KinKal {

  typedef ROOT::Math::Cartesian3D<double> ThreeV;
  typedef ROOT::Math::XYZTVector FourV; // double precision. can used for momentum+energy and spacetime vectors

}
#endif
