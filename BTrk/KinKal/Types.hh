#ifndef KinKal_Types_hh
#define KinKal_Types_hh
// typedefs used throughout the KinKal package
#include "Math/Cartesian3D.h"
#include "Math/Vector4D.h"
namespace KinKal {

  typedef ROOT::Math::Cartesian3D<double> Pos3; // spatial only vector
  typedef ROOT::Math::XYZTVector Pos4; // double precision spacetime vector
  typedef ROOT::Math::PxPyPzMVector Mom4; // double precision 4-momentum with explicit mass

}
#endif
