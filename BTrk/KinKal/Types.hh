#ifndef KinKal_Types_hh
#define KinKal_Types_hh
// typedefs used throughout the KinKal package
#include "Math/Cartesian3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/Vector2D.h"
namespace KinKal {
  enum TDir{forwards=0,backwards}; // time direction; forwards = increasing time, backwards = decreasing time
  typedef ROOT::Math::XYZVector Vec3; // double precision spatial only vector
  typedef ROOT::Math::XYZTVector Vec4; // double precision spacetime vector, 4th component = time or energy
  typedef ROOT::Math::PxPyPzMVector Mom4; // double precision 4-momentum with 4th component = mass
  typedef ROOT::Math::Cylindrical3D<double> Cyl3; // Cylindrical vector.  Context-dependent z axis definition
  typedef ROOT::Math::Polar2DVector Pol2; // polar vector.  Context-dependent z axis definition
}
#endif
