#ifndef KinKal_Vectors_hh
#define KinKal_Vectors_hh
// vetor typedefs used throughout the KinKal package
#include "Math/Cartesian3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/Vector2D.h"
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include <stdexcept>
namespace KinKal {
  constexpr size_t NParams() { return 6; } // kinematic fit parameter space and phase space dimension
// Physical vectors (space + spacetime) in GenVector format
  using Vec3 = ROOT::Math::XYZVector; // spatial only vector
  using Cyl3 = ROOT::Math::Cylindrical3D<double>; // Cylindrical vector.  Context-dependent z axis definition
  using Pol3 = ROOT::Math::Polar3D<double>; // 3D polar vector
  using Vec4 = ROOT::Math::XYZTVector; // double precision spacetime vector, 4th component = time or energy
  using Mom4 = ROOT::Math::PxPyPzMVector; // 4-momentum with 4th component = mass
  using Pol2 = ROOT::Math::Polar2D<double>; // 2D polar vector.  Context-dependent z axis definition
  // Algebraic representations of spatial vectors: these are not miscible so must be manually cross-constructed from the above
  using VMAT = ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3>>;  // matrix type for spatial vector covariance
  using DPDV = ROOT::Math::SMatrix<double,NParams(),3,ROOT::Math::MatRepStd<double,NParams(),3>>; // parameter derivatives WRT space dimension type
  using DVDP = ROOT::Math::SMatrix<double,3,NParams(),ROOT::Math::MatRepStd<double,3,NParams()>>; // space dimension derivatives WRT parameter type
  using SVec3 = ROOT::Math::SVector<double,3>;
  using RMAT = ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepStd<double,3,3>>; // algebraic rotation matrix
  // purely algebraic vectors
  using DVEC = ROOT::Math::SVector<double,NParams()>; // data vector for parameters and weights
  using DMAT = ROOT::Math::SMatrix<double,NParams(),NParams(),ROOT::Math::MatRepSym<double,NParams()>>;  // associated matrix

}
#endif
