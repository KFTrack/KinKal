#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Annulus.hh"
#include "KinKal/Geometry/Rectangle.hh"
#include "KinKal/Geometry/Frustrum.hh"
#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

using VEC3 = ROOT::Math::XYZVectorD;
using KinKal::Ray;
using KinKal::Cylinder;
using KinKal::Frustrum;
using KinKal::Annulus;
using KinKal::Rectangle;
using KinKal::IntersectFlag;
static struct option long_options[] = {
  {"x",     required_argument, 0, 'x' },
  {"y",     required_argument, 0, 'y'  },
  {"z",     required_argument, 0, 'z'  },
  {"costheta",     required_argument, 0, 't'  },
  {"phi",        required_argument, 0, 'p'  },
  {"tol",     required_argument, 0, 'T'  },
  {"forwards",     required_argument, 0, 'f'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: GeometryTest --x --y --z (point) --costheta, --phi (direction angles) --tol (tolerance) --forwards \n ");
}

int main(int argc, char** argv) {

  int opt;
  int long_index =0;
  VEC3 point(0.0,0.0,-1.0);
  double cost(0.5), phi(0.0), tol(1e-8);
  bool forwards(true);
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'x' : point.SetX(atof(optarg));
                 break;
      case 'y' : point.SetY(atof(optarg));
                 break;
      case 'z' : point.SetZ(atof(optarg));
                 break;
      case 't' : cost = atof(optarg);
                 break;
      case 'p' : phi = atof(optarg);
                 break;
      case 'f' : forwards = atoi(optarg);
                 break;
      case 'T' : tol = atof(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }

  double sint = sqrt(1.0-cost*cost);
  VEC3 rdir(sint*cos(phi), sint*sin(phi), cost);
  Ray ray(rdir,point);

  VEC3 zdir(0.0,0.0,1.0);
  VEC3 origin(0.0,0.0,0.0);
  VEC3 udir(1.0,0.0,0.0);
  Annulus ann(zdir,udir,origin,1.0,2.0);
  Rectangle rect(zdir,udir,origin,1.0,2.0);
  Cylinder cyl(zdir,origin,2.0,10.0);
  Frustrum fru(zdir,origin,3,1,10);

  std::cout << "Test " << ray << std::endl;
  std::cout << "Test " << ann << std::endl;
  std::cout << "Test " << rect << std::endl;
  std::cout << "Test " << fru << std::endl;
  std::cout << "Test " << cyl << std::endl;

  if(ann.onSurface(point,tol))
    std::cout << "On Annulus "<< std::endl;
  else
    std::cout <<  "Not On Annulus "<< std::endl;
  if(rect.onSurface(point,tol))
    std::cout << "On Rectangle "<< std::endl;
  else
    std::cout << "Not On Rectangle "<< std::endl;
  if(cyl.onSurface(point,tol))
    std::cout << "On Cylinder "<< std::endl;
  else
    std::cout << "Not On Cylinder "<< std::endl;
  if(fru.onSurface(point,tol))
    std::cout << "On Frustrum "<< std::endl;
  else
    std::cout << "Not On Frustrum "<< std::endl;

  double dist;
  auto iflag = ann.intersect(ray,dist,forwards,tol);
  if(iflag.onsurface_){
    std::cout << "Annulus intersect " << iflag.inbounds_ << " at distance " << dist << " point " << ray.position(dist) << std::endl;
    if(!ann.onSurface(ray.position(dist),tol)){
      std::cout << "Intersection point not on surface " << std::endl;
      return -2;
    }
  } else
    std::cout << "No Annulus intersection" << std::endl;

  iflag = rect.intersect(ray,dist,forwards,tol);
  if(iflag.onsurface_)
    std::cout << "Rectangle intersect " << iflag.inbounds_ << " at distance " << dist << " point " << ray.position(dist) << std::endl;
  else
    std::cout << "No Rectangle intersection" << std::endl;

  iflag = cyl.intersect(ray,dist,forwards,tol);
  if(iflag.onsurface_){
    std::cout << "Cylinder intersect " << iflag.inbounds_ << " at distance " << dist << " point " << ray.position(dist) << std::endl;
    if(!cyl.onSurface(ray.position(dist),tol)){
      std::cout << "Intersection point not on surface " << std::endl;
      return -2;
    }
    auto tplane = cyl.tangentPlane(ray.position(dist));
    std::cout << "tangent " << tplane << std::endl;
  } else
    std::cout << "No Cylinder intersection" << std::endl;

  iflag = fru.intersect(ray,dist,forwards,tol);
  if(iflag.onsurface_){
    std::cout << "Frustrum intersect " << iflag.inbounds_ << " at distance " << dist << " point " << ray.position(dist) << std::endl;
    if(!fru.onSurface(ray.position(dist),tol)){
      std::cout << "Intersection point not on surface " << std::endl;
      return -2;
    }
    auto tplane = fru.tangentPlane(ray.position(dist));
    std::cout << "tangent " << tplane << std::endl;
    if(fabs(tplane.vDirection().Dot(fru.axis()) -cos(fru.halfAngle())) > tol/fru.halfLength()){
      std::cout << "Tangent angle doesn't match" << std::endl;
      return -3;
    }

  } else
    std::cout << "No Frustrum intersection" << std::endl;

  return 0;
}
