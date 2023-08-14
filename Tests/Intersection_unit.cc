//
// Test intersections with KinKal objects
//
#include "KinKal/General/ParticleState.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/Intersect.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Disk.hh"
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
using KinKal::Annulus;
using KinKal::Frustrum;
using KinKal::Rectangle;
using KinKal::Disk;
using KinKal::IntersectFlag;
using KinKal::ParticleState;
using KinKal::TimeRange;
static struct option long_options[] = {
  {"scost",      required_argument, 0, 'c'  },
  {"sphi",       required_argument, 0, 'p'  },
  {"slen1",      required_argument, 0, 'r'  },
  {"slen2",      required_argument, 0, 'l'  },
  {"slen3",      required_argument, 0, 'L'  },
  {"pcost",      required_argument, 0, 'C'  },
  {"pphi",       required_argument, 0, 'P'  },
  {"pmom",       required_argument, 0, 'M'  },
  {"zpos",         required_argument, 0, 'z'  },
  {"ypos",         required_argument, 0, 'y'  },
  {"xpos",         required_argument, 0, 'x'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: IntersectionTest  --slen1 --slen2 -(cylinder radius, half length, disk u, v 1/2 size, or Annulus inner/outer radius), -scost --sphi (surface axis direction)  --pmom, --pcost, --pphi  --zpos (particle momentum in MeV/c, direction angles, z position) ");
}

int main(int argc, char** argv) {

  int opt;
  int long_index =0;
  VEC3 point(0.0,0.0,0.0);
  double scost(1.0), sphi(0.0), slen1(400), slen2(1000), slen3(500);
  double pcost(0.5), pphi(1.0), pmom(100);
  VEC3 ppos(0.0,0.0,0.0);
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'c' : scost = atof(optarg);
                 break;
      case 'p' : sphi = atof(optarg);
                 break;
      case 'r' : slen1 = atof(optarg);
                 break;
      case 'l' : slen2 = atof(optarg);
                 break;
      case 'L' : slen3 = atof(optarg);
                 break;
      case 'C' : pcost = atof(optarg);
                 break;
      case 'P' : pphi = atof(optarg);
                 break;
      case 'M' : pmom = atof(optarg);
                 break;
      case 'x' : ppos.SetX(atof(optarg));
                 break;
      case 'y' : ppos.SetY(atof(optarg));
                 break;
      case 'z' : ppos.SetZ(atof(optarg));
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  VEC3 origin(0.0,0.0,0.0);

  double ssint = sqrt(1.0-scost*scost);
  VEC3 axis(ssint*cos(sphi), ssint*sin(sphi), scost);
  VEC3 udir(scost*cos(sphi), scost*sin(sphi), -ssint);

  double psint = sqrt(1.0-pcost*pcost);
  VEC3 momvec(psint*cos(pphi), psint*sin(pphi), pcost);
  momvec *= pmom;
  ParticleState pstate(ppos,momvec,0.0,0.5,-1);
//  std::cout << "Test " << pstate << std::endl;
  std::cout << "Test particle position " << ppos << " momentum " << pstate.momentum3() << std::endl;
  double speed = pstate.speed();
  double tmax = 2*sqrt(slen1*slen1 + slen2*slen2)/speed;

  double tol(1.0e-8);
  VEC3 bnom(0.0,0.0,1.0);
  TimeRange trange(0.0,tmax);
  KinKal::KinematicLine ktraj(pstate,bnom,trange);
  // intersect with various surfaces
  Cylinder cyl(axis,origin,slen1,slen2);
  std::cout << "Test " << cyl << std::endl;
  auto kc_inter = intersect(ktraj,cyl, trange, tol);
  std::cout << "KinematicLine Cylinder Intersection status " << kc_inter << std::endl;
  if(kc_inter.inbounds_){
    auto iplane = cyl.tangentPlane(kc_inter.pos_);
    auto dist = cyl.distance(kc_inter.pos_);
    std::cout << "distance " << dist  << " tangent plane at intersection " << iplane << std::endl;
    if(fabs(dist) > tol) return -1;
  }

  Frustrum fru(axis,origin,slen1,slen3,slen2);
  std::cout << "Test " << fru << std::endl;
  auto kf_inter = intersect(ktraj,fru, trange, tol);
  std::cout << "KinematicLine Frustrum Intersection status " << kf_inter << std::endl;
  if(kf_inter.inbounds_){
    auto iplane = fru.tangentPlane(kf_inter.pos_);
    auto dist = fru.distance(kf_inter.pos_);
    std::cout << "distance " << dist  << " tangent plane at intersection " << iplane << std::endl;
    if(fabs(dist) > tol) return -1;
  }

  Disk disk(axis,udir,origin,slen1);
  std::cout << "Test " << disk << std::endl;

  auto kd_inter = intersect(ktraj,disk, trange, tol);
  std::cout << "KinematicLine Disk Intersection status " << kd_inter << std::endl;

  Annulus ann(axis,udir,origin,slen1, slen2);
  std::cout << "Test " << ann << std::endl;

  auto ka_inter = intersect(ktraj,ann, trange, tol);
  std::cout << "KinematicLine Annulus Intersection status " << ka_inter << std::endl;

  Rectangle rect(axis, udir, origin, slen1, slen2);
  std::cout << "Test " << rect << std::endl;

  auto kr_inter = intersect(ktraj,rect, trange, tol);
  std::cout << "KinematicLine Rectangle Intersection status " << kr_inter << std::endl;

  // now try with helices
  KinKal::LoopHelix lhelix(pstate,bnom,trange);
  auto ld_inter = intersect(lhelix,disk, trange, tol);
  std::cout << "LoopHelix Disk Intersection status " << ld_inter << std::endl;

  auto lc_inter = intersect(lhelix,cyl, trange, tol);
  std::cout << "loophelix cylinder intersection status " << lc_inter << std::endl;

  auto lf_inter = intersect(lhelix,fru, trange, tol);
  std::cout << "loophelix frustrum intersection status " << lf_inter << std::endl;

  // test generic surface intersection
  KinKal::Surface const& psurf = static_cast<KinKal::Surface const&>(disk);
  auto ls_inter = intersect(lhelix,psurf,trange,tol);
   std::cout << "loophelix surface (plane) intersection status " << ls_inter << std::endl;
 if(ls_inter.onsurface_ != ld_inter.onsurface_ || (ls_inter.pos_-ld_inter.pos_).R() > tol){
    std::cout << "Generic plane intersection failed" << std::endl;
    return -1;
  }
  // test generic surface intersection
  KinKal::Surface const& csurf = static_cast<KinKal::Surface const&>(cyl);
  auto ls2_inter = intersect(lhelix,csurf,trange,tol);
   std::cout << "loophelix surface (cylinder) intersection status " << ls2_inter << std::endl;
 if(ls2_inter.onsurface_ != lc_inter.onsurface_ || (ls2_inter.pos_-lc_inter.pos_).R() > tol){
    std::cout << "Generic cylinder intersection failed" << std::endl;
    return -1;
  }

  return 0;
}
