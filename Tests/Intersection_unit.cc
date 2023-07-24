//
// Test intersections with KinKal objects
//
#include "KinKal/General/ParticleState.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersect.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Annulus.hh"
#include "KinKal/Geometry/Rectangle.hh"
#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

using VEC3 = ROOT::Math::XYZVectorD;
using KinKal::Ray;
using KinKal::Cylinder;
using KinKal::Annulus;
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
  {"pcost",      required_argument, 0, 'C'  },
  {"pphi",       required_argument, 0, 'P'  },
  {"pmom",       required_argument, 0, 'M'  },
  {"zpos",         required_argument, 0, 'z'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: IntersectionTest  --slen1 --slen2 -(cylinder radius, half length, disk u, v 1/2 size, or Annulus inner/outer radius), -scost --sphi (surface axis direction)  --pmom, --pcost, --pphi  --zpos (particle momentum in MeV/c, direction angles, z position) ");
}

int main(int argc, char** argv) {

  int opt;
  int long_index =0;
  VEC3 point(0.0,0.0,0.0);
  double scost(1.0), sphi(0.0), slen1(400), slen2(1000);
  double pcost(0.5), pphi(1.0), pmom(100);
  double zpos(0.0);
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
      case 'C' : pcost = atof(optarg);
                 break;
      case 'P' : pphi = atof(optarg);
                 break;
      case 'M' : pmom = atof(optarg);
                 break;
      case 'z' : zpos = atof(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  VEC3 origin(0.0,0.0,0.0);
  VEC3 ppos(0.0,0.0,zpos);

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

  VEC3 bnom(0.0,0.0,1.0);
  TimeRange trange(0.0,tmax);
  KinKal::KinematicLine ktraj(pstate,bnom,trange);
  // intersect with various surfaces
  Cylinder cyl(axis,origin,slen1,slen2);
  std::cout << "Test " << cyl << std::endl;
  auto kc_inter = intersect(ktraj,cyl, trange, 1.0e-8);
  std::cout << "KinematicLine Cylinder Intersection status " << kc_inter.flag_ << " position " << kc_inter.pos_ << " time " << kc_inter.time_ << std::endl;

  if(kc_inter.flag_.inbounds_){
    auto iplane = cyl.tangentRectangle(kc_inter.pos_);
    std::cout << "tangent plane at intersection " << iplane << std::endl;
  }

  Disk disk(axis,udir,origin,slen1);
  std::cout << "Test " << disk << std::endl;

  auto kd_inter = intersect(ktraj,disk, trange, 1.0e-8);
  std::cout << "KinematicLine Disk Intersection status " << kd_inter.flag_ << " position " << kd_inter.pos_ << " time " << kd_inter.time_ << std::endl;

  Annulus ann(axis,udir,origin,slen1, slen2);
  std::cout << "Test " << ann << std::endl;

  auto ka_inter = intersect(ktraj,ann, trange, 1.0e-8);
  std::cout << "KinematicLine Annulus Intersection status " << ka_inter.flag_ << " position " << ka_inter.pos_ << " time " << ka_inter.time_ << std::endl;

  Rectangle rect(axis, udir, origin, slen1, slen2);
  std::cout << "Test " << rect << std::endl;

  auto kr_inter = intersect(ktraj,rect, trange, 1.0e-8);
  std::cout << "KinematicLine Rectangle Intersection status " << kr_inter.flag_ << " position " << kr_inter.pos_ << " time " << kr_inter.time_ << std::endl;

  // now try with helices
  KinKal::LoopHelix lhelix(pstate,bnom,trange);
  auto ld_inter = intersect(lhelix,disk, trange, 1.0e-8);
  std::cout << "LoopHelix Disk Intersection status " << ld_inter.flag_ << " position " << ld_inter.pos_ << " time " << ld_inter.time_ << std::endl;

  return 0;
}
