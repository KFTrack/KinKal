/*

To instantiate KKTrk on KTLine we need to specialize TPoca on the pair
<KTLine, TLine>. This class takes care of that.

S Middleton 2020

*/

#include "KinKal/TPoca.hh"
#include "KinKal/KTLine.hh"
#include "KinKal/TLine.hh"
#include "KinKal/KTLine.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/POCAUtil.hh"
#include <limits>

// specializations for TPoca
using namespace std;
namespace KinKal {

//1) Specialization for KTLine:
  template<> TPoca<KTLine,TLine>::TPoca(KTLine const& ktline, TLine const& tline, TPocaHint const& hint, double precision) : TPocaBase(precision),ktraj_(&ktline), straj_(&tline) {
    // reset status
    reset();
    double htoca= ktline.ztime(tline.z0());
    double stoca=  tline.t0();
    double doca(0.0);
    double ktspeed = ktline.speed(ktline.t0());

    Vec3 ktpos = ktline.position(htoca); 
    Vec3 ktdir = ktline.direction(htoca); 

    POCAUtil *poca = new POCAUtil(ktpos,ktdir, tline.pos0(), tline.dir());
    htoca = poca->s1()/ktspeed;
    stoca = tline.t0() - poca->s2()/tline.speed(stoca) - stoca;
    double dd2 = poca->dca(); 

    if(dd2 < 0.0 ){
      status_ = TPoca::pocafailed;
    }
    doca = sqrt(dd2);
    if(isnan(doca)){
      status_ = pocafailed;
    }
    if(status_ != TPoca::pocafailed){
      status_ = TPoca::converged;
    }else{
      status_ = TPoca::unconverged;
    }

    partPoca_.SetE(htoca);
    ktline.position(partPoca_);
    sensPoca_.SetE(stoca);
    tline.position(sensPoca_);

    // sign doca by angular momentum projected onto difference vector (same as helix)
    double lsign =   tline.dir().Cross(ktline.direction(partPoca_.T())).Dot(sensPoca_.Vect()-partPoca_.Vect());
    double dsign = copysign(1.0,lsign);
    doca_ = doca*dsign;
    std::cout<<"DOCA: "<<doca<<std::endl;
    Vec3 ddir;
    ddir = delta().Vect().Unit();// direction vector along D(POCA) from traj 2 to 1 (line to ktline)
    ktline.direction(particlePoca().T());
    double d0 = ktline.d0();
    double cosPhi0 = ktline.cosPhi0();
    double sinPhi0 = ktline.sinPhi0();
    double sinT = ktline.sinTheta();
    std::cout << "In KTLine. Params set to: " <<ktline.params()<< std::endl;
    double f = CLHEP::c_light * ktline.beta() * (particlePoca().T() - ktline.t0());
    dDdP_[KTLine::phi0_] = -dsign*((d0*sinPhi0) + f*sinT*sinPhi0)*ddir.x()+(d0*cosPhi0 - f*sinT*sinPhi0)*ddir.y();
    dDdP_[KTLine::cost_] = f*ddir.z();
    dDdP_[KTLine::d0_] = -dsign*(-1*cosPhi0*ddir.x()+sinPhi0*ddir.y());
    dDdP_[KTLine::z0_] = -dsign*ddir.z();
    dTdP_[KTLine::t0_] = -1.0; 

    // propagate parameter covariance to variance on doca and toca
    docavar_ = ROOT::Math::Similarity(dDdP(),ktline.params().covariance());
    tocavar_ = ROOT::Math::Similarity(dTdP(),ktline.params().covariance());
    ddot_ = ktline.direction(particleToca()).Dot(tline.direction(sensorToca()));

  }

   // specialization between a piecewise KTLine and a line
  typedef PKTraj<KTLine> PKTLINE;
  template<> TPoca<PKTLINE,TLine>::TPoca(PKTLINE const& pktline, TLine const& tline, TPocaHint const& hint, double precision) : TPocaBase(precision), ktraj_(&pktline), straj_(&tline)  {

    static const unsigned maxiter=10; 
    unsigned niter=0;
    size_t oldindex= pktline.pieces().size(); 
    size_t index;
    index = size_t(rint(oldindex/2.0));
    status_ = converged; 
    while(status_ == converged and index != oldindex){

      KTLine const& piece = pktline.pieces()[index];
      TPoca<KTLine,TLine> tpoca(piece,tline,hint,precision);
      status_ = tpoca.status();
      if(tpoca.usable()){
	      partPoca_ = tpoca.particlePoca();
	      sensPoca_ = tpoca.sensorPoca();
	      doca_ = tpoca.doca();
	      dDdP_ = tpoca.dDdP();
	      dTdP_ = tpoca.dTdP();
	      docavar_ = tpoca.docaVar();
	      tocavar_ = tpoca.tocaVar();
	      ddot_ = tpoca.dirDot();
      }
     oldindex = index;
     index = pktline.nearestIndex(tpoca.particlePoca().T());
    }
    if(status_ == converged && niter >= maxiter) status_ = unconverged;
  }

}
