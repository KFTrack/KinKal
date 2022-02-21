#ifndef KinKal_StrawXingConfig_hh
#define KinKal_StrawXingConfig_hh
namespace KinKal {
  // simple struct to hold crossing calculation configuration parameters
  struct StrawXingConfig {
    double minsigdoca_; // minimum doca sigma to integrate
    double nsig_; // number of sigma past wall to consider 'inside' the straw
    StrawXingConfig(double ddmax, double nsig) : minsigdoca_(ddmax), nsig_(nsig) {}
  };
}
#endif

