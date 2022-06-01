#ifndef KinKal_StrawXingConfig_hh
#define KinKal_StrawXingConfig_hh
namespace KinKal {
  // simple struct to hold crossing calculation configuration parameters
  struct StrawXingConfig {
    bool average_; // use the average effect no matter what DOCA value is
    double minsigdoca_; // minimum doca error to use non-averaged value on
    double maxdoca_; // maximum DOCA to consider this straw hit: otherwise set no path
    double maxddoca_; // maximum DOCA to use 'exact' calculation, otherwise average
    // default constructor is functional but will always use the average correction
    StrawXingConfig() : average_(true), minsigdoca_(-1.0), maxdoca_(0.0), maxddoca_(0.0) {}
    StrawXingConfig(double minsigdoca, double maxdoca, double maxddoca) : average_(false), minsigdoca_(minsigdoca),
    maxdoca_(maxdoca), maxddoca_(maxddoca){}
  };
}
#endif

