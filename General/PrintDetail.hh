#ifndef KinKal_PrintDetail_hh
#define KinKal_PrintDetail_hh
// Define standard verbosity levels for printout
namespace KinKal {
  struct PrintDetail {
    enum printLevel{minimal=0,basic,extended,complete,extreme};
  };
}
#endif

