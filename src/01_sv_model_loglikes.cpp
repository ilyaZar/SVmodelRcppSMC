#include "SVmodelLogLike.h"

namespace SVmodelLogLike {
    const double resampleFreq = 0.5; // resampling frequency
}

using namespace std;
using namespace SVmodelLogLike;

//' Implementing PMMH for the toy SV model.
//'
//' The main cpp-function powering the R wrapper.
//'
//' @param measurements arma::vec providing the measurements (or y or measurements)
//' @param lNumber number of particles
//' @param initVals arma::vec giving the three starting values for phiX,
//'   sigmaX and betaY
//' @return double; value of the log-likelihood estimated via the BPF
//'
//' @export
// [[Rcpp::export]]
double svModelBpfLogLike(arma::vec measurements,
                         unsigned long lNumber,
                         arma::vec initVals) {
    // Some other housekeeping
    // General:
    yReturn = measurements;
    long lIterates = yReturn.n_rows;
    // Set starting values:
    thetaProp.phiX   = initVals(0);
    thetaProp.sigmaX = initVals(1);
    thetaProp.betaY  = initVals(2);

    double loglike;
    try {
        //Initialize and run the sampler
        myMoveLogLike = new SVmodelLogLikeMove;
        smc::sampler<double,smc::nullParams> Sampler(lNumber,
                                                     HistoryType::NONE,
                                                     myMoveLogLike);
        // Getting a particle filtering estimate of the log likelihood.
        Sampler.SetResampleParams(ResampleType::MULTINOMIAL, resampleFreq);
        Sampler.Initialise();
        Sampler.IterateUntil(lIterates-1);
        loglike = Sampler.GetLogNCPath();

        loglike -= (lIterates - 1) * log(lNumber);
        
        delete myMoveLogLike;
        return(loglike);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;
    }
    return 0.0; // to provide a return
}
