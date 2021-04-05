#include "SVmodelLogLike.h"

namespace SVmodelLogLike {
    const double resample_freq = 0.5; // resampling frequency
}

using namespace std;
using namespace SVmodelLogLike;

//' Implementing PMMH for the toy SV model.
//'
//' The main cpp-function powering the R wrapper.
//'
//' @param measurements arma::vec providing the measurements (or y or measurements)
//' @param lNumber number of particles
//' @param starting_vals arma::vec giving the three starting values for phi_x,
//'   sigma_x and beta_y
//' @return double; value of the log-likelihood estimated via the BPF
//'
//' @export
// [[Rcpp::export]]
double bpf_loglike_sv(arma::vec measurements,
                      unsigned long lNumber,
                      arma::vec starting_vals) {
    // Some other housekeeping
    // General:
    y_returns = measurements;
    long lIterates = y_returns.n_rows;
    // Set starting values:
    theta_prop.phi_x   = starting_vals(0);
    theta_prop.sigma_x = starting_vals(1);
    theta_prop.beta_y  = starting_vals(2);

    double loglike;
    try {
        //Initialize and run the sampler
        my_move_loglike = new SVmodelPMMH_LogLike;
        smc::sampler<double,smc::nullParams> Sampler(lNumber,
                                                     HistoryType::NONE,
                                                     my_move_loglike);
        // Getting a particle filtering estimate of the log likelihood.
        Sampler.SetResampleParams(ResampleType::MULTINOMIAL, resample_freq);
        Sampler.Initialise();
        Sampler.IterateUntil(lIterates-1);
        loglike = Sampler.GetLogNCPath();

        loglike -= (lIterates - 1) * log(lNumber);
        
        delete my_move_loglike;
        return(loglike);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;
    }
    return 0.0; // to provide a return
}