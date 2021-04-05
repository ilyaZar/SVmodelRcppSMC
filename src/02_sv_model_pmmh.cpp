#include "SVmodelPMMH.h"

namespace SVmodelPMMH {
const double prior_a = 0.01;
const double prior_b = 0.01;

const double resample_freq = 0.5; // resampling frequency
arma::vec y;
}

using namespace std;
using namespace SVmodelPMMH;

//' Implementing PMMH for the toy SV model.
//'
//' The main cpp-function powering the R wrapper.
//'
//' @param measurements arma::vec providing the measurements (or y or measurements)
//' @param lNumber number of particles
//' @param lMCMCits number of PMMH iterations
//' @param starting_vals arma::vec giving the three starting values for phi_x,
//'   sigma_x and beta_y
//' @param rw_mh_var standard deviations for the RW-MH proposal step
//' @param num_progress_outputs int giving the number of progress outputs i.e.
//'   if set to 10, then progress output occurs for every additional 10% of 
//'   completion
//' @return Rcpp::List containing the results: parameter samples (sigma_x, 
//'   beta_y) and log-prior and log-likelihoood estimates
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame sv_model_pmmh_cpp(arma::vec measurements,
                                  unsigned long lNumber,
                                  unsigned long lMCMCits,
                                  arma::vec starting_vals,
                                  arma::vec rw_mh_var,
                                  const int num_progress_outputs = 10) {
    // Initializing data containers
    arma::vec sigma_x  = arma::zeros(lMCMCits+1);
    arma::vec beta_y   = arma::zeros(lMCMCits+1);
    arma::vec loglike  = arma::zeros(lMCMCits+1);
    arma::vec logprior = arma::zeros(lMCMCits+1);
    
    // Some other housekeeping
    // General:
    y = measurements;
    long lIterates = y.n_rows;
    // Set starting values:
    theta_prop.phi_x   = starting_vals(0);
    theta_prop.sigma_x = starting_vals(1);
    theta_prop.beta_y  = starting_vals(2);
    sigma_x(0) = theta_prop.sigma_x;
    beta_y(0)  = theta_prop.beta_y;
    // Set variables related to the MH-step:
    double mh_ratio      = 0.0;
    double loglike_prop  = 0.0;
    double logprior_prop = 0.0;
    double rw_var_x = rw_mh_var(0);
    double rw_var_y = rw_mh_var(1);
    // Set variables related to progress monitoring of the overall procedure:
    int progress_intervall_num = std::round(lMCMCits/num_progress_outputs);
    double check_prog       = 0.0;
    double progress_ratio   = 0.0;
    double acceptance_rate  = 0.0;
    int acceptance_rate_add = 0;
    try {
        //Initialize and run the sampler
        myMove = new SVmodelPMMH_move;
        smc::sampler<double,smc::nullParams> Sampler(lNumber,
                                                     HistoryType::NONE,
                                                     myMove);
        
        Rcpp::NumericVector innovations_x = Rcpp::rnorm(lMCMCits, 0, rw_var_x);
        Rcpp::NumericVector innovations_y = Rcpp::rnorm(lMCMCits, 0, rw_var_y);
        Rcpp::NumericVector mh_unif_rand  = Rcpp::runif(lMCMCits);
        
        // Getting a particle filtering estimate of the log likelihood.
        Sampler.SetResampleParams(ResampleType::MULTINOMIAL, resample_freq);
        Sampler.Initialise();
        Sampler.IterateUntil(lIterates-1);
        loglike(0) = Sampler.GetLogNCPath();
        
        // Inverse gamma prior calculations
        logprior(0) = get_logprior_param(theta_prop);
        Rcpp::Rcout << "Starting PMMH ..." << std::endl;
        for (unsigned int i = 1; i < lMCMCits + 1; ++i){
            // RW proposal for parameters
            theta_prop.sigma_x = sigma_x(i - 1) + innovations_x(i - 1);
            theta_prop.beta_y = beta_y(i - 1) + innovations_y(i - 1);
            
            // Getting a particle filtering estimate of the log likelihood.
            Sampler.Initialise();
            Sampler.IterateUntil(lIterates - 1);
            loglike_prop = Sampler.GetLogNCPath();
            
            // Inverse gamma prior
            logprior_prop = get_logprior_param(theta_prop);
            
            // Rcpp::Rcout << logprior_prop << std::endl;
            mh_ratio  = loglike_prop - loglike(i - 1);
            mh_ratio += logprior_prop - logprior(i - 1);
            mh_ratio  = exp(mh_ratio);
            
            if (mh_ratio > mh_unif_rand(i - 1)){
                sigma_x(i)  = theta_prop.sigma_x;
                beta_y(i)   = theta_prop.beta_y;
                loglike(i)  = loglike_prop;
                logprior(i) = logprior_prop;
                
                acceptance_rate_add++;
            } else {
                sigma_x(i)  = sigma_x(i -1);
                beta_y(i)   = beta_y(i -1);
                loglike(i)  = loglike(i -1);
                logprior(i) = logprior(i -1);
            }
            check_prog = i % progress_intervall_num;
            if (check_prog == 0) {
                Rcpp::Rcout << "#################################" << std::endl;
                progress_ratio = (double(i)/lMCMCits)*100.0;
                progress_ratio = std::round(progress_ratio);
                Rcpp::Rcout << "Percentage completed: "
                            << progress_ratio << "%." << std::endl;
                
                Rcpp::Rcout << "Current mean sigma_x: " <<
                    arma::mean(sigma_x.head(i)) << std::endl;
                Rcpp::Rcout << "Current mean beta_y:  " <<
                    arma::mean(beta_y.head(i)) << std::endl;
                
                acceptance_rate = (double(acceptance_rate_add)/i)*100.0;
                acceptance_rate = std::round(acceptance_rate);
                Rcpp::Rcout << "Current MH acceptance rate: "
                            << acceptance_rate << "%." << std::endl;
            }
        }
        delete myMove;
        
        return Rcpp::List::create(Rcpp::Named("samples_sigma_x") = sigma_x,
                                  Rcpp::Named("samples_beta_y") = beta_y,
                                  Rcpp::Named("loglike") = loglike,
                                  Rcpp::Named("logprior") = logprior);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;
    }
    return R_NilValue; // to provide a return
}
namespace SVmodelPMMH {
//' A function to calculate the log prior for a proposal.
//' 
//' The prior is IG(0.01,0.01).
//'
//' @param proposal a constant reference to proposed values of the parameters
//' 
//' @return the log-prior ratio (part of the PMMH-ratio) to be combined with the
//'   likelihood ratio
double get_logprior_param(const parameters& proposal)
{
    double out = 0;
    out += 2*prior_a*log(prior_b)-2*lgamma(prior_a);
    out += -(prior_a+1)*log(proposal.sigma_x)-prior_b/proposal.sigma_x;
    out += -(prior_a+1)*log(proposal.beta_y)-prior_b/proposal.beta_y;
    return(out);
}
//' A function to initialize a particle 
//' 
//' Used implicitly by myMove at the corresponding particle class.
//'
//' @param X a reference to the (empty) particle value
//' @param logweight a reference to the empty particle log weight
//' @param param additional algorithm parameters
//' 
//' @return no return; modify in place
void SVmodelPMMH_move::pfInitialise(double& X,
                                    double& logweight,
                                    smc::nullParams& param)
{
    X = R::rnorm(0.0, sqrt(5.0));
    double sd = std::pow(theta_prop.beta_y, 0.5) * exp(0.5 * X);
    logweight = R::dnorm(y(0), 0.0, sd,TRUE);
}
//' The proposal function.
//' 
//' Used implicitly by myMove at the corresponding particle class.
//' 
//' @param lTim     The sampler iteration.
//' @param X            A reference to the current particle value
//' @param logweight    A reference to the current particle log weight
//' @param param additional algorithm parameters
//' 
//' @return no return; modify in place
void SVmodelPMMH_move::pfMove(long lTime,
                              double& X,
                              double& logweight,
                              smc::nullParams& param)
{
    X  = theta_prop.phi_x * X;
    X += R::rnorm(0.0, std::pow(theta_prop.sigma_x, 0.5));
    double sd = std::pow(theta_prop.beta_y, 0.5) * exp(0.5 * X);
    logweight += R::dnorm(y(lTime), 0.0, sd, TRUE);
}
}