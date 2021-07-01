#include "SVmodelPG.h"

using namespace std;
using namespace SVmodelPG;

namespace SVmodelPG {
const double resampleFreq = 0.5; // resampling frequency

}



//' Implementing PG for the toy SV model.
//'
//' The main cpp-function powering the R wrapper.
//'
//' @param measurements arma::vec providing the measurements (or y or measurements)
//' @param lNumber number of particles
//' @param lMCMCits number of PG iterations
//' @param initVals arma::vec giving the three starting values for phiX,
//'   sigmaX and betaY
//' @param initReferenceTrajectory standard deviations for the RW-MH proposal step
//' @param num_progress_outputs int giving the number of progress outputs i.e.
//'   if set to 10, then progress output occurs for every additional 10% of 
//'   completion
//' @return Rcpp::List containing the results: parameter samples (sigmaX, 
//'   betaY) and log-prior and log-likelihoood estimates
//'
//' @export
// [[Rcpp::export]]
Rcpp::List svmodelPGimpl(arma::vec measurements,
                         unsigned long lNumber,
                         unsigned long lMCMCits,
                         arma::vec initVals,
                         arma::vec initReferenceTrajectory,
                         const int num_progress_outputs = 10) {
    // Initializing data containers
    arma::vec sigmaX  = arma::zeros(lMCMCits+1);
    arma::vec betaY   = arma::zeros(lMCMCits+1);
    
    // Some other housekeeping
    // General:
    y = measurements;
    long lIterates = y.n_rows;
    PriorAB.a = 0.01;
    PriorAB.b= 0.01;
    // Set starting values:
    Theta.phiX   = initVals(0);
    Theta.sigmaX = initVals(1);
    Theta.betaY  = initVals(2);
    sigmaX(0) = Theta.sigmaX;
    betaY(0)  = Theta.betaY;
    // Set variables related to progress monitoring of the overall procedure:
    int progress_intervall_num = round(lMCMCits/num_progress_outputs);
    double check_prog       = 0.0;
    double progress_ratio   = 0.0;
    try {
        //Initialize and run the sampler
        myMovePG = new SVmodelPGmove;
        smc::sampler<double,smc::nullParams> Sampler(lNumber,
                                                     HistoryType::AL,
                                                     myMovePG);
        
        arma::vec referenceTrajectoryX(lNumber, arma::fill::zeros);
        referenceTrajectoryX = initReferenceTrajectory;
        Theta = sampleGibbsParams(y, referenceTrajectoryX, Theta, PriorAB);
        Rcpp::Rcout << "aSigma ist" << Theta.sigmaX << endl;
        Rcpp::Rcout << "bSigma ist" << Theta.betaY << endl;
        sigmaX(0) = Theta.sigmaX;
        betaY(0) = Theta.betaY;
        // Get output of conditional SMC
        Sampler.SetResampleParams(ResampleType::MULTINOMIAL, lNumber+1);
        // Sampler.Initialise();
        // Sampler.IterateUntil(lIterates-1);
        
        Rcpp::Rcout << "Starting PG ..." << std::endl;
        for (unsigned int i = 1; i < lMCMCits + 1; ++i){

            Theta = sampleGibbsParams(y, referenceTrajectoryX, Theta, PriorAB);
            sigmaX(i) = Theta.sigmaX;
            betaY(i) = Theta.betaY;
            // Getting a particle filtering estimate of the log likelihood.
            // Sampler.Initialise();
            // Sampler.IterateUntil(lIterates - 1);

            check_prog = i % progress_intervall_num;
            if (check_prog == 0) {
                Rcpp::Rcout << "#################################" << std::endl;
                progress_ratio = (double(i)/lMCMCits)*100.0;
                progress_ratio = round(progress_ratio);
                Rcpp::Rcout << "Percentage completed: "
                            << progress_ratio << "%." << std::endl;

                Rcpp::Rcout << "Current mean sigmaX: " <<
                    arma::mean(sigmaX.head(i)) << std::endl;
                Rcpp::Rcout << "Current mean betaY:  " <<
                    arma::mean(betaY.head(i)) << std::endl;

            }
        }
        delete myMovePG;
        
        return Rcpp::List::create(Rcpp::Named("samplesSigmaX") = sigmaX,
                                  Rcpp::Named("samplesBetaY") = betaY);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;
    }
    return R_NilValue; // to provide a return
}
namespace SVmodelPG {
//' A function to sample model parameters from full conditional Gibbs blocks.
//'
//' The model parameters are \code{sigmaX} and \code{betaY}, passed via a 
//' parameter class. The prior is set to \code{IG(0.01,0.01)}
//'
//' @param proposal a constant reference to proposed values of the parameters
//'
//' @return a parameter class 
Parameters sampleGibbsParams(arma::vec y, arma::vec referenceTrajectory, 
                             const Parameters& GibbsParameters,
                             const Priors& GibbsPriors){ 
    Parameters outTheta;
    
    int T = y.n_rows;
    
    double aSigma = GibbsPriors.a + static_cast<double>(T + 1)/2.0;
    double bSigma = GibbsPriors.b + 0.5*arma::sum(arma::pow(referenceTrajectory.tail(T) - GibbsParameters.phiX*referenceTrajectory.head(T), 2));
    // std::cout << "referenceTraj head" <<arma::sum(arma::pow(referenceTrajectory.tail(T) - GibbsParameters.phiX*referenceTrajectory.head(T), 2))<< std::endl;
    outTheta.sigmaX = std::sqrt(1.0/R::rgamma(aSigma, 1/bSigma));
        
    double aBeta = aSigma;
    double bBeta = GibbsPriors.b + 0.5 * arma::sum(arma::exp((-1)*referenceTrajectory.tail(T)) % arma::pow(y, 2));
    outTheta.betaY =  std::sqrt(1.0/R::rgamma(aBeta, 1/bBeta));
    
    // no update of phi in Gibbs blocks at the moment for the SV model
    outTheta.phiX = GibbsParameters.phiX; 

    return(outTheta);
}

// sample_gibbs_params <- function(y, x, param, a, b) {
//     TT <- (length(y) + 1)
//     a_sigma <- a + TT/2
//     b_sigma <- b + 0.5*sum((x[2:TT] - param*x[1:(TT - 1)])^2)
//     sig_sq_x <- 1/rgamma(1, a_sigma, b_sigma)
//     a_beta <- a + TT/2
//     b_beta <- b + 0.5*sum(exp(-x[2:TT])*(y^2))
//     bet_sq_x = 1/rgamma(1, a_beta, b_beta)
// #
//     sig_sq_x <- sqrt(sig_sq_x)
//     bet_sq_x <- sqrt(bet_sq_x)
//     return(c(sig_sq_x, bet_sq_x))
// }

// double get_logprior_param(const parameters& proposal)
// {
//     double out = 0;
//     out += 2*prior_a*log(prior_b)-2*lgamma(prior_a);
//     out += -(prior_a+1)*log(proposal.sigmaX)-prior_b/proposal.sigmaX;
//     out += -(prior_a+1)*log(proposal.betaY)-prior_b/proposal.betaY;
//     return(out);
// }
}