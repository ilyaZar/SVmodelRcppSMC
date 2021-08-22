#include "SVmodelPMMH.h"

namespace SVmodelPMMH {
    const double resampleFreq = 0.5; // resampling frequency
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
//' @param initVals arma::vec giving the three starting values for phiX,
//'   sigmaX and betaY
//' @param rwMHvar standard deviations for the RW-MH proposal step
//' @param numProgressOutputs int giving the number of progress outputs i.e.
//'   if set to 10, then progress output occurs for every additional 10% of 
//'   completion
//' @return Rcpp::List containing the results: parameter samples (sigmaX, 
//'   betaY) and log-prior and log-likelihoood estimates
//'
//' @export
// [[Rcpp::export]]
Rcpp::List svModelPMMHimpl(arma::vec measurements,
                             unsigned long lNumber,
                             unsigned long lMCMCits,
                             arma::vec initVals,
                             arma::vec rwMHvar,
                             const int numProgressOutputs = 10) {
    // Initializing data containers
    arma::vec sigmaX  = arma::zeros(lMCMCits+1);
    arma::vec betaY   = arma::zeros(lMCMCits+1);
    arma::vec loglike  = arma::zeros(lMCMCits+1);
    arma::vec logprior = arma::zeros(lMCMCits+1);
    
    // Some other housekeeping
    // General:
    yPMMHsimul = measurements;
    long lIterates = yPMMHsimul.n_rows;
    // Set starting values:
    thetaProp.phiX   = initVals(0);
    thetaProp.sigmaX = initVals(1);
    thetaProp.betaY  = initVals(2);
    sigmaX(0) = thetaProp.sigmaX;
    betaY(0)  = thetaProp.betaY;
    // Set variables related to the MH-step:
    double mhRatio      = 0.0;
    double loglikeProp  = 0.0;
    double logpriorProp = 0.0;
    double rwVarX = rwMHvar(0);
    double rwVarY = rwMHvar(1);
    // Set variables related to progress monitoring of the overall procedure:
    int progressIntervallNum = round(lMCMCits/numProgressOutputs);
    double checkProg       = 0.0;
    double progressRatio   = 0.0;
    double acceptanceRate  = 0.0;
    int acceptanceRateAdd = 0;
    try {
        //Initialize and run the sampler
        myMovePMMH = new SVmodelPMMHmove;
        smc::sampler<double,smc::nullParams> Sampler(lNumber,
                                                     HistoryType::NONE,
                                                     myMovePMMH);
        
        Rcpp::NumericVector innovationsX = Rcpp::rnorm(lMCMCits, 0, rwVarX);
        Rcpp::NumericVector innovationsY = Rcpp::rnorm(lMCMCits, 0, rwVarY);
        Rcpp::NumericVector mhUnifRand  = Rcpp::runif(lMCMCits);
        
        // Getting a particle filtering estimate of the log likelihood.
        Sampler.SetResampleParams(ResampleType::MULTINOMIAL, resampleFreq);
        Sampler.Initialise();
        Sampler.IterateUntil(lIterates-1);
        loglike(0) = Sampler.GetLogNCPath();
        
        // Inverse gamma prior calculations
        logprior(0) = getLogPriorParam(thetaProp);
        Rcpp::Rcout << "Starting PMMH ..." << std::endl;
        for (unsigned int i = 1; i < lMCMCits + 1; ++i){
            // RW proposal for parameters
            thetaProp.sigmaX = sigmaX(i - 1) + innovationsX(i - 1);
            thetaProp.betaY = betaY(i - 1) + innovationsY(i - 1);
            
            // Getting a particle filtering estimate of the log likelihood.
            Sampler.Initialise();
            Sampler.IterateUntil(lIterates - 1);
            loglikeProp = Sampler.GetLogNCPath();
            
            // Inverse gamma prior
            logpriorProp = getLogPriorParam(thetaProp);
            
            mhRatio  = loglikeProp - loglike(i - 1);
            mhRatio += logpriorProp - logprior(i - 1);
            mhRatio  = exp(mhRatio);
            
            if (mhRatio > mhUnifRand(i - 1)){
                sigmaX(i)  = thetaProp.sigmaX;
                betaY(i)   = thetaProp.betaY;
                loglike(i)  = loglikeProp;
                logprior(i) = logpriorProp;
                
                acceptanceRateAdd++;
            } else {
                sigmaX(i)  = sigmaX(i -1);
                betaY(i)   = betaY(i -1);
                loglike(i)  = loglike(i -1);
                logprior(i) = logprior(i -1);
            }
            checkProg = i % progressIntervallNum;
            if (checkProg == 0) {
                Rcpp::Rcout << "#################################" << std::endl;
                progressRatio = (double(i)/lMCMCits)*100.0;
                progressRatio = round(progressRatio);
                Rcpp::Rcout << "Percentage completed: "
                            << progressRatio << "%." << std::endl;
                
                Rcpp::Rcout << "Current mean sigmaX: " <<
                    arma::mean(sigmaX.head(i)) << std::endl;
                Rcpp::Rcout << "Current mean betaY:  " <<
                    arma::mean(betaY.head(i)) << std::endl;
                
                acceptanceRate = (double(acceptanceRateAdd)/i)*100.0;
                acceptanceRate = round(acceptanceRate);
                Rcpp::Rcout << "Current MH acceptance rate: "
                            << acceptanceRate << "%." << std::endl;
            }
        }
        delete myMovePMMH;
        
        return Rcpp::List::create(Rcpp::Named("samplesSigmaX") = sigmaX,
                                  Rcpp::Named("samplesBetaY") = betaY,
                                  Rcpp::Named("loglike") = loglike,
                                  Rcpp::Named("logprior") = logprior);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;
    }
    return R_NilValue; // to provide a return
}
namespace SVmodelPMMH {
const double priorA = 0.01;
const double priorB = 0.01;
//' A function to calculate the log prior for a proposal.
//' 
//' The prior is IG(0.01,0.01).
//'
//' @param proposal a constant reference to proposed values of the parameters
//' 
//' @return the log-prior ratio (part of the PMMH-ratio) to be combined with the
//'   likelihood ratio
double getLogPriorParam(const parameters& proposal)
{
    double out = 0;
    out += 2*priorA*log(priorB)-2*lgamma(priorA);
    out += -(priorA+1)*log(proposal.sigmaX)-priorB/proposal.sigmaX;
    out += -(priorA+1)*log(proposal.betaY)-priorB/proposal.betaY;
    return(out);
}
}
