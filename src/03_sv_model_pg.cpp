#include "SVmodelPG.h"

using namespace std;
using namespace SVmodelPG;

namespace SVmodelPG {
const double resampleFreq = 0.5; // resampling frequency

}
//' Implementing PG for the toy SV model.
//'
//' The main cpp-function powering the R wrapper for the PG sampler.
//'
//' @param measurements arma::vec providing the measurements
//' @param lNumber number of particles
//' @param lMCMCits number of PG iterations
//' @param initVals arma::vec giving the three starting values for phiX,
//'   sigmaX and betaY
//' @param initReferenceTrajectory initial reference trajectory values to 
//'   condition on
//' @param numProgressOutputs int giving the number of progress outputs i.e.
//'   if set to 10, then progress output occurs for every additional 10% of 
//'   completion
//' @return Rcpp::List containing the results: parameter samples sigmaX, 
//'   betaY
//'
//' @export
// [[Rcpp::export]]
Rcpp::List svModelPGimpl(arma::vec measurements,
                         unsigned long lNumber,
                         unsigned long lMCMCits,
                         arma::vec initVals,
                         const std::vector<double> & initReferenceTrajectory,
                         const int numProgressOutputs = 10) {
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
    int progressIntervallNum = round(lMCMCits/numProgressOutputs);
    double checkProg       = 0.0;
    double progressRatio   = 0.0;
    
    std::vector<double> referenceTrajectory(initReferenceTrajectory.begin(),
                                            initReferenceTrajectory.end());
    arma::vec weightsFinalIteration(lNumber, arma::fill::zeros);
    double indexNewReferenceTrajectory = 0;
    try {
        //Initialize and run for T=0 the conditional sampler
        myMovePG = new SVmodelPGmove;
        //Initialization of the reference trajectory within the class 
        //constructor
        smc::conditionalSampler<double,smc::nullParams> conditionalSampler(lNumber,
                                                                           HistoryType::AL,
                                                                           myMovePG,
                                                                           referenceTrajectory);
        conditionalSampler.SetResampleParams(ResampleType::MULTINOMIAL, 0.5);
        // conditionalSampler.SetResampleParams(ResampleType::STRATIFIED, 0.5);
        // conditionalSampler.SetResampleParams(ResampleType::SYSTEMATIC, 0.5);
        // conditionalSampler.SetResampleParams(ResampleType::RESIDUAL, 0.5);
        
        //Initial Gibbs samples
        Theta = sampleGibbsParams(y, referenceTrajectory, Theta, PriorAB);
        sigmaX(0) = Theta.sigmaX;
        betaY(0) = Theta.betaY;
        //Conditional on Gibbs output (and initial reference trajectory), 
        //run conditional SMC; initialization using void-type 'Initialise' 
        //member function
        conditionalSampler.Initialise();
        //This requires to set the reference trajectory to condition on manually
        //Only after reference trajectory is set, run the conditional SMC
        conditionalSampler.IterateUntil(lIterates - 1);
        // conditionalSampler.SetDigitsPrint(6);
        // Rcpp::Rcout << conditionalSampler << std::endl;
        // for(int i = 0; i < lIterates - 1; ++i) {
        //     conditionalSampler.Iterate();
        //     Rcpp::Rcout << conditionalSampler << std::endl;
        // }
        // for(int i = 0; i < lNumber; ++i) {
        //     // for(int t = 0; t < lIterates; ++t){
        //         tmp = conditionalSampler.GetALineSpace(i)[0];
        //         Rcpp::Rcout << tmp << " weights " << 10 << " weights 2" << 20 << " end!"<< std::endl;
        //     // }
        // }
        // double tmp = conditionalSampler.GetMcmcRepeats();
        // tmp = 0;
        
        //Obtain particle weights
        weightsFinalIteration = conditionalSampler.GetParticleWeight();
        // Rcpp::Rcout << weightsFinalIteration << std::endl;
        //Sample new reference trajectory index uniformly from particle weights
        indexNewReferenceTrajectory = Rcpp::sample(lNumber, 1, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(weightsFinalIteration)))[0] - 1;
        // indexNewReferenceTrajectory = Rcpp::RcppArmadillo::sample(arma::linspace(0, lNumber - 1, lNumber), 1, false, weightsFinalIteration)(0);
        //Update reference trajectory: ancestral line of the sampled index
        referenceTrajectory = conditionalSampler.GetALineSpace(indexNewReferenceTrajectory);
        
        Rcpp::Rcout << "Starting PG ..." << std::endl;
        for (unsigned int i = 1; i < lMCMCits + 1; ++i){

            Theta = sampleGibbsParams(y, referenceTrajectory, Theta, PriorAB);
            sigmaX(i) = Theta.sigmaX;
            betaY(i) = Theta.betaY;
            
            //Conditional on Gibbs output (and initial reference trajectory), 
            //run conditional SMC; initialization using 'Initialise' member 
            //function that takes previous reference trajectory as input 
            //(testing the overloaded version)
            conditionalSampler.Initialise(referenceTrajectory);
            //Only after reference trajectory is set, run the conditional SMC
            conditionalSampler.IterateUntil(lIterates - 1);
            //Obtain particle weights
            weightsFinalIteration = conditionalSampler.GetParticleWeight();
            //Sample new reference trajectory index uniformly from particle 
            //weights
            indexNewReferenceTrajectory = Rcpp::sample(lNumber, 1, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(weightsFinalIteration)))[0] - 1;
            // indexNewReferenceTrajectory = Rcpp::RcppArmadillo::sample(arma::linspace(0, lNumber - 1, lNumber), 1, false, weightsFinalIteration)(0);
            //Update reference trajectory: ancestral line of the sampled index
            referenceTrajectory = conditionalSampler.GetALineSpace(indexNewReferenceTrajectory);

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
//' @param y an arma::vec of data measurements
//' @param referenceTrajectoryX a std::vector<double> passing the 
//'   reference trajectory
//' @param GibbsParameters a class storing the Gibbs samples
//' @param GibbsPriors a class storing the Gibbs sampler prior settings
//'
//' @return a parameter class 
Parameters sampleGibbsParams(arma::vec y,
                             const std::vector<double> & referenceTrajectoryX, 
                             const Parameters& GibbsParameters,
                             const Priors& GibbsPriors){ 
    
    arma::vec referenceTrajectory = arma::conv_to<arma::colvec>::from(referenceTrajectoryX);
    Parameters outTheta;
    
    int T = y.n_rows;
    
    double aSigma = GibbsPriors.a + static_cast<double>(T + 1)/2.0;
    double bSigma = GibbsPriors.b + 0.5*arma::sum(arma::pow(referenceTrajectory.tail(T - 1) - GibbsParameters.phiX*referenceTrajectory.head(T - 1), 2));
    outTheta.sigmaX = std::sqrt(1.0/R::rgamma(aSigma, 1/bSigma));
        
    double aBeta = aSigma;
    double bBeta = GibbsPriors.b + 0.5 * arma::sum(arma::exp((-1)*referenceTrajectory) % arma::pow(y, 2));
    outTheta.betaY =  std::sqrt(1.0/R::rgamma(aBeta, 1/bBeta));
    
    // no update of phi in Gibbs blocks at the moment for the SV model
    outTheta.phiX = GibbsParameters.phiX; 

    return(outTheta);
}
}
