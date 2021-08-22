#include "RcppSMC.h"
#include <cmath>
namespace SVmodelPG {

// Adding classes
// I. Define class for parameters
class Parameters {
public:
  double phiX, sigmaX, betaY;
};
class Priors {
public:
  double a, b;
};
// II. Derived class for the proposal/moves
class SVmodelPGmove : public smc::moveset<double, smc::nullParams> {
public:
  void pfInitialise(double& value,
                    double& logweight,
                    smc::nullParams& param);
  void pfMove(long lTime,
              double& value,
              double& logweight,
              smc::nullParams& param);
  void pfWeight(long lTime,
                double& value,
                double& logweight,
                smc::nullParams& param);
  
  ~SVmodelPGmove() {};
};
Parameters Theta;
Priors PriorAB;
arma::vec y;

smc::moveset<double, smc::nullParams>* myMovePG;

Parameters sampleGibbsParams(arma::vec y,
                             const std::vector<double> & referenceTrajectoryX, 
                             const Parameters& currentParameterValues,
                             const Priors& GibbsPriors); 

//' A function to initialize a particle 
//' 
//' Used implicitly by myMove at the corresponding particle class.
//'
//' @param X a reference to the (empty) particle value
//' @param logweight a reference to the empty particle log weight
//' @param param additional algorithm parameters
//' 
//' @return no return; modify in place
void SVmodelPGmove::pfInitialise(double& X,
                                 double& logweight,
                                 smc::nullParams& param)
{
  X = R::rnorm(0.0, Theta.sigmaX / pow((1 - pow(Theta.phiX, 2)), 0.5)); 
  double sd = std::pow(Theta.betaY, 0.5) * exp(0.5 * X);
  logweight = R::dnorm(y(0), 0.0, sd,TRUE);
}
//' The proposal function.
//' 
//' Used implicitly by myMove at the corresponding particle class.
//' 
//' @param lTime     The sampler iteration.
//' @param X            A reference to the current particle value
//' @param logweight    A reference to the current particle log weight
//' @param param additional algorithm parameters
//' 
//' @return no return; modify in place
void SVmodelPGmove::pfMove(long lTime,
                           double& X,
                           double& logweight,
                           smc::nullParams& param)
{
  X  = Theta.phiX * X;
  X += R::rnorm(0.0, Theta.sigmaX);
  double sd = Theta.betaY * exp(0.5 * X);
  logweight += R::dnorm(y(lTime), 0.0, sd, TRUE);
}
//' Weighting function of the conditional/reference particle coordinate.
//' 
//' Used implicitly by myMove at the corresponding (derived) conditional SMC
//' class.
//' 
//' @param lTime     The sampler iteration.
//' @param X            A reference to the current particle value
//' @param logweight    A reference to the current particle log weight
//' @param param additional algorithm parameters
//' 
//' @return no return; modify in place
void SVmodelPGmove::pfWeight(long lTime,
                             double& X,
                             double& logweight,
                             smc::nullParams& param)
{
  double sd = Theta.betaY * exp(0.5 * X);
  logweight += R::dnorm(y(lTime), 0.0, sd, TRUE);
}
}
