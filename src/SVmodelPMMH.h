#include "RcppSMC.h"
#include <cmath>
namespace SVmodelPMMH {

// Adding classes
// I. Define class for parameters
class parameters {
public:
  double phiX, sigmaX, betaY;
};
// II. Derived class for the proposal/moves
class SVmodelPMMHmove : public smc::moveset<double, smc::nullParams> {
public:
  void pfInitialise(double& value,
                    double& logweight,
                    smc::nullParams& param);
  void pfMove(long lTime,
              double& value,
              double& logweight,
              smc::nullParams& param);
  
  ~SVmodelPMMHmove() {};
};
double getLogPriorParam(const parameters& proposal);
parameters thetaProp;
smc::moveset<double, smc::nullParams>* myMovePMMH;

arma::vec yPMMHsimul;

//' A function to initialize a particle 
//' 
//' Used implicitly by myMove at the corresponding particle class.
//'
//' @param X a reference to the (empty) particle value
//' @param logweight a reference to the empty particle log weight
//' @param param additional algorithm parameters
//' 
//' @return no return; modify in place
void SVmodelPMMHmove::pfInitialise(double& X,
                                    double& logweight,
                                    smc::nullParams& param)
{
  // X = R::rnorm(0.0, std::pow(thetaProp.sigmaX, 0.5) / pow((1 - pow(thetaProp.phiX, 2)), 0.5)); 
  X = R::rnorm(0.0, thetaProp.sigmaX / pow((1 - pow(thetaProp.phiX, 2)), 0.5)); 
  double sd = thetaProp.betaY * exp(0.5 * X);
  logweight = R::dnorm(yPMMHsimul(0), 0.0, sd,TRUE);
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
void SVmodelPMMHmove::pfMove(long lTime,
                              double& X,
                              double& logweight,
                              smc::nullParams& param)
{
  X  = thetaProp.phiX * X;
  X += R::rnorm(0.0, thetaProp.sigmaX);
  double sd = thetaProp.betaY * exp(0.5 * X);
  logweight += R::dnorm(yPMMHsimul(lTime), 0.0, sd, TRUE);
}
}
