#include "RcppSMC.h"
#include <cmath>
namespace SVmodelLogLike {

// Adding classes
// I. Define class for parameters
class parameters {
public:
  double phi_x, sigma_x, beta_y;
};
// II. Derived class for the proposal/moves
class SVmodelPMMH_LogLike : public smc::moveset<double, smc::nullParams> {
public:
  void pfInitialise(double& value,
                    double& logweight,
                    smc::nullParams& param);
  void pfMove(long lTime,
              double& value,
              double& logweight,
              smc::nullParams& param);
  
  ~SVmodelPMMH_LogLike() {};
};
parameters theta_prop;
smc::moveset<double, smc::nullParams>* my_move_loglike;

arma::vec y_returns;

//' A function to initialize a particle 
//' 
//' Used implicitly by myMove at the corresponding particle class.
//'
//' @param X a reference to the (empty) particle value
//' @param logweight a reference to the empty particle log weight
//' @param param additional algorithm parameters
//' 
//' @return no return; modify in place
void SVmodelPMMH_LogLike::pfInitialise(double& X,
                                       double& logweight,
                                       smc::nullParams& param)
{
  X = R::rnorm(0.0, theta_prop.sigma_x / pow((1 - pow(theta_prop.phi_x, 2)), 0.5));
  double sd = theta_prop.beta_y * exp(0.5 * X);
  logweight = R::dnorm(y_returns(0), 0.0, sd,TRUE);
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
void SVmodelPMMH_LogLike::pfMove(long lTime,
                              double& X,
                              double& logweight,
                              smc::nullParams& param)
{
  X  = theta_prop.phi_x * X;
  X += R::rnorm(0.0, theta_prop.sigma_x);
  double sd = theta_prop.beta_y * exp(0.5 * X);
  logweight += R::dnorm(y_returns(lTime), 0.0, sd, TRUE);
}
}
