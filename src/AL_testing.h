#include <RcppSMC.h>
#include <cmath>


namespace ALtesting {

class parameters {
public:
  double phi_x, sigma_x, beta_y;
};

class ALtracking_move: public smc::moveset<double, smc::nullParams> {
public:
  void pfInitialise(double& value,
                    double& logweight,
                    smc::nullParams& param);
  void pfMove(long lTime,
              double& value,
              double& logweight,
              smc::nullParams& param);
  
  ~ALtracking_move() {};
};

parameters theta_prop;
smc::moveset<double, smc::nullParams>* my_move_al;

arma::vec y_pmmh_simul;

//' A function to initialize a particle
//' 
//' Used implicitly by myMove at the corresponding particle class.
//'
//' @param X a reference to the (empty) particle value
//' @param logweight a reference to the empty particle log weight
//' @param param additional algorithm parameters
//' 
//' @return no return; modify in place
void ALtracking_move::pfInitialise(double& X,
                                   double& logweight,
                                   smc::nullParams& param)
{
  X = R::rnorm(0.0, sqrt(5.0));
  double sd = std::pow(theta_prop.beta_y, 0.5) * exp(0.5 * X);
  logweight = R::dnorm(y_pmmh_simul(0), 0.0, sd,TRUE);
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
void ALtracking_move::pfMove(long lTime,
                             double& X,
                             double& logweight,
                             smc::nullParams& param)
{
  X  = theta_prop.phi_x * X;
  X += R::rnorm(0.0, std::pow(theta_prop.sigma_x, 0.5));
  double sd = std::pow(theta_prop.beta_y, 0.5) * exp(0.5 * X);
  logweight += R::dnorm(y_pmmh_simul(lTime), 0.0, sd, TRUE);
}
}
