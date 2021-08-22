#include <RcppSMC.h>
#include <cmath>


namespace ALtesting {

class parameters {
public:
  double phiX, sigmaX, betaY;
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

parameters thetaProp;
smc::moveset<double, smc::nullParams>* my_move_al;

arma::vec ySVsimul;

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
  double sd = std::pow(thetaProp.betaY, 0.5) * exp(0.5 * X);
  logweight = R::dnorm(ySVsimul(0), 0.0, sd,TRUE);
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
  X  = thetaProp.phiX * X;
  X += R::rnorm(0.0, std::pow(thetaProp.sigmaX, 0.5));
  double sd = std::pow(thetaProp.betaY, 0.5) * exp(0.5 * X);
  logweight += R::dnorm(ySVsimul(lTime), 0.0, sd, TRUE);
}
}
