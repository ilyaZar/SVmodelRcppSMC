#include "RcppSMC.h"

namespace SVmodelPMMH {

    // Adding classes
    // I. Define class for parameters
    class parameters {
    public:
        double phi_x, sigma_x, beta_y;
    };
    // II. Derived class for the proposal/moves
    class SVmodelPMMH_move : public smc::moveset<double, smc::nullParams> {
    public:
        void pfInitialise(double& value,
                          double& logweight,
                          smc::nullParams& param);
        void pfMove(long lTime,
                    double& value,
                    double& logweight,
                    smc::nullParams& param);

        ~SVmodelPMMH_move() {};
    };
    double get_logprior_param(const parameters& proposal);
    parameters theta_prop;
	smc::moveset<double, smc::nullParams>* myMove;
}
