#include "AL_testing.h"

using namespace std;
using namespace ALtesting;

//' Implementing PMMH for the toy SV model.
//'
//' The main cpp-function powering the R wrapper.
//'
//' @param measurements arma::vec providing the measurements (or y or
//' measurements)
//' @param lNumber number of particles
//' @param initVals arma::vec giving the three starting values for phiX,
//'   sigmaX and betaY
//' @param resampleFreq  frequency at which resampling is performed; if
//'   negative, then resampling is never performed; if between [0,1), then
//'   resampling is performed when  the ESS falls below that proportion of the
//'   number of particles and when it is greater than or equal to 1, resampling
//'   is carried out when the ESS falls below that value (note: if this
//'   parameter is larger than the total number of particles, then resampling
//'   will always be performed!)
//' @printALinfo boolean; if FALSE, no information about the ancestral lines is
//'   printed to the screen and just the plot is returned
//' @return Rcpp::List containing the results: parameter samples (sigmaX, 
//'   betaY) and log-prior and log-likelihoood estimates
//'
//' @export
// [[Rcpp::export]]
Rcpp::List svModelALtrackingImp(arma::vec measurements,
                                arma::vec initVals,
                                unsigned long lNumber,
                                const double resampleFreq = 0.5,
                                bool printALinfo = false) {
  // Set variables related to progress monitoring of the overall procedure:
  ySVsimul = measurements;
  long lIterates = ySVsimul.n_rows;
  // Set starting values:
  thetaProp.phiX   = initVals(0);
  thetaProp.sigmaX = initVals(1);
  thetaProp.betaY  = initVals(2);
  try {
    //Initialize and run the sampler
    my_move_al = new ALtracking_move;
    smc::sampler<double,smc::nullParams> Sampler(lNumber,
                                                 HistoryType::AL,
                                                 my_move_al);
    Sampler.SetResampleParams(ResampleType::MULTINOMIAL, resampleFreq);
    Sampler.Initialise();
    // Sampler.IterateUntil(lIterates - 1);
    for (unsigned int i = 0; i < lIterates; ++i) {
      if (printALinfo) {
        Rcpp::Rcout << i + 1 << " of " << lIterates << std::endl;
        smc::population<double> current_population = Sampler.GetHistoryPopulation(i);
        std::vector<double> current_vals = current_population.GetValue();
        arma::Col<unsigned int> current_population_al = Sampler.GetuRSIndices();
        for(int j = 0; j < static_cast<int>(current_vals.size()); j++)
        {
          Rcpp::Rcout << "Particle no: " << j + 1 << ". Ancestor index "
                      << current_population_al[j] <<  " with value " << current_vals[j]
                      << std::endl;
        }
      }
      if(i != lIterates - 1) Sampler.Iterate();
    }
    arma::Mat<unsigned int> ancestral_lines_indeces(lNumber, lIterates);
    arma::Mat<double> ancestral_lines_values(lNumber, lIterates);
    for(unsigned long i = 0; i < lNumber; ++i) {
      arma::Col<unsigned int> test_ancestral_line_index = Sampler.GetALineInd(i);
      ancestral_lines_indeces.row(i) = test_ancestral_line_index.t();
      arma::Col<double> test_ancestral_line_value = Sampler.GetALineSpace(i);
      ancestral_lines_values.row(i) = test_ancestral_line_value.t();
      if(printALinfo) {
        Rcpp::Rcout << "Particle no: " << i + 1 
                    << " ancestor line: re-mapped index and value." 
                    << std::endl;
        Rcpp::Rcout << test_ancestral_line_index.t() << std::endl;
        Rcpp::Rcout << test_ancestral_line_value.t() << std::endl;
      }
    }
    delete my_move_al;
    return Rcpp::List::create(Rcpp::Named("remaped_ai") = ancestral_lines_indeces,
                              Rcpp::Named("ancestorLines") = ancestral_lines_values);
  }
  catch(smc::exception  e) {
    Rcpp::Rcout << e;
  }
  return R_NilValue; // to provide a return
}
