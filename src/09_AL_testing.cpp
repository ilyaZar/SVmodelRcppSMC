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
//' @param starting_vals arma::vec giving the three starting values for phi_x,
//'   sigma_x and beta_y
//' @param resample_freq int giving the number of progress outputs i.e.
//'   if set to 10, then progress output occurs for every additional 10% of 
//'   completion
//' @return Rcpp::List containing the results: parameter samples (sigma_x, 
//'   beta_y) and log-prior and log-likelihoood estimates
//'
//' @export
// [[Rcpp::export]]
Rcpp::List sv_model_al_tracking_impl(arma::vec measurements,
                                     arma::vec starting_vals,
                                     unsigned long lNumber,
                                     const double resample_freq = 0.5) {
  // Set variables related to progress monitoring of the overall procedure:
  long lIterates = y_pmmh_simul.n_rows;
  // int progress_intervall_num = round(lIterates/num_progress_outputs);
  // int check_prog        = 0;
  // double progress_ratio = 0.0;
  // Some other housekeeping
  y_pmmh_simul = measurements;
  // Set starting values:
  theta_prop.phi_x   = starting_vals(0);
  theta_prop.sigma_x = starting_vals(1);
  theta_prop.beta_y  = starting_vals(2);
  
  
  try {
    //Initialize and run the sampler
    my_move_al = new ALtracking_move;
    smc::sampler<double,smc::nullParams> Sampler(lNumber,
                                                 HistoryType::AL,
                                                 my_move_al);
    Sampler.SetResampleParams(ResampleType::MULTINOMIAL, resample_freq);
    Sampler.Initialise();
    // Sampler.IterateUntil(lIterates - 1);
    for (unsigned int i = 0; i < lIterates; ++i) {
      std::cout << i + 1 << " of " << lIterates << std::endl;
      smc::population<double> current_population = Sampler.GetHistoryPopulation(i);
      std::vector<double> current_vals = current_population.GetValue();
      arma::Col<unsigned int> current_population_al = Sampler.GetuRSIndices();
      for(int j = 0; j < static_cast<int>(current_vals.size()); j++)
      {
      std::cout << "Particle no: " << j + 1 << ". Ancestor index " 
        << current_population_al[j] <<  " with value " << current_vals[j] 
        << std::endl;
      }
      if(i != lIterates - 1) Sampler.Iterate();
    }
    arma::Mat<unsigned int> ancestral_lines_indeces(lNumber, lIterates);
    arma::Mat<double> ancestral_lines_values(lNumber, lIterates);
    for(long i = 0; i < lNumber; ++i) {
      arma::Col<unsigned int> test_ancestral_line_index = Sampler.GetALineInd(i);
      ancestral_lines_indeces.row(i) = test_ancestral_line_index.t();
      arma::Col<double> test_ancestral_line_value = Sampler.GetALineSpace(i);
      ancestral_lines_values.row(i) = test_ancestral_line_value.t();
      std::cout << "Particle no: " << i + 1 << " ancestor line: re-mapped index and value." << std::endl;
      std::cout << test_ancestral_line_index.t() << std::endl;
      std::cout << test_ancestral_line_value.t() << std::endl;
    }
    delete my_move_al;
    return Rcpp::List::create(Rcpp::Named("remaped_ai") = ancestral_lines_indeces,
                              Rcpp::Named("ancestor_lines") = ancestral_lines_values);
  }
  catch(smc::exception  e) {
    Rcpp::Rcout << e;
  }
  return R_NilValue; // to provide a return
}