#ifndef FITFUNCS_HPP
#define FITFUNCS_HPP

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "common_variables.hpp"
//#include "Stats.hpp"
#include "inverted_hsum.hpp"
#include "Transition_Probability.hpp"

double A_fit(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, int stype, struct Common_Variables * common_variables);

double Fitness_rate_indep(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, struct Common_Variables * common_variables);

double Fitness_rate_dep(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, Eigen::MatrixXd* kd, struct Common_Variables * common_variables);

#endif
