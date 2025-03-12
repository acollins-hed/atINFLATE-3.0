#ifndef GET_MUT_TRANS_HPP
#define GET_MUT_TRANS_HPP

#include<iostream>
#include<algorithm>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "common_variables.hpp"
#include "Mutant.hpp"
#include "fitness_functions.hpp"
#include "Transition_Probability.hpp"
#include "Population.hpp"

void get_trans_prob_rd(std::vector<Mutant*> mut_vec, struct Population * population, struct Common_Variables * common_variables, long double * sum);

void get_trans_prob_ri(std::vector<Mutant*> mut_vec, struct Population * population, struct Common_Variables * common_variables, long double * sum);

#endif
