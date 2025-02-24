#ifndef MUTANT_HPP
#define MUTANT_HPP

//#include<iostream>
//#include<Eigen/Dense>
//#include<Eigen/Core>
#include "Mutation.hpp"
//#include "common_variables.hpp"
//#include "generate_initial_trans_machinery.hpp"

struct Mutant{
  struct Mutation mutation_1;
  struct Mutation mutation_2;
  struct Mutant *next;
  long double trans_prob;
  Mutant();

  Mutant(char tm, int which_one, int position, int mask);

  Mutant(char tm1, int wo1, int pos1, int mask1, char tm2, int wo2, int pos2, int mask2);
};

#endif
