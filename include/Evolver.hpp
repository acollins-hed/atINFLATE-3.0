#ifndef EVOLVER_HPP
#define EVOLVER_HPP

#include<iostream>
#include<thread>
#include<algorithm>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Documents.hpp"
#include "Mutant.hpp"
#include "fitness_functions.hpp"
#include "inverted_hsum.hpp"
#include "Transition_Probability.hpp"
#include "Population.hpp"
#include "get_mutant_trans_probs.hpp"

struct Evolver{
  struct Population population;
  std::vector<Mutant> mutant_vector;
  struct Documents documents;

  Evolver(struct Population pop);

  int Hoarse_Partition(int low, int high);

  void Quick_Sort(int low, int high);

  void initialize_mutants(struct Common_Variables * common_variables);
  
  void Get_Mutants_rd(struct Common_Variables * common_variables);

  void Get_Mutants_ri(struct Common_Variables * common_variables);
  
  void Get_Next_TransMach(int trajectory, struct Common_Variables * common_variables);

  void Record_Initial_State(struct Common_Variables * common_variables);
  
  void Record_Data(int trajectory, int fixation, struct Common_Variables * common_variables);

  void Record_Final_State(struct Common_Variables * common_variables);

  void Set_New_Population_rd(struct Common_Variables * common_variables);

  void Set_New_Population_ri(struct Common_Variables * common_variables);
  
  void Fix_rd(struct Common_Variables * common_variables);

  void Fix_ri(struct Common_Variables * common_variables);
  
  void Run_Simulation(struct Common_Variables * common_variables);
 
};

#endif
