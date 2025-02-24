#ifndef EVOLVER_HPP
#define EVOLVER_HPP

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Documents.hpp"
#include "Mutant.hpp"
#include "fitness_functions.hpp"
//#include "Stats.hpp"
#include "inverted_hsum.hpp"
#include "Transition_Probability.hpp"
#include "Population.hpp"

struct Evolver{
  struct Population population;
  std::vector<Mutant> mutant_vector;
  struct Documents documents;

  Evolver(struct Population pop);

  int Hoarse_Partition(int low, int high);

  void Quick_Sort(int low, int high);
  
  void Get_Mutants_rdm(struct Common_Variables * common_variables);

  void Get_Mutants_rim(struct Common_Variables * common_variables);

  void Get_Mutants_rdu(struct Common_Variables * common_variables);

  void Get_Mutants_riu(struct Common_Variables * common_variables);
  
  void Get_Next_TransMach(int trajectory, struct Common_Variables * common_variables);

  void Record_Initial_State(struct Common_Variables * common_variables);
  
  void Record_Data(int trajectory, int fixation, struct Common_Variables * common_variables);

  void Record_Final_State(struct Common_Variables * common_variables);

  void Set_New_Population_rdm(struct Common_Variables * common_variables);

  void Set_New_Population_rim(struct Common_Variables * common_variables);
  
  void Fix_rdm(struct Common_Variables * common_variables);

  void Fix_rim(struct Common_Variables * common_variables);

  void Fix_rdu(struct Common_Variables * common_variables);
 
  void Fix_riu(struct Common_Variables * common_variables);
  
  void Run_Simulation(struct Common_Variables * common_variables);
 
};

#endif
