#include<iostream>
#include<fstream>
#include<string>
#include<cfloat>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<unsupported/Eigen/KroneckerProduct>
#include<cmath>
#include<chrono>
#include<random>
#include<map>
#include<utility>
#include "Documents.hpp"
#include "convertToDouble.hpp"
#include "BadConversion.hpp"
#include "common_variables.hpp"
#include "inverted_hsum.hpp"
#include "Transition_Probability.hpp"
#include "fitness_functions.hpp"
#include "tRNAs.hpp"
#include "aaRSs.hpp"
#include "Genotype.hpp"
#include "Get_Mutant_Vector_Size.hpp"
#include "generate_initial_trans_machinery.hpp"
#include "Mutation.hpp"
#include "Mutant.hpp"
#include "Population.hpp"
#include "Evolver.hpp"
#include "Read_Input_File.hpp"
#include "initialize_variables.hpp"

int main(int argc, char* argv[]){

  auto start = std::chrono::high_resolution_clock::now();

  struct Common_Variables common_variables;

  Eigen::initParallel();
  if(!initialize_variables(argc, argv,&common_variables))    
    return 0;
  common_variables.N_total_mutants = Get_Mutant_Vector_Size(&common_variables);
  if(common_variables.mask)
    common_variables.N_single_mutants = 2*(common_variables.N_tRNA+common_variables.N_aaRS)*common_variables.N_int_interface;
  else
    common_variables.N_single_mutants = (common_variables.N_tRNA+common_variables.N_aaRS)*common_variables.N_int_interface;

  if(!common_variables.bl_input_filename)
    gen_trans_machry(&common_variables);

  struct tRNAs trna(0,&common_variables);

  struct aaRSs aars(0,&common_variables);

  struct Genotype genotype(trna,aars,&common_variables);

  struct Population population(genotype,&common_variables);

  if(common_variables.rate)
    population.fitness = Fitness_rate_dep(&population.codon_frequency,&population.genotype.code,&population.genotype.kd,&common_variables);
  else
    population.fitness = Fitness_rate_indep(&population.codon_frequency,&population.genotype.code,&common_variables);

  struct Evolver evolver(population);

  evolver.mutant_vector.resize(common_variables.N_total_mutants);

  evolver.initialize_mutants(&common_variables);

  evolver.Run_Simulation(&common_variables);

  auto end = std::chrono::high_resolution_clock::now();
  double Time_taken = (end - start).count() * ((double) std::chrono::high_resolution_clock::period::num/std::chrono::high_resolution_clock::period::den);

  std::cout << "Run time: "<<Time_taken<<" seconds"<<std::endl;
  std::cout << "Run time: "<<Time_taken/60<<" minutes"<<std::endl;
  std::cout << "Run time: "<<Time_taken/3600<<" hours"<<std::endl;

  return 0;
}
