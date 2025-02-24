#include<iostream>
#include<fstream>
#include<string>
#include<cfloat>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<unsupported/Eigen/KroneckerProduct>
#include<cmath>
#include<time.h>
#include<random>
#include<map>
#include<utility>
#include "Documents.hpp"
#include "convertToDouble.hpp"
#include "BadConversion.hpp"
#include "common_variables.hpp"
//#include "Stats.hpp"
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
//#include "atinflate_parcer.hpp"
#include "Read_Input_File.hpp"
#include "initialize_variables.hpp"
//compile with g++ -std=c++11 -O2 -I /PATH/TO/Eigen/ atinflate_3.cpp -o atinflate

//using namespace std;

int main(int argc, char* argv[]){
  clock_t start_time,end_time;
  start_time = clock();
  
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

  std::cout<<"The code is \n"<<genotype.code<<std::endl;

  struct Population population(genotype,&common_variables);

  if(common_variables.rate)
    population.fitness = Fitness_rate_dep(&population.codon_frequency,&population.genotype.code,&population.genotype.kd,&common_variables);
  else
    population.fitness = Fitness_rate_indep(&population.codon_frequency,&population.genotype.code,&common_variables);

  std::cout<<"codon frequencies are\n"<<population.codon_frequency<<std::endl;
  std::cout<<"population fitness is "<<population.fitness<<std::endl;

  struct Evolver evolver(population);//,&common_variables);
  //evolver.mutant_array = (struct Mutant *) malloc(common_variables.N_total_mutants*sizeof(struct Mutant));
  evolver.mutant_vector.resize(common_variables.N_total_mutants);

  std::cout<<"The total number of mutants is "<<common_variables.N_total_mutants<<std::endl;

  evolver.Run_Simulation(&common_variables);
  
  //free(evolver.mutant_array);
  end_time = clock();
  std::cout<<(double)(end_time-start_time)/CLOCKS_PER_SEC<<" seconds\n";
  std::cout<<(double)(end_time-start_time)/CLOCKS_PER_SEC/60<<" minutes\n";
  std::cout<<(double)(end_time-start_time)/CLOCKS_PER_SEC/3600<<" hours\n";
  return 0;
}
