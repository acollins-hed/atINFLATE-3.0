#include<iostream>
#include<thread>
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

Evolver::Evolver(struct Population pop)
  :population(pop){
}

//int Hoarse_Partition(struct Mutant * mutant_array, int low){
int Evolver::Hoarse_Partition(int low, int high){
  int i = low-1, j = high+1;
  int mid = low + (high - low)/2;
  struct Mutant pivot = mutant_vector[mid];
    
  while(true){
    do{
      i++;
    }while(mutant_vector[i].trans_prob > pivot.trans_prob);
    
    do{
      j--;
    }while(mutant_vector[j].trans_prob < pivot.trans_prob);
    
    if(i>=j)
      return j;
    
    std::swap(mutant_vector[i],mutant_vector[j]);
  }
  
}
  
//void Quick_Sort(struct Mutant * mutant_array, int low, int high){
void Evolver::Quick_Sort(int low, int high){
  if(low < high){
    int part_index = Hoarse_Partition(low,high);
    
    Quick_Sort(low,part_index);
    Quick_Sort(part_index+1,high);
  }
}

void Evolver::initialize_mutants(struct Common_Variables * common_variables){
  mutant_vector.resize(common_variables->N_total_mutants);
  //int N_int_interface = common_variables->N_int_interface;
  int N_sm = common_variables->N_single_mutants;
  int index = 0;
  
  //Finding single mutation mutants
  for(int row=0; row<2; row++){
    for(int col=0; col<population.genotype.trnas.iis.cols(); col++){
      for(int i=0;i<population.genotype.trnas.N_int_interface;i++){
	mutant_vector[index].mutation_1.kind_trans_machinery = 't';
	mutant_vector[index].mutation_2.kind_trans_machinery = '0';
	mutant_vector[index].mutation_1.which_trans_machinery=col;
	mutant_vector[index].mutation_1.mutation_position=i;
	mutant_vector[index].mutation_1.mask=row;
	
	index++;
      }
    }
  }
  
  for(int row=0; row<2; row++){
    for(int col=0; col<population.genotype.aarss.iis.cols(); col++){
      for(int i=0;i<population.genotype.aarss.N_int_interface;i++){
	mutant_vector[index].mutation_1.kind_trans_machinery = 'a';
	mutant_vector[index].mutation_2.kind_trans_machinery = '0';
	mutant_vector[index].mutation_1.which_trans_machinery=col;
	mutant_vector[index].mutation_1.mutation_position=i;
	mutant_vector[index].mutation_1.mask=row;

	index++;
      }
    }
  }


  for(int i=0;i<N_sm;i++){
    for(int j=i+1;j<N_sm;j++){
      
      mutant_vector[index].mutation_1.kind_trans_machinery = mutant_vector[i].mutation_1.kind_trans_machinery;
      mutant_vector[index].mutation_1.which_trans_machinery = mutant_vector[i].mutation_1.which_trans_machinery;
      mutant_vector[index].mutation_1.mutation_position = mutant_vector[i].mutation_1.mutation_position;
      mutant_vector[index].mutation_1.mask = mutant_vector[i].mutation_1.mask;
      
      mutant_vector[index].mutation_2.kind_trans_machinery = mutant_vector[j].mutation_1.kind_trans_machinery;
      mutant_vector[index].mutation_2.which_trans_machinery = mutant_vector[j].mutation_1.which_trans_machinery;
      mutant_vector[index].mutation_2.mutation_position = mutant_vector[j].mutation_1.mutation_position;
      mutant_vector[index].mutation_2.mask = mutant_vector[i].mutation_1.mask;
      
      index++;
    }
  }
  
}

void Evolver::Get_Mutants_rdm(struct Common_Variables * common_variables){
  long double sum=0;//for normalizing the sum of transition probabilities
  struct Mutant mutant;
  int index = 0;
  common_variables->prob_of_1_mutation=0;
  common_variables->prob_of_2_mutations=0;
  
  //Finding single mutation mutants
  
  for(int row=0; row<2; row++){
    for(int col=0; col<population.genotype.trnas.iis.cols(); col++){
      for(int i=0;i<population.genotype.trnas.N_int_interface;i++){
	struct Population mutant_pop = population;
	mutant_pop.genotype.trnas.iis(row,col) ^= 1<<i;
	mutant_pop.genotype.Get_Code();
	mutant_pop.fitness = Fitness_rate_dep(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd,common_variables);
	mutant_vector[index].mutation_1.kind_trans_machinery = 't';
	mutant_vector[index].mutation_2.kind_trans_machinery = '0';
	mutant_vector[index].mutation_1.which_trans_machinery=col;
	mutant_vector[index].mutation_1.mutation_position=i;
	mutant_vector[index].mutation_1.mask=row;
	long double trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,1,common_variables);
	mutant_vector[index].trans_prob = trans_prob;
	sum += trans_prob;
	
	index++;
      }
    }
  }
  /*
  std::ofstream mut_file;
  mut_file.open("mut1.dat",std::ios_base::app);
  for(int i=0;i<common_variables->N_single_mutants/2;i++){
    mut_file<<i<<" "<<mutant_vector[i].trans_prob<<std::endl;
  }
  mut_file.close();
  */
  for(int row=0; row<2; row++){
    for(int col=0; col<population.genotype.aarss.iis.cols(); col++){
      for(int i=0;i<population.genotype.aarss.N_int_interface;i++){
	struct Population mutant_pop = population;
	mutant_pop.genotype.aarss.iis(row,col) ^= 1<<i;
	mutant_pop.genotype.Get_Code();
	mutant_pop.fitness = Fitness_rate_dep(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd,common_variables);
	mutant_vector[index].mutation_1.kind_trans_machinery = 'a';
	mutant_vector[index].mutation_2.kind_trans_machinery = '0';
	mutant_vector[index].mutation_1.which_trans_machinery=col;
	mutant_vector[index].mutation_1.mutation_position=i;
	mutant_vector[index].mutation_1.mask=row;
	long double trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,1,common_variables);
	mutant_vector[index].trans_prob = trans_prob;
	sum += trans_prob;
	
	index++;
      }
    }
  }
  
  common_variables->prob_of_1_mutation = sum;
  
  //Finding double mutation mutants
  if(common_variables->bl_double_mutants){
    for(int i=0;i<common_variables->N_single_mutants;i++){
      
      for(int j=i+1;j<common_variables->N_single_mutants;j++){
	struct Population mutant_pop = population;
	long double trans_prob=0;
	if(mutant_vector[i].mutation_1.kind_trans_machinery == 't')
	  mutant_pop.genotype.trnas.iis(mutant_vector[i].mutation_1.mask,mutant_vector[i].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[i].mutation_1.mutation_position);
	else
	  if(mutant_vector[i].mutation_1.kind_trans_machinery == 'a')
	    mutant_pop.genotype.aarss.iis(mutant_vector[i].mutation_1.mask,mutant_vector[i].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[i].mutation_1.mutation_position);
	else
	  std::cout<<"Error, some mutation is neither tRNA nor aaRS.\n";
	
	
	mutant_vector[index].mutation_1.kind_trans_machinery = mutant_vector[i].mutation_1.kind_trans_machinery;
	mutant_vector[index].mutation_1.which_trans_machinery = mutant_vector[i].mutation_1.which_trans_machinery;
	mutant_vector[index].mutation_1.mutation_position = mutant_vector[i].mutation_1.mutation_position;
	mutant_vector[index].mutation_1.mask = mutant_vector[i].mutation_1.mask;
	
	if(mutant_vector[j].mutation_1.kind_trans_machinery == 't')
	  mutant_pop.genotype.trnas.iis(mutant_vector[j].mutation_1.mask,mutant_vector[j].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[j].mutation_1.mutation_position);
	else
	  if(mutant_vector[j].mutation_1.kind_trans_machinery == 'a')
	    mutant_pop.genotype.aarss.iis(mutant_vector[j].mutation_1.mask,mutant_vector[j].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[j].mutation_1.mutation_position);
	  else
	    std::cout<<"Error, some mutation is neither tRNA nor aaRS.\n";
	
	mutant_vector[index].mutation_2.kind_trans_machinery = mutant_vector[j].mutation_1.kind_trans_machinery;
	mutant_vector[index].mutation_2.which_trans_machinery = mutant_vector[j].mutation_1.which_trans_machinery;
	mutant_vector[index].mutation_2.mutation_position = mutant_vector[j].mutation_1.mutation_position;
	mutant_vector[index].mutation_2.mask = mutant_vector[j].mutation_1.mask;
	
	mutant_pop.genotype.Get_Code();
	mutant_pop.fitness = Fitness_rate_dep(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd,common_variables);
	trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,2,common_variables);
	mutant_vector[index].trans_prob = trans_prob; 
	
	index++;
	common_variables->prob_of_2_mutations += mutant_vector[index].trans_prob;
	sum += trans_prob;
      }
    }
  }
  /*
  std::cout<<"mutant_vector[last].trans_prob = "<<mutant_vector[mutant_vector.size()-1].trans_prob<<std::endl;
  std::ofstream mut_file;
  mut_file.open("mut_file.dat",std::ios_base::app);
  for(unsigned int i=0; i<mutant_vector.size();i++){
    mut_file<<mutant_vector[i].mutation_1.kind_trans_machinery<<" "<<mutant_vector[i].mutation_1.which_trans_machinery<<" "<<mutant_vector[i].mutation_1.mask<<" "<<mutant_vector[i].mutation_1.mutation_position<<" "<<mutant_vector[i].mutation_2.kind_trans_machinery<<" "<<mutant_vector[i].mutation_2.which_trans_machinery<<" "<<mutant_vector[i].mutation_2.mask<<" "<<mutant_vector[i].mutation_2.mutation_position<<" "<<mutant_vector[i].trans_prob<<std::endl;
  }
  mut_file.close();
  
  std::cout<<"sum is "<<sum<<std::endl;
  */
  for(int i=0;i<index;i++){
    mutant_vector[i].trans_prob /= sum;
  }
  
  common_variables->prob_of_1_mutation = common_variables->prob_of_1_mutation/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_2_mutations = common_variables->prob_of_2_mutations/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_1_mutation = common_variables->prob_of_1_mutation/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_2_mutations = common_variables->prob_of_2_mutations/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  
  Quick_Sort(0,common_variables->N_total_mutants - 1);
}

void Evolver::Get_Mutants_rim(struct Common_Variables * common_variables){
  long double sum=0;//for normalizing the sum of transition probabilities
  struct Mutant mutant;
  int index = 0;
  common_variables->prob_of_1_mutation=0;
  common_variables->prob_of_2_mutations=0;
  
  //Finding single mutation mutants
  for(int row=0; row<2; row++){
    for(int col=0; col<population.genotype.trnas.iis.cols(); col++){
      for(int i=0;i<population.genotype.trnas.N_int_interface;i++){
	struct Population mutant_pop = population;
	mutant_pop.genotype.trnas.iis(row,col) ^= 1<<i;
	mutant_pop.genotype.Get_Code();
	mutant_pop.fitness = Fitness_rate_indep(&population.codon_frequency,&mutant_pop.genotype.code,common_variables);
	mutant_vector[index].mutation_1.kind_trans_machinery = 't';
	mutant_vector[index].mutation_2.kind_trans_machinery = '0';
	mutant_vector[index].mutation_1.which_trans_machinery=col;
	mutant_vector[index].mutation_1.mutation_position=i;
	mutant_vector[index].mutation_1.mask=row;
	long double trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,1,common_variables);
	mutant_vector[index].trans_prob = trans_prob;
	sum += trans_prob;
	
	index++;
      }
    }
  }
  
  for(int row=0; row<2; row++){
    for(int col=0; col<population.genotype.aarss.iis.cols(); col++){
      for(int i=0;i<population.genotype.aarss.N_int_interface;i++){
	struct Population mutant_pop = population;
	mutant_pop.genotype.aarss.iis(row,col) ^= 1<<i;
	mutant_pop.genotype.Get_Code();
	mutant_pop.fitness = Fitness_rate_indep(&population.codon_frequency,&mutant_pop.genotype.code,common_variables);
	mutant_vector[index].mutation_1.kind_trans_machinery = 'a';
	mutant_vector[index].mutation_2.kind_trans_machinery = '0';
	mutant_vector[index].mutation_1.which_trans_machinery=col;
	mutant_vector[index].mutation_1.mutation_position=i;
	mutant_vector[index].mutation_1.mask=row;
	long double trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,1,common_variables);
	mutant_vector[index].trans_prob = trans_prob;
	sum += trans_prob;
	
	index++;
      }
    }
  }
  
  common_variables->prob_of_1_mutation = sum;
  
  //Finding double mutation mutants
  if(common_variables->bl_double_mutants){
    for(int i=0;i<common_variables->N_single_mutants;i++){
      
      for(int j=i+1;j<common_variables->N_single_mutants;j++){

	struct Population mutant_pop = population;
	long double trans_prob=0;
	if(mutant_vector[i].mutation_1.kind_trans_machinery == 't')
	  mutant_pop.genotype.trnas.iis(mutant_vector[i].mutation_1.mask,mutant_vector[i].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[i].mutation_1.mutation_position);
	else
	  if(mutant_vector[i].mutation_1.kind_trans_machinery == 'a')
	    mutant_pop.genotype.aarss.iis(mutant_vector[i].mutation_1.mask,mutant_vector[i].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[i].mutation_1.mutation_position);
	  else
	    std::cout<<"Error, some mutation is neither tRNA nor aaRS.\n";
	
	mutant_vector[index].mutation_1.kind_trans_machinery = mutant_vector[i].mutation_1.kind_trans_machinery;
	mutant_vector[index].mutation_1.which_trans_machinery = mutant_vector[i].mutation_1.which_trans_machinery;
	mutant_vector[index].mutation_1.mutation_position = mutant_vector[i].mutation_1.mutation_position;
	mutant_vector[index].mutation_1.mask = mutant_vector[i].mutation_1.mask;
	
	if(mutant_vector[j].mutation_1.kind_trans_machinery == 't')
	  mutant_pop.genotype.trnas.iis(mutant_vector[j].mutation_1.mask,mutant_vector[j].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[j].mutation_1.mutation_position);
	else
	  if(mutant_vector[j].mutation_1.kind_trans_machinery == 'a')
	    mutant_pop.genotype.aarss.iis(mutant_vector[j].mutation_1.mask,mutant_vector[j].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[j].mutation_1.mutation_position);
	  else
	    std::cout<<"Error, some mutation is neither tRNA nor aaRS.\n";
	
	mutant_vector[index].mutation_2.kind_trans_machinery = mutant_vector[j].mutation_1.kind_trans_machinery;
	mutant_vector[index].mutation_2.which_trans_machinery = mutant_vector[j].mutation_1.which_trans_machinery;
	mutant_vector[index].mutation_2.mutation_position = mutant_vector[j].mutation_1.mutation_position;
	mutant_vector[index].mutation_2.mask = mutant_vector[i].mutation_1.mask;
	
	mutant_pop.genotype.Get_Code();
	mutant_pop.fitness = Fitness_rate_indep(&population.codon_frequency,&mutant_pop.genotype.code,common_variables);
	trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,2,common_variables);
	mutant_vector[index].trans_prob = trans_prob; 
	
	index++;
	common_variables->prob_of_2_mutations += mutant_vector[index].trans_prob;
	sum += trans_prob;
      }
    }
  }
  
  for(int i=0;i<index;i++){
    mutant_vector[i].trans_prob /= sum;
  }
  
  common_variables->prob_of_1_mutation = common_variables->prob_of_1_mutation/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_2_mutations = common_variables->prob_of_2_mutations/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_1_mutation = common_variables->prob_of_1_mutation/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_2_mutations = common_variables->prob_of_2_mutations/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  
  Quick_Sort(0,common_variables->N_total_mutants - 1);
}

void Evolver::Get_Mutants_rdu(struct Common_Variables * common_variables){
  long double sum=0;//for normalizing the sum of transition probabilities
  struct Mutant mutant;
  int index = 0;
  common_variables->prob_of_1_mutation=0;
  common_variables->prob_of_2_mutations=0;
  
  //Finding single mutation mutants
  for(int col=0; col<population.genotype.trnas.iis.cols(); col++){
    for(int i=0;i<population.genotype.trnas.N_int_interface;i++){
      struct Population mutant_pop = population;
      mutant_pop.genotype.trnas.iis(0,col) ^= 1<<i;
      mutant_pop.genotype.Get_Code();
      mutant_pop.fitness = Fitness_rate_dep(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd,common_variables);
      mutant_vector[index].mutation_1.kind_trans_machinery = 't';
      mutant_vector[index].mutation_2.kind_trans_machinery = '0';
      mutant_vector[index].mutation_1.which_trans_machinery=col;
      mutant_vector[index].mutation_1.mutation_position=i;
      mutant_vector[index].mutation_1.mask=0;
      long double trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,1,common_variables);
      mutant_vector[index].trans_prob = trans_prob;
      sum += trans_prob;
      
      index++;
    }
  }
  
  
  for(int col=0; col<population.genotype.aarss.iis.cols(); col++){
    for(int i=0;i<population.genotype.aarss.N_int_interface;i++){
      struct Population mutant_pop = population;
      mutant_pop.genotype.aarss.iis(0,col) ^= 1<<i;
      mutant_pop.genotype.Get_Code();
      mutant_pop.fitness = Fitness_rate_dep(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd,common_variables);
      mutant_vector[index].mutation_1.kind_trans_machinery = 'a';
      mutant_vector[index].mutation_2.kind_trans_machinery = '0';
      mutant_vector[index].mutation_1.which_trans_machinery=col;
      mutant_vector[index].mutation_1.mutation_position=i;
      mutant_vector[index].mutation_1.mask=0;
      long double trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,1,common_variables);
      mutant_vector[index].trans_prob = trans_prob;
      sum += trans_prob;
      
      index++;
    }
  }
  
  
  common_variables->prob_of_1_mutation = sum;
  
  //Finding double mutation mutants
  if(common_variables->bl_double_mutants){
    for(int i=0;i<common_variables->N_single_mutants;i++){
      
      for(int j=i+1;j<common_variables->N_single_mutants;j++){

	struct Population mutant_pop = population;
	long double trans_prob=0;
	if(mutant_vector[i].mutation_1.kind_trans_machinery == 't')
	  mutant_pop.genotype.trnas.iis(mutant_vector[i].mutation_1.mask,mutant_vector[i].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[i].mutation_1.mutation_position);
	else
	  if(mutant_vector[i].mutation_1.kind_trans_machinery == 'a')
	    mutant_pop.genotype.aarss.iis(mutant_vector[i].mutation_1.mask,mutant_vector[i].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[i].mutation_1.mutation_position);
	  else
	    std::cout<<"Error, some mutation is neither tRNA nor aaRS.\n";
	
	mutant_vector[index].mutation_1.kind_trans_machinery = mutant_vector[i].mutation_1.kind_trans_machinery;
	mutant_vector[index].mutation_1.which_trans_machinery = mutant_vector[i].mutation_1.which_trans_machinery;
	mutant_vector[index].mutation_1.mutation_position = mutant_vector[i].mutation_1.mutation_position;
	mutant_vector[index].mutation_1.mask = mutant_vector[i].mutation_1.mask;
	
	if(mutant_vector[j].mutation_1.kind_trans_machinery == 't')
	  mutant_pop.genotype.trnas.iis(mutant_vector[j].mutation_1.mask,mutant_vector[j].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[j].mutation_1.mutation_position);
	else
	  if(mutant_vector[j].mutation_1.kind_trans_machinery == 'a')
	    mutant_pop.genotype.aarss.iis(mutant_vector[j].mutation_1.mask,mutant_vector[j].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[j].mutation_1.mutation_position);
	  else
	    std::cout<<"Error, some mutation is neither tRNA nor aaRS.\n";
	
	mutant_vector[index].mutation_2.kind_trans_machinery = mutant_vector[j].mutation_1.kind_trans_machinery;
	mutant_vector[index].mutation_2.which_trans_machinery = mutant_vector[j].mutation_1.which_trans_machinery;
	mutant_vector[index].mutation_2.mutation_position = mutant_vector[j].mutation_1.mutation_position;
	mutant_vector[index].mutation_2.mask = mutant_vector[i].mutation_1.mask;
	
	mutant_pop.genotype.Get_Code();
	mutant_pop.fitness = Fitness_rate_dep(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd,common_variables);
	trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,2,common_variables);
	mutant_vector[index].trans_prob = trans_prob; 
	
	index++;
	common_variables->prob_of_2_mutations += mutant_vector[index].trans_prob;
	sum += trans_prob;
      }
    }
  }
  
  for(int i=0;i<index;i++){
    mutant_vector[i].trans_prob /= sum;
  }
  
  common_variables->prob_of_1_mutation = common_variables->prob_of_1_mutation/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_2_mutations = common_variables->prob_of_2_mutations/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_1_mutation = common_variables->prob_of_1_mutation/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_2_mutations = common_variables->prob_of_2_mutations/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  
  Quick_Sort(0,common_variables->N_total_mutants - 1);
}

void Evolver::Get_Mutants_riu(struct Common_Variables * common_variables){
  long double sum=0;//for normalizing the sum of transition probabilities
  struct Mutant mutant;
  int index = 0;
  common_variables->prob_of_1_mutation=0;
  common_variables->prob_of_2_mutations=0;
  
  //Finding single mutation mutants
  for(int col=0; col<population.genotype.trnas.iis.cols(); col++){
    for(int i=0;i<population.genotype.trnas.N_int_interface;i++){
      struct Population mutant_pop = population;
      mutant_pop.genotype.trnas.iis(0,col) ^= 1<<i;
      mutant_pop.genotype.Get_Code();
      mutant_pop.fitness = Fitness_rate_indep(&population.codon_frequency,&mutant_pop.genotype.code,common_variables);
      mutant_vector[index].mutation_1.kind_trans_machinery = 't';
      mutant_vector[index].mutation_2.kind_trans_machinery = '0';
      mutant_vector[index].mutation_1.which_trans_machinery=col;
      mutant_vector[index].mutation_1.mutation_position=i;
      mutant_vector[index].mutation_1.mask=0;
      long double trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,1,common_variables);
      mutant_vector[index].trans_prob = trans_prob;
      sum += trans_prob;
      
      index++;
    }
  }
  
  
  for(int col=0; col<population.genotype.aarss.iis.cols(); col++){
    for(int i=0;i<population.genotype.aarss.N_int_interface;i++){
      struct Population mutant_pop = population;
      mutant_pop.genotype.aarss.iis(0,col) ^= 1<<i;
      mutant_pop.genotype.Get_Code();
      mutant_pop.fitness = Fitness_rate_indep(&population.codon_frequency,&mutant_pop.genotype.code,common_variables);
      mutant_vector[index].mutation_1.kind_trans_machinery = 'a';
      mutant_vector[index].mutation_2.kind_trans_machinery = '0';
      mutant_vector[index].mutation_1.which_trans_machinery=col;
      mutant_vector[index].mutation_1.mutation_position=i;
      mutant_vector[index].mutation_1.mask=0;
      long double trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,1,common_variables);
      mutant_vector[index].trans_prob = trans_prob;
      sum += trans_prob;
      
      index++;
    }
  }
  
  
  common_variables->prob_of_1_mutation = sum;
  
  //Finding double mutation mutants
  if(common_variables->bl_double_mutants){
    for(int i=0;i<common_variables->N_single_mutants;i++){
      
      for(int j=i+1;j<common_variables->N_single_mutants;j++){

	struct Population mutant_pop = population;
	long double trans_prob=0;
	if(mutant_vector[i].mutation_1.kind_trans_machinery == 't')
	  mutant_pop.genotype.trnas.iis(mutant_vector[i].mutation_1.mask,mutant_vector[i].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[i].mutation_1.mutation_position);
	else
	  if(mutant_vector[i].mutation_1.kind_trans_machinery == 'a')
	    mutant_pop.genotype.aarss.iis(mutant_vector[i].mutation_1.mask,mutant_vector[i].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[i].mutation_1.mutation_position);
	  else
	    std::cout<<"Error, some mutation is neither tRNA nor aaRS.\n";
	
	mutant_vector[index].mutation_1.kind_trans_machinery = mutant_vector[i].mutation_1.kind_trans_machinery;
	mutant_vector[index].mutation_1.which_trans_machinery = mutant_vector[i].mutation_1.which_trans_machinery;
	mutant_vector[index].mutation_1.mutation_position = mutant_vector[i].mutation_1.mutation_position;
	mutant_vector[index].mutation_1.mask = mutant_vector[i].mutation_1.mask;
	
	if(mutant_vector[j].mutation_1.kind_trans_machinery == 't')
	  mutant_pop.genotype.trnas.iis(mutant_vector[j].mutation_1.mask,mutant_vector[j].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[j].mutation_1.mutation_position);
	else
	  if(mutant_vector[j].mutation_1.kind_trans_machinery == 'a')
	    mutant_pop.genotype.aarss.iis(mutant_vector[j].mutation_1.mask,mutant_vector[j].mutation_1.which_trans_machinery) ^= 1<<(mutant_vector[j].mutation_1.mutation_position);
	  else
	    std::cout<<"Error, some mutation is neither tRNA nor aaRS.\n";
	
	mutant_vector[index].mutation_2.kind_trans_machinery = mutant_vector[j].mutation_1.kind_trans_machinery;
	mutant_vector[index].mutation_2.which_trans_machinery = mutant_vector[j].mutation_1.which_trans_machinery;
	mutant_vector[index].mutation_2.mutation_position = mutant_vector[j].mutation_1.mutation_position;
	mutant_vector[index].mutation_2.mask = mutant_vector[i].mutation_1.mask;
	
	mutant_pop.genotype.Get_Code();
	mutant_pop.fitness = Fitness_rate_indep(&population.codon_frequency,&mutant_pop.genotype.code,common_variables);
	trans_prob = Transition_Probability(population.fitness,mutant_pop.fitness,2,common_variables);
	mutant_vector[index].trans_prob = trans_prob; 
	
	index++;
	common_variables->prob_of_2_mutations += mutant_vector[index].trans_prob;
	sum += trans_prob;
      }
    }
  }
  
  for(int i=0;i<index;i++){
    mutant_vector[i].trans_prob /= sum;
  }
  
  common_variables->prob_of_1_mutation = common_variables->prob_of_1_mutation/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_2_mutations = common_variables->prob_of_2_mutations/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_1_mutation = common_variables->prob_of_1_mutation/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_2_mutations = common_variables->prob_of_2_mutations/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  
  Quick_Sort(0,common_variables->N_total_mutants - 1);
}

  
void Evolver::Get_Next_TransMach(int trajectory, struct Common_Variables * common_variables){
  population.genotype.trnas.iis.row(0) = common_variables->tRNA_State_bits.row(trajectory);
  population.genotype.trnas.iis.row(1) = common_variables->tRNA_Mask_bits.row(trajectory);
  population.genotype.aarss.iis.row(0) = common_variables->aaRS_State_bits.row(trajectory);
  population.genotype.aarss.iis.row(1) = common_variables->aaRS_Mask_bits.row(trajectory);
}

void Evolver::Record_Initial_State(struct Common_Variables * common_variables){
  //records variables to a checkpoint file
  
  documents.ocheckpoint_file<<common_variables->output_filename<<std::endl;
  documents.ocheckpoint_file<<common_variables->N_int_interface<<std::endl;
  documents.ocheckpoint_file<<common_variables->N_tRNA<<std::endl;
  documents.ocheckpoint_file<<common_variables->N_aaRS<<std::endl;
  documents.ocheckpoint_file<<common_variables->N_site_type<<std::endl;
  documents.ocheckpoint_file<<common_variables->phi<<std::endl;
  documents.ocheckpoint_file<<common_variables->transition_bias<<std::endl;
  documents.ocheckpoint_file<<common_variables->N_population<<std::endl;
  for(int i=0;i<common_variables->N_site_type;i++)
    documents.ocheckpoint_file<<common_variables->site_type_freqs(i)<<" ";
  documents.ocheckpoint_file<<std::endl;
  documents.ocheckpoint_file<<common_variables->rate<<std::endl;
  documents.ocheckpoint_file<<common_variables->mask<<std::endl;
  documents.ocheckpoint_file<<common_variables->proofreading<<std::endl;
  documents.ocheckpoint_file<<common_variables->rate_constant<<std::endl;
  documents.ocheckpoint_file<<common_variables->mu_id_feat<<std::endl;
  documents.ocheckpoint_file<<common_variables->mu_per_codon<<std::endl;
  if(common_variables->proofreading){
    documents.ocheckpoint_file<<sqrt(common_variables->kmax)<<std::endl;
    documents.ocheckpoint_file<<sqrt(common_variables->kmin)<<std::endl;
  } else{
    documents.ocheckpoint_file<<common_variables->kmax<<std::endl;
    documents.ocheckpoint_file<<common_variables->kmin<<std::endl;
  }
  for(int i = 0; i<common_variables->N_site_type;i++)
    documents.ocheckpoint_file<<common_variables->site_types(i)<<" ";
  documents.ocheckpoint_file<<std::endl;
  documents.ocheckpoint_file<<common_variables->codon_ring_space<<std::endl;
  documents.ocheckpoint_file<<common_variables->N_base<<std::endl;
  documents.ocheckpoint_file<<common_variables->binom_p<<std::endl;
  documents.ocheckpoint_file<<common_variables->N_trajectory<<std::endl;
  documents.ocheckpoint_file<<common_variables->halting_fixation<<std::endl;
  documents.ocheckpoint_file<<common_variables->halting_fitness<<std::endl;
}
  
void Evolver::Record_Data(int trajectory, int fixation, struct Common_Variables * common_variables){
  double frac_on=0;
  //trajectory file
  for(int i=0; i<common_variables->N_tRNA;i++)
    frac_on += __builtin_popcount(population.genotype.trnas.iis(1,i));
  for(int i=0;i<common_variables->N_aaRS;i++)
    frac_on += __builtin_popcount(population.genotype.aarss.iis(1,i));
  frac_on = frac_on/(common_variables->N_int_interface*(common_variables->N_tRNA+common_variables->N_aaRS));
  documents.traj_file<<trajectory<<" "<<fixation<<" "<<population.fitness<<" "<<frac_on<<" "<<common_variables->mutation_type<<" "<<common_variables->prob_of_1_mutation<<" "<<common_variables->prob_of_2_mutations<<std::endl;
  
  //code file
  for(int trna_no=0;trna_no<common_variables->N_tRNA;trna_no++){
    for(int aars_no=0;aars_no<common_variables->N_aaRS;aars_no++){
      documents.code_file<<trajectory<<" "<<fixation<<" "<<trna_no<<" "<<common_variables->amino_acids(aars_no)<<" "<<population.genotype.code(trna_no,aars_no)<<" "<<population.genotype.matches(trna_no,aars_no)<<" "<<population.genotype.kd(trna_no,aars_no)<<std::endl;
    }
  }
  
  //interface file
  for(int i = 0; i<common_variables->N_tRNA;i++){
    documents.int_file<<trajectory<<" "<<fixation<<" "<<"tRNA "<<i<<" State "<<population.genotype.trnas.get_trna_int(0,i)<<" "<<population.genotype.trnas.iis(0,i)<<std::endl;
  }
  for(int i = 0; i<common_variables->N_tRNA;i++){
    documents.int_file<<trajectory<<" "<<fixation<<" "<<"tRNA "<<i<<" Mask "<<population.genotype.trnas.get_trna_int(1,i)<<" "<<population.genotype.trnas.iis(1,i)<<std::endl;
  }
  for(int i = 0; i<common_variables->N_aaRS;i++){
    documents.int_file<<trajectory<<" "<<fixation<<" "<<"aaRS "<<common_variables->amino_acids(i)<<" State "<<population.genotype.aarss.get_aars_int(0,i)<<" "<<population.genotype.trnas.iis(0,i)<<std::endl;
  }
  for(int i = 0; i<common_variables->N_aaRS;i++){
    documents.int_file<<trajectory<<" "<<fixation<<" "<<"aaRS "<<common_variables->amino_acids(i)<<" Mask "<<population.genotype.aarss.get_aars_int(1,i)<<" "<<population.genotype.aarss.iis(1,i)<<std::endl;
  }
  
  //codon file
  for(int stype_no=0;stype_no<common_variables->N_site_type;stype_no++){
    for(int codon_no=0;codon_no<common_variables->N_tRNA;codon_no++){
      documents.codon_file<<trajectory<<" "<<fixation<<" "<<common_variables->site_types(stype_no)<<" "<<codon_no<<" "<<population.codon_frequency(stype_no,codon_no)<<std::endl;
    }
  }
  
  //probability file
  Eigen::MatrixXd Mod(common_variables->N_site_type,common_variables->N_aaRS);
  Mod = population.codon_frequency*population.genotype.code;
  for(int s=0;s<common_variables->N_site_type;s++){
    for(int alpha=0;alpha<common_variables->N_aaRS;alpha++){
      documents.prob_file<<trajectory<<" "<<fixation<<" "<<common_variables->site_types(s)<<" "<<common_variables->amino_acids(alpha)<<" "<<Mod(s,alpha)<<"\n";
    }
  }
    
}

void Evolver::Record_Final_State(struct Common_Variables * common_variables){
  for(int i=0; i<common_variables->N_tRNA;i++)
    documents.ocheckpoint_file<<population.genotype.trnas.iis(0,i)<<" ";
  documents.ocheckpoint_file<<std::endl;
  for(int i=0; i<common_variables->N_tRNA;i++)
    documents.ocheckpoint_file<<population.genotype.trnas.iis(1,i)<<" ";
  documents.ocheckpoint_file<<std::endl;
  for(int i=0; i<common_variables->N_aaRS;i++)
    documents.ocheckpoint_file<<population.genotype.aarss.iis(0,i)<<" ";
  documents.ocheckpoint_file<<std::endl;
  for(int i=0; i<common_variables->N_aaRS;i++)
    documents.ocheckpoint_file<<population.genotype.aarss.iis(1,i)<<" ";
  documents.ocheckpoint_file<<std::endl;
}

void Evolver::Set_New_Population_rdm(struct Common_Variables * common_variables){
  population.genotype.Get_Code();
  population.Get_Codon_Freq();
  population.fitness = Fitness_rate_dep(&population.codon_frequency,&population.genotype.code,&population.genotype.kd,common_variables);
}

void Evolver::Set_New_Population_rim(struct Common_Variables * common_variables){
  population.genotype.Get_Code();
  population.Get_Codon_Freq();
  population.fitness = Fitness_rate_indep(&population.codon_frequency,&population.genotype.code,common_variables);
}

void Evolver::Fix_rdm(struct Common_Variables * common_variables){
  long double rand_n = common_variables->gillespie_rand(common_variables->mersene_twister), sum = 0;
  long unsigned int mutant_num=0;
  
  common_variables->mutation_type = "single";
  
  Get_Mutants_rdm(common_variables);
  
  for(long unsigned int nth_mutant=0; nth_mutant<mutant_vector.size();nth_mutant++){
    sum += mutant_vector[nth_mutant].trans_prob;
    if(sum > rand_n){
      mutant_num=nth_mutant;
      break;
    }
  }
  
  if(mutant_vector[mutant_num].mutation_1.kind_trans_machinery == 't'){
    population.genotype.trnas.iis(mutant_vector[mutant_num].mutation_1.mask,mutant_vector[mutant_num].mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_1.mutation_position;
  }
  else{
    population.genotype.aarss.iis(mutant_vector[mutant_num].mutation_1.mask,mutant_vector[mutant_num].mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_1.mutation_position;
  }
  
  if(mutant_vector[mutant_num].mutation_2.kind_trans_machinery != '0'){
    if(mutant_vector[mutant_num].mutation_2.kind_trans_machinery == 't'){
      population.genotype.trnas.iis(mutant_vector[mutant_num].mutation_2.mask,mutant_vector[mutant_num].mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_2.mutation_position;
    }
    else{
      population.genotype.aarss.iis(mutant_vector[mutant_num].mutation_2.mask,mutant_vector[mutant_num].mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_2.mutation_position;
    }
    common_variables->mutation_type = "double";
  }
  
  Set_New_Population_rdm(common_variables);
}  

void Evolver::Fix_rim(struct Common_Variables * common_variables){
  long double rand_n = common_variables->gillespie_rand(common_variables->mersene_twister), sum = 0;
  long unsigned int mutant_num=0;
  
  common_variables->mutation_type = "single";
  
  Get_Mutants_rim(common_variables);
  
  for(long unsigned int nth_mutant=0; nth_mutant<mutant_vector.size();nth_mutant++){
    sum += mutant_vector[nth_mutant].trans_prob;
    if(sum > rand_n){
      mutant_num=nth_mutant;
      break;
    }
  }
  
  if(mutant_vector[mutant_num].mutation_1.kind_trans_machinery == 't'){
    population.genotype.trnas.iis(mutant_vector[mutant_num].mutation_1.mask,mutant_vector[mutant_num].mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_1.mutation_position;
  }
  else{
    population.genotype.aarss.iis(mutant_vector[mutant_num].mutation_1.mask,mutant_vector[mutant_num].mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_1.mutation_position;
  }
  
  if(mutant_vector[mutant_num].mutation_2.kind_trans_machinery != '0'){
    if(mutant_vector[mutant_num].mutation_2.kind_trans_machinery == 't'){
      population.genotype.trnas.iis(mutant_vector[mutant_num].mutation_2.mask,mutant_vector[mutant_num].mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_2.mutation_position;
    }
    else{
      population.genotype.aarss.iis(mutant_vector[mutant_num].mutation_2.mask,mutant_vector[mutant_num].mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_2.mutation_position;
    }
    common_variables->mutation_type = "double";
  }
  
  Set_New_Population_rim(common_variables);
}  

void Evolver::Fix_rdu(struct Common_Variables * common_variables){
  long double rand_n = common_variables->gillespie_rand(common_variables->mersene_twister), sum = 0;
  long unsigned int mutant_num=0;
  
  common_variables->mutation_type = "single";
  
  Get_Mutants_rdu(common_variables);
  
  for(long unsigned int nth_mutant=0; nth_mutant<mutant_vector.size();nth_mutant++){
    sum += mutant_vector[nth_mutant].trans_prob;
    if(sum > rand_n){
      mutant_num=nth_mutant;
      break;
    }
  }
  
  if(mutant_vector[mutant_num].mutation_1.kind_trans_machinery == 't'){
    population.genotype.trnas.iis(mutant_vector[mutant_num].mutation_1.mask,mutant_vector[mutant_num].mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_1.mutation_position;
  }
  else{
    population.genotype.aarss.iis(mutant_vector[mutant_num].mutation_1.mask,mutant_vector[mutant_num].mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_1.mutation_position;
  }
  
  if(mutant_vector[mutant_num].mutation_2.kind_trans_machinery != '0'){
    if(mutant_vector[mutant_num].mutation_2.kind_trans_machinery == 't'){
      population.genotype.trnas.iis(mutant_vector[mutant_num].mutation_2.mask,mutant_vector[mutant_num].mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_2.mutation_position;
    }
    else{
      population.genotype.aarss.iis(mutant_vector[mutant_num].mutation_2.mask,mutant_vector[mutant_num].mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_2.mutation_position;
    }
    common_variables->mutation_type = "double";
  }
  
  Set_New_Population_rdm(common_variables);
}  

void Evolver::Fix_riu(struct Common_Variables * common_variables){
  long double rand_n = common_variables->gillespie_rand(common_variables->mersene_twister), sum = 0;
  long unsigned int mutant_num=0;
  
  common_variables->mutation_type = "single";
  
  Get_Mutants_riu(common_variables);
  
  for(long unsigned int nth_mutant=0; nth_mutant<mutant_vector.size();nth_mutant++){
    sum += mutant_vector[nth_mutant].trans_prob;
    if(sum > rand_n){
      mutant_num=nth_mutant;
      break;
    }
  }
  
  if(mutant_vector[mutant_num].mutation_1.kind_trans_machinery == 't'){
    population.genotype.trnas.iis(mutant_vector[mutant_num].mutation_1.mask,mutant_vector[mutant_num].mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_1.mutation_position;
  }
  else{
    population.genotype.aarss.iis(mutant_vector[mutant_num].mutation_1.mask,mutant_vector[mutant_num].mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_1.mutation_position;
  }
  
  if(mutant_vector[mutant_num].mutation_2.kind_trans_machinery != '0'){
    if(mutant_vector[mutant_num].mutation_2.kind_trans_machinery == 't'){
      population.genotype.trnas.iis(mutant_vector[mutant_num].mutation_2.mask,mutant_vector[mutant_num].mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_2.mutation_position;
    }
    else{
      population.genotype.aarss.iis(mutant_vector[mutant_num].mutation_2.mask,mutant_vector[mutant_num].mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num].mutation_2.mutation_position;
    }
    common_variables->mutation_type = "double";
  }
  
  Set_New_Population_rim(common_variables);
}  

  
  
void Evolver::Run_Simulation(struct Common_Variables * common_variables){
  double halting_fitness = common_variables->halting_fitness;
  int halting_fixation = common_variables->halting_fixation;
  int fixation;
  int N_trajectory = common_variables->N_trajectory;
  
  if(common_variables->bl_input_filename)
    documents.Open_files(common_variables->input_filename,common_variables->bl_input_filename);
  else{
    documents.Open_files(common_variables->output_filename,common_variables->bl_input_filename);
    documents.Write_headers();
  }
  
  Record_Initial_State(common_variables);
  std::cout.flush();
  std::cout<<"[0/"<<N_trajectory<<" trajectories completed]";
  std::cout.flush();
  
  //rate dependent and masking is applied
  if(common_variables->rate && common_variables->mask){
    for(int trajectory=0; trajectory<N_trajectory;trajectory++){
      Get_Next_TransMach(trajectory,common_variables);
      Set_New_Population_rdm(common_variables);
      common_variables->prob_of_1_mutation=0;
      common_variables->prob_of_2_mutations=0;
      
      fixation = common_variables->end_fixation;
      if(!common_variables->bl_input_filename){
	common_variables->mutation_type="none";
	Record_Data(trajectory,fixation,common_variables);
      }
      
      while(population.fitness < halting_fitness && fixation < halting_fixation){
	Fix_rdm(common_variables);
	fixation++;
	Record_Data(trajectory,fixation,common_variables);
	
      }
      Record_Final_State(common_variables);
      std::cout<<"\r["<<trajectory+1<<"/"<<N_trajectory<<" trajectories completed]";
      std::cout.flush();
    }
  }
  //rate independent and masking is applied
  if(!common_variables->rate && common_variables->mask){
    for(int trajectory=0; trajectory<N_trajectory;trajectory++){
      Get_Next_TransMach(trajectory,common_variables);
      Set_New_Population_rim(common_variables);
      common_variables->prob_of_1_mutation=0;
      common_variables->prob_of_2_mutations=0;
      
      fixation = common_variables->end_fixation;
      if(!common_variables->bl_input_filename){
	common_variables->mutation_type="none";
	Record_Data(trajectory,fixation,common_variables);
      }
      
      while(population.fitness < halting_fitness && fixation < halting_fixation){
	Fix_rim(common_variables);
	fixation++;
	Record_Data(trajectory,fixation,common_variables);
	
      }
      Record_Final_State(common_variables);
      std::cout<<"\r["<<trajectory+1<<"/"<<N_trajectory<<" trajectories completed]";
      std::cout.flush();
    }
  }
  
  if(common_variables->rate && !common_variables->mask){
    for(int trajectory=0; trajectory<N_trajectory;trajectory++){
      Get_Next_TransMach(trajectory,common_variables);
      Set_New_Population_rdm(common_variables);
      common_variables->prob_of_1_mutation=0;
      common_variables->prob_of_2_mutations=0;
      
      fixation = common_variables->end_fixation;
      if(!common_variables->bl_input_filename){
	common_variables->mutation_type="none";
	Record_Data(trajectory,fixation,common_variables);
      }
      
      while(population.fitness < halting_fitness && fixation < halting_fixation){
	Fix_rdu(common_variables);
	fixation++;
	Record_Data(trajectory,fixation,common_variables);
	
      }
      Record_Final_State(common_variables);
      std::cout<<"\r["<<trajectory+1<<"/"<<N_trajectory<<" trajectories completed]";
      std::cout.flush();
    }
  }
  
  if(!common_variables->rate && !common_variables->mask){
    for(int trajectory=0; trajectory<N_trajectory;trajectory++){
      Get_Next_TransMach(trajectory,common_variables);
      Set_New_Population_rim(common_variables);
      common_variables->prob_of_1_mutation=0;
      common_variables->prob_of_2_mutations=0;
      
      fixation = common_variables->end_fixation;
      if(!common_variables->bl_input_filename){
	common_variables->mutation_type="none";
	Record_Data(trajectory,fixation,common_variables);
      }
      
      while(population.fitness < halting_fitness && fixation < halting_fixation){
	Fix_riu(common_variables);
	fixation++;
	Record_Data(trajectory,fixation,common_variables);
	
      }
      Record_Final_State(common_variables);
      std::cout<<"\r["<<trajectory+1<<"/"<<N_trajectory<<" trajectories completed]";
      std::cout.flush();
    }
  }

  std::cout<<std::endl;
  documents.Close_files();
}
 

