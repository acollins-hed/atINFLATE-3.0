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
  std::vector<Mutant*> mutant_vector;
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

Evolver::Evolver(struct Population pop)
  :population(pop){
}

//int Hoarse_Partition(struct Mutant * mutant_array, int low){
int Evolver::Hoarse_Partition(int low, int high){
  int i = low-1, j = high+1;
  int mid = low + (high - low)/2;
  struct Mutant* pivot = mutant_vector[mid];
    
  while(true){
    do{
      i++;
    }while(mutant_vector[i]->trans_prob > pivot->trans_prob);
    
    do{
      j--;
    }while(mutant_vector[j]->trans_prob < pivot->trans_prob);
    
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
  //mutant_vector.resize(common_variables->N_total_mutants);
  //int N_int_interface = common_variables->N_int_interface;
  mutant_vector.clear();
  int N_sm = common_variables->N_single_mutants;
  
  int maxrow = 1;
  char tmach[2];
  tmach[0] = 't';
  tmach[1] = 'a';

  if(common_variables->mask)
    maxrow = 2;
  
  //Finding single mutation mutants
  for(int row=0; row<maxrow; row++){
    for(int col=0; col<population.genotype.trnas.iis.cols(); col++){
      for(int i=0;i<population.genotype.trnas.N_int_interface;i++){
	mutant_vector.push_back(new struct Mutant(tmach[0],col,i,row));
      }
    }
  }
  
  for(int row=0; row<maxrow; row++){
    for(int col=0; col<population.genotype.aarss.iis.cols(); col++){
      for(int i=0;i<population.genotype.aarss.N_int_interface;i++){
	mutant_vector.push_back(new struct Mutant(tmach[1],col,i,row));
      }
    }
  }

  if(common_variables->bl_double_mutants){
    for(int i=0;i<N_sm;i++){
      for(int j=i+1;j<N_sm;j++){
	char t1 = mutant_vector[i]->mutation_1.kind_trans_machinery;
	char t2 = mutant_vector[j]->mutation_1.kind_trans_machinery;
	int w1 = mutant_vector[i]->mutation_1.which_trans_machinery;
	int w2 = mutant_vector[j]->mutation_1.which_trans_machinery;
	int m1 = mutant_vector[i]->mutation_1.mask;
	int m2 = mutant_vector[j]->mutation_1.mask;
	int p1 = mutant_vector[i]->mutation_1.mutation_position;
	int p2 = mutant_vector[j]->mutation_1.mutation_position;
	mutant_vector.push_back(new struct Mutant(t1,w1,p1,m1,t2,w2,p2,m2));
      }
    }
  }
}

void Evolver::Get_Mutants_rd(struct Common_Variables * common_variables){
  long double sum=0;//for normalizing the sum of transition probabilities
  struct Mutant mutant;
  //common_variables->prob_of_1_mutation=0;
  //common_variables->prob_of_2_mutations=0;
  
  unsigned int N_threads = common_variables->N_threads;
  unsigned int mutants_per_thread = mutant_vector.size()/N_threads;
  
  std::vector<long double*> sum_vec;
  sum_vec.clear();
  std::vector<std::vector<struct Mutant*>> mut_subvec;
  std::vector<std::thread*> thrd;
  
  for(unsigned int i=0, start_val=0;start_val<mutant_vector.size();start_val+=mutants_per_thread,i++){
    sum_vec.push_back(new long double (0));
    unsigned int mutants_to_do = mutants_per_thread;
    if(start_val+mutants_per_thread < mutant_vector.size() && start_val + mutants_per_thread*2 > mutant_vector.size())
      mutants_to_do = ((int) mutant_vector.size()) - start_val;
    std::vector<struct Mutant*>::const_iterator first = mutant_vector.begin()+start_val;
    std::vector<struct Mutant*>::const_iterator last = mutant_vector.begin()+start_val + mutants_to_do;
    std::vector<struct Mutant*> v(first,last);
    mut_subvec.push_back(v);
    thrd.push_back(new std::thread(get_trans_prob_rd,mut_subvec[i],&population,common_variables,sum_vec[i]));
    if(mutants_to_do != mutants_per_thread)
      break;
  }
  
  for(unsigned int i=0;i<N_threads;i++)
    thrd[i]->join();

  std::for_each(sum_vec.begin(),sum_vec.end(), [&sum] (long double *subtotal) {sum += *subtotal;});

  for(unsigned int i=0;i<mutant_vector.size();i++)
    mutant_vector[i]->trans_prob /= sum;

  for(unsigned int i=0;i<N_threads;i++){
    delete thrd[i];
    delete sum_vec[i];
  }
  /*
  common_variables->prob_of_1_mutation = sum;

  common_variables->prob_of_1_mutation = common_variables->prob_of_1_mutation/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_2_mutations = common_variables->prob_of_2_mutations/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  */
  Quick_Sort(0,common_variables->N_total_mutants - 1);
}

void Evolver::Get_Mutants_ri(struct Common_Variables * common_variables){
  long double sum=0;//for normalizing the sum of transition probabilities
  struct Mutant mutant;
  
  //common_variables->prob_of_1_mutation=0;
  //common_variables->prob_of_2_mutations=0;

  unsigned int N_threads = common_variables->N_threads;
  unsigned int mutants_per_thread = mutant_vector.size()/N_threads;
  
  std::vector<long double*> sum_vec;
  sum_vec.clear();
  std::vector<std::vector<struct Mutant*>> mut_subvec;
  std::vector<std::thread*> thrd;
  
  for(unsigned int i=0, start_val=0;start_val<mutant_vector.size();start_val+=mutants_per_thread,i++){
    sum_vec.push_back(new long double (0));
    unsigned int mutants_to_do = mutants_per_thread;
    if(start_val+mutants_per_thread < mutant_vector.size() && start_val + mutants_per_thread*2 > mutant_vector.size())
      mutants_to_do = ((int) mutant_vector.size()) - start_val;
    std::vector<struct Mutant*>::const_iterator first = mutant_vector.begin()+start_val;
    std::vector<struct Mutant*>::const_iterator last = mutant_vector.begin()+start_val + mutants_to_do;
    std::vector<struct Mutant*> v(first,last);
    mut_subvec.push_back(v);
    thrd.push_back(new std::thread(get_trans_prob_ri,mut_subvec[i],&population,common_variables,sum_vec[i]));
    if(mutants_to_do != mutants_per_thread)
      break;
  }
  
  for(unsigned int i=0;i<N_threads;i++)
    thrd[i]->join();

  std::for_each(sum_vec.begin(),sum_vec.end(), [&sum] (long double *subtotal) {sum += *subtotal;});

  
  for(unsigned int i=0;i<mutant_vector.size();i++)
    mutant_vector[i]->trans_prob /= sum;
  

  for(unsigned int i=0;i<N_threads;i++){
    delete thrd[i];
    delete sum_vec[i];
  }
  /*
  common_variables->prob_of_1_mutation = sum;

  common_variables->prob_of_1_mutation = common_variables->prob_of_1_mutation/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  common_variables->prob_of_2_mutations = common_variables->prob_of_2_mutations/(common_variables->prob_of_1_mutation + common_variables->prob_of_2_mutations);
  */
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
  documents.ocheckpoint_file<<common_variables->N_threads<<std::endl;
}
  
void Evolver::Record_Data(int trajectory, int fixation, struct Common_Variables * common_variables){
  double frac_on=0;
  //trajectory file
  for(int i=0; i<common_variables->N_tRNA;i++)
    frac_on += __builtin_popcount(population.genotype.trnas.iis(1,i));
  for(int i=0;i<common_variables->N_aaRS;i++)
    frac_on += __builtin_popcount(population.genotype.aarss.iis(1,i));
  frac_on = frac_on/(common_variables->N_int_interface*(common_variables->N_tRNA+common_variables->N_aaRS));
  documents.traj_file<<trajectory<<" "<<fixation<<" "<<population.fitness<<" "<<frac_on<<" "<<common_variables->mutation_type<<" "<<common_variables->fixed_mutant.mutation_1.kind_trans_machinery<<" "<<common_variables->fixed_mutant.mutation_1.which_trans_machinery<<" "<<common_variables->fixed_mutant.mutation_1.mutation_position<<" "<<common_variables->fixed_mutant.mutation_1.mask<<" "<<common_variables->fixed_mutant.mutation_2.kind_trans_machinery<<" "<<common_variables->fixed_mutant.mutation_2.which_trans_machinery<<" "<<common_variables->fixed_mutant.mutation_2.mutation_position<<" "<<common_variables->fixed_mutant.mutation_2.mask<<std::endl;//" "<<common_variables->prob_of_1_mutation<<" "<<common_variables->prob_of_2_mutations<<std::endl;
  
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
    documents.int_file<<trajectory<<" "<<fixation<<" "<<"aaRS "<<common_variables->amino_acids(i)<<" State "<<population.genotype.aarss.get_aars_int(0,i)<<" "<<population.genotype.aarss.iis(0,i)<<std::endl;
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

void Evolver::Set_New_Population_rd(struct Common_Variables * common_variables){
  population.genotype.Get_Code();
  population.Get_Codon_Freq();
  population.fitness = Fitness_rate_dep(&population.codon_frequency,&population.genotype.code,&population.genotype.kd,common_variables);
}

void Evolver::Set_New_Population_ri(struct Common_Variables * common_variables){
  population.genotype.Get_Code();
  population.Get_Codon_Freq();
  population.fitness = Fitness_rate_indep(&population.codon_frequency,&population.genotype.code,common_variables);
}

void Evolver::Fix_rd(struct Common_Variables * common_variables){
  long double rand_n = common_variables->gillespie_rand(common_variables->mersene_twister), sum = 0;
  long unsigned int mutant_num=0;
  
  common_variables->mutation_type = "single";
  
  Get_Mutants_rd(common_variables);
  
  for(long unsigned int nth_mutant=0; nth_mutant<mutant_vector.size();nth_mutant++){
    sum += mutant_vector[nth_mutant]->trans_prob;
    if(sum > rand_n){
      mutant_num=nth_mutant;
      break;
    }
  }

  common_variables->fixed_mutant = *mutant_vector[mutant_num];
  
  if(mutant_vector[mutant_num]->mutation_1.kind_trans_machinery == 't'){
    population.genotype.trnas.iis(mutant_vector[mutant_num]->mutation_1.mask,mutant_vector[mutant_num]->mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num]->mutation_1.mutation_position;
  }
  else{
    population.genotype.aarss.iis(mutant_vector[mutant_num]->mutation_1.mask,mutant_vector[mutant_num]->mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num]->mutation_1.mutation_position;
  }
  
  if(mutant_vector[mutant_num]->mutation_2.kind_trans_machinery != '0'){
    if(mutant_vector[mutant_num]->mutation_2.kind_trans_machinery == 't'){
      population.genotype.trnas.iis(mutant_vector[mutant_num]->mutation_2.mask,mutant_vector[mutant_num]->mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num]->mutation_2.mutation_position;
    }
    else{
      population.genotype.aarss.iis(mutant_vector[mutant_num]->mutation_2.mask,mutant_vector[mutant_num]->mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num]->mutation_2.mutation_position;
    }
    common_variables->mutation_type = "double";
  }
  
  Set_New_Population_rd(common_variables);
}  

void Evolver::Fix_ri(struct Common_Variables * common_variables){
  long double rand_n = common_variables->gillespie_rand(common_variables->mersene_twister), sum = 0;
  long unsigned int mutant_num=0;
  
  common_variables->mutation_type = "single";
  
  Get_Mutants_ri(common_variables);
  
  for(long unsigned int nth_mutant=0; nth_mutant<mutant_vector.size();nth_mutant++){
    sum += mutant_vector[nth_mutant]->trans_prob;
    if(sum > rand_n){
      mutant_num=nth_mutant;
      break;
    }
  }

  common_variables->fixed_mutant = *mutant_vector[mutant_num];
  
  if(mutant_vector[mutant_num]->mutation_1.kind_trans_machinery == 't'){
    population.genotype.trnas.iis(mutant_vector[mutant_num]->mutation_1.mask,mutant_vector[mutant_num]->mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num]->mutation_1.mutation_position;
  }
  else{
    population.genotype.aarss.iis(mutant_vector[mutant_num]->mutation_1.mask,mutant_vector[mutant_num]->mutation_1.which_trans_machinery) ^= 1<<mutant_vector[mutant_num]->mutation_1.mutation_position;
  }
  
  if(mutant_vector[mutant_num]->mutation_2.kind_trans_machinery != '0'){
    if(mutant_vector[mutant_num]->mutation_2.kind_trans_machinery == 't'){
      population.genotype.trnas.iis(mutant_vector[mutant_num]->mutation_2.mask,mutant_vector[mutant_num]->mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num]->mutation_2.mutation_position;
    }
    else{
      population.genotype.aarss.iis(mutant_vector[mutant_num]->mutation_2.mask,mutant_vector[mutant_num]->mutation_2.which_trans_machinery) ^= 1<<mutant_vector[mutant_num]->mutation_2.mutation_position;
    }
    common_variables->mutation_type = "double";
  }
  
  Set_New_Population_ri(common_variables);
}

void Evolver::Run_Simulation(struct Common_Variables * common_variables){
  double halting_fitness = common_variables->halting_fitness;
  int halting_fixation = common_variables->halting_fixation;
  int fixation;
  int N_trajectory = common_variables->N_trajectory;
  std::vector<struct Mutant*> mut_vec_ptr(mutant_vector);
  
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
  if(common_variables->rate){
    for(int trajectory=0; trajectory<N_trajectory;trajectory++){
      Get_Next_TransMach(trajectory,common_variables);
      Set_New_Population_rd(common_variables);
      //common_variables->prob_of_1_mutation=0;
      //common_variables->prob_of_2_mutations=0;
      
      fixation = common_variables->end_fixation;
      if(!common_variables->bl_input_filename){
	common_variables->mutation_type="none";
	Record_Data(trajectory,fixation,common_variables);
      }
      
      while(population.fitness < halting_fitness && fixation < halting_fixation){
	Fix_rd(common_variables);
	fixation++;
	Record_Data(trajectory,fixation,common_variables);
	mutant_vector = mut_vec_ptr;
      }
      Record_Final_State(common_variables);
      std::cout<<"\r["<<trajectory+1<<"/"<<N_trajectory<<" trajectories completed]";
      std::cout.flush();
    }
  } else{
    for(int trajectory=0; trajectory<N_trajectory;trajectory++){
      Get_Next_TransMach(trajectory,common_variables);
      Set_New_Population_ri(common_variables);
      //common_variables->prob_of_1_mutation=0;
      //common_variables->prob_of_2_mutations=0;
      
      fixation = common_variables->end_fixation;
      if(!common_variables->bl_input_filename){
	common_variables->mutation_type="none";
	Record_Data(trajectory,fixation,common_variables);
      }
      
      while(population.fitness < halting_fitness && fixation < halting_fixation){
	Fix_ri(common_variables);
	fixation++;
	Record_Data(trajectory,fixation,common_variables);
	mutant_vector = mut_vec_ptr;
      }
      Record_Final_State(common_variables);
      std::cout<<"\r["<<trajectory+1<<"/"<<N_trajectory<<" trajectories completed]";
      std::cout.flush();
    }
  }

  for(unsigned int i=0;i<mutant_vector.size();i++){
    delete mutant_vector[i];
  }
  
  std::cout<<std::endl;
  documents.Close_files();
}
 

