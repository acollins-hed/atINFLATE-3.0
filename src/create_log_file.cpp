#include<iostream>
#include<fstream>
#include "common_variables.hpp"

void create_log_file(struct Common_Variables * common_variables, bool bl_binom_p_0, bool bl_binom_p_1, bool bl_seed, unsigned int random_number_from_random_device){
  std::ofstream log_file(common_variables->output_filename+".log");
  time_t begin = time(0);
  tm *ltm = localtime(&begin);
  log_file<<"Started at "<<1900+ltm->tm_year<<"-"<<1+ltm->tm_mon<<"-"<<ltm->tm_mday<<" "<<ltm->tm_hour<<":"<<ltm->tm_min<<":"<<ltm->tm_sec<<std::endl;
  log_file<<"Saved in files: "<<common_variables->output_filename<<std::endl;
  std::cout<<"Copyleft (2025) A.I. Collins-Hed\nAll Wrongs Reversed.\nPlease cite Collins-Hed et al. 2019 in published works using this software.\n";
  std::cout<<"\n--------------------------------------\nn (number of sites per interface) is "<<common_variables->N_int_interface<<std::endl;
  log_file<<"Interaction interface width (bits): n = "<<common_variables->N_int_interface<<std::endl;
  //std::cout<<"k (number of sites for full energy) is "<<k<<std::endl;
  //log_file<<"Number of sites for full energy: k = "<<k<<std::endl;
  std::cout<<"T (number of tRNAs) is "<<common_variables->N_tRNA<<std::endl;
  log_file<<"Number of distinct tRNAs: T = "<<common_variables->N_tRNA<<std::endl;
  //std::cout<<"T starting (beginning number of tRNAs) is "<<Tcap<<std::endl;
  //log_file<<"Number of beginning tRNAs: Tcap = "<<Tcap<<std::endl;
  std::cout<<"A (number of aaRSs) is "<<common_variables->N_aaRS<<std::endl;
  log_file<<"Number of distinct aaRSs: A = "<<common_variables->N_aaRS<<std::endl;
  //std::cout<<"A starting (beginning number of aaRSs) is "<<Acap<<std::endl;
  //log_file<<"Number of beginning aaRSs: Acap = "<<Acap<<std::endl;
  std::cout<<"S (number of site-types) is "<<common_variables->N_site_type<<std::endl;
  log_file<<"Number of distinct site-types: S = "<<common_variables->N_site_type<<std::endl;
  std::cout<<"Site-Type frequencies are ";
  log_file<<"Site-Type frequencies: ";
  for(int i=0;i<common_variables->N_site_type;i++){
    std::cout<<common_variables->site_type_freqs(i);
    log_file<<common_variables->site_type_freqs(i);
    if(i != common_variables->N_site_type - 1)
      {
	log_file<<", ";
	std::cout<<",";
      }
  }
  std::cout<<std::endl;
  log_file<<std::endl;
  std::cout<<"The site-type physicochemical values are "<<common_variables->site_types(0);
  log_file<<common_variables->site_types(0);
  for(int i=1;i<common_variables->N_site_type;i++){
    std::cout<<", "<<common_variables->site_types(i);
    log_file<<", "<<common_variables->site_types(i);
  }
  std::cout<<std::endl;
  log_file<<std::endl;
  std::cout<<"N (population size) is "<<common_variables->N_population<<std::endl;
  log_file<<"Population size: N = "<<common_variables->N_population<<std::endl;
  if(common_variables->rate)
    {
      log_file<<"Rate dependent\n";
      std::cout<<"Rate dependent\n";
    }
  else
    {
      log_file<<"Rate independent\n";
      std::cout<<"Rate indepdent\n";
    }
  if(common_variables->mask)
    {
      log_file<<"Masking applied\n";
      std::cout<<"Masking applied\n";
    }
  else
    {
      log_file<<"No masking applied\n";
      std::cout<<"No masking applied\n";
    }
  if(common_variables->proofreading)
    {
      log_file<<"Proofreading applied\n";
      std::cout<<"Proofreading applied\n";
    }
  else
    {
      log_file<<"No proofreading\n";
      std::cout<<"No proofreading\n";
    }
  
  std::cout<<"The rate selection is "<<common_variables->rate_constant<<std::endl;
  log_file<<"Rate selection: "<<common_variables->rate_constant<<std::endl;
  
  if(bl_binom_p_0){
    log_file<<"Initialized at all 0 genotype.\n";
    std::cout<<"Initialized at all 0 genotype.\n";
  }
  else{
    if(bl_binom_p_1){
      log_file<<"Initialized at all 1 genotype.\n";
      std::cout<<"Initialized at all 1 genotype.\n";
    }
    else{
      log_file<<"Initialized at random genotype.\n";
      std::cout<<"Initialized at random genotype.\n";
      }
  }

  log_file<<"Mutation rate of identifying features and epsitatic modifiers: mu = "<<common_variables->mu_id_feat<<std::endl;
  std::cout<<"mu (mutation rate of identifying features and epistatic modifiers)is "<<common_variables->mu_id_feat<<std::endl;
  log_file<<"Codon substitution rate: Mu = "<<common_variables->mu_per_codon<<std::endl;
  std::cout<<"Mu (codon substitution rate) is "<<common_variables->mu_per_codon<<std::endl;
  log_file<<"Transition bias: "<<common_variables->transition_bias<<std::endl;
  std::cout<<"Transition bias: "<<common_variables->transition_bias<<std::endl;
  if(common_variables->codon_ring_space){
    log_file<<"Codon ring space.\n";
    std::cout<<"Codon ring space.\n";
  }
  if(common_variables->codonspace2){
    log_file<<"Codon space with 2 bases.\n";
    std::cout<<"Codon space with 2 bases.\n";
  }
  if(common_variables->codonspace4){
    log_file<<"Codon space with 4 bases.\n";
    std::cout<<"Codon space with 4 bases.\n";
  }

  log_file<<"Baseline missense fitness: phi = "<<common_variables->phi<<std::endl;
  std::cout<<"phi (the baseline fitness) is "<<common_variables->phi<<std::endl;
  log_file<<"The halting fitness: "<<common_variables->halting_fitness<<std::endl;
  std::cout<<"The halting fitness is "<<common_variables->halting_fitness<<std::endl;
  log_file<<"The halting fixation: "<<common_variables->halting_fixation<<std::endl;
  std::cout<<"The halting fixation is "<<common_variables->halting_fixation<<std::endl;
  log_file<<"The binomial parameter for genotype initialization: p = "<<common_variables->binom_p<<std::endl;
  std::cout<<"p for the binomial is "<<common_variables->binom_p<<std::endl;
  log_file<<"The number of trajectories: "<<common_variables->N_trajectory<<std::endl;
  std::cout<<"Number of trajectories is "<<common_variables->N_trajectory<<std::endl;
  log_file<<"The maximum dissociation rate: kmax = "<<common_variables->kmax<<std::endl;
  std::cout<<"kmax is "<<common_variables->kmax<<std::endl;
  log_file<<"The minimum dissociation rate: kmin = "<<common_variables->kmin<<std::endl;
  std::cout<<"kmin is "<<common_variables->kmin<<std::endl;
  //log_file<<"Fixations from the Moran Process.\n";
  //std::cout<<"Fixations from the Moran Process.\n";
  if(bl_seed){
    std::cout<<"PRNG seed is "<<common_variables->seed<<std::endl;
    log_file<<"PRNG seed: "<<common_variables->seed<<std::endl;
  }
  else{
    std::cout<<"PRNG seed is "<<random_number_from_random_device<<std::endl;
    log_file<<"PRNG seed: "<<random_number_from_random_device<<std::endl;
  }
  std::cout<<"Number of threads is "<<common_variables->N_threads<<std::endl;
  log_file<<"Number of threads: Nthread = "<<common_variables->N_threads<<std::endl;
  std::cout<<"Output file name is \""<<common_variables->output_filename<<"\"\n--------------------------------------\n\n";
  
  log_file.close();
}
