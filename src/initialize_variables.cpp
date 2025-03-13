#include<iostream>
#include<fstream>
#include<random>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Read_Input_File.hpp"
#include "convertToDouble.hpp"
#include "BadConversion.hpp"
#include "common_variables.hpp"
#include "create_log_file.hpp"
#include<cfloat>
#include<unsupported/Eigen/KroneckerProduct>

int initialize_variables(int argc, char* argv[], struct Common_Variables * common_variables){

  bool bl_N_int_interface=false, bl_N_population=false, bl_mu_id_feat=false, bl_mu_per_codon=false,bl_binom_p=false, bl_N_tRNA=false, bl_N_aaRS=false, bl_N_trajectory=false, bl_transition_bias=false, bl_seed=false, bl_kmax=false, bl_kmin=false, bl_halting_fitness=false, bl_halting_fixation=false, bl_phi=false, bl_output_filename=false, bl_N_site_type=false,bl_binom_p_1=false,bl_binom_p_0=false,bl_rate_constant=false,bl_uniform_amino_acids=false,bl_stype_pcv=false;//bl_rate,bl_proofreading=false, bl_mask=false,bl_codon_space=false;
  
  std::string str_arg,sitestring,codonstring;

  std::uniform_real_distribution<double> aadist(0,std::nextafter(1,DBL_MAX));

  common_variables->codonspace2 = false;
  common_variables->codonspace4 = false;
  common_variables->codon_ring_space=true;
  common_variables->proofreading=false;
  common_variables->bl_input_filename=false;
  common_variables->output_filename="run";
  common_variables->end_fixation=0;
  common_variables->bl_double_mutants=true;
  //common_variables->prob_of_1_mutation = 0;
  //common_variables->prob_of_2_mutations = 0;
  common_variables->rate=true;
  common_variables->mask=true;
  common_variables->N_threads = 1;
  
  if(argc == 1)
    std::cout<<"\n\nEnter \"atinflate -h\" or \"atinflate --help\" for help\n\n";

  for(int i=1;i<argc;i++){
    str_arg = std::string(argv[i]);
    if(str_arg == "-h" || str_arg == "--help")
      {
	int help = system("echo -e 'ATINFLATE(1)			User Commands			ATINFLATE(1) \n \nNAME \n	atinflate - simulate aaRS-tRNA interaction network evolution. \n \nSYNOPSIS \n	atinflate [OPTIONS] \n \nDESCRIPTION \n	Initializes aaRS and tRNA interaction interfaces and evolve them through \n	a Moran Process \n	 \n	-0\n\t\tInitial genotype is all 0 genotype. Equivalent to\n\t\tatinflate --bp=0\n\t\tRandom by default.\n\n\t-1 \n\t\tInitial genotype is all 1 genotype. Equivalent to\n\t\tatinflate --bp=1\n\t\tRandom by default.\n \n\t-A=x, --aaRSs=x\n\t\twhere x is a positive integer greater than 1 representing the\n\t\tnumber of aaRSs. Default is 4.\n\n\t--bp=x\n\t\twhere x is in [0,1] representing the parameter p for a\n\t\tbinomial. Default is 0.5.\n\n\t--codon-space-x\n\t\twhere x is in {2,4} representing the number of bases.\n\t\tThis also switches from the default space which is a ring space\n\t\twhere each codon can only mutate to one of two neighbors\n\t\ton the ring.\n\t\tDefault is ring space.\n\n\t-f=x\n\t\tsets the rate selection parameter to x. Default: x = 1/44 for no\n\t\tproofreading and x = 1/9680 for proofreading.\n\n\t--halting-fitness=x \n\t\twhere x is in (0,1] representing the halting fitness. Default \n\t\tis 0.001.\n\n\t--halting-fixation=x\n\t\twhere x is in [1,infinity) representing the halting fixation.\n\t\tDefault is 1.\n\n\t-i=file_name, --ifile=file_name\n\t\twhere file_name is the name of the input checkpoint file\n\t\tDefault: none.\n\t\tOnly used with --codon-space-4\n\t\tDefault is 1.\n\n\t--kmax=x\n\t\twhere x is a positive float representing the maximum\n\t\tdissociation rate. Default is 10000\n\n\t--kmin=x \n\t\twhere x is a positive float representing the minimum\n\t\tdissociation rate. Default is 220\n\n\t--mu=x \n\t\twhere x is in (0,1) representing substitution rate in the\n\t\ttranslation machinery. Default is 1e-6.\n\n\t--Mu=x\n\t\twhere x is in (0,1) representing the mutation rate between\n\t\tcodons in the entire genome. Default is 1e-4.\n\n\t-n=x\n\t\twhere x is a positive integer representing interface size.\n\t\tDefault is 4.\n\n\t-N=x, --popsize=x\n\t\twhere x is a positive integer representing population size.\n\t\tDefault is 100\n\n\t--Nthread=x\n\t\tSets number of threads for multithreading. Default is 1.\n\n\t--no-double-mutants\n\t\tOnly checks single mutants. Checks up to double mutants by\n\t\tdefault.\n\n\t--no-mask\n\t\tTurns off masking. On by default.\n\n\t--no-rate\n\t\tMakes fitness rate independent. Fitness is rate dependent by\n\t\tdefault.\n\n\t--num-traj=x\n\t\twhere x is a positive integer representing the number of\n\t\tgenotypes to sample from the binomial for the run. Default is 1.\n\n\t-o=file_name, --ofile=file_name\n\t\twhere file_name is the name of the output file without \n\t\textension. Default is \"run.\"\n\n\t--phi=x\n\t\twhere x is in (0,1) representing missense tolerance, phi.\n\t\tDefault is 0.99.\n\n\t--proofreading\n\t\tapplies proofreading. Default: no proofreading.\n\n\t--seed=x\n\t\twhere x is the seed for the PRNG. Default is a random device.\n\n\t--trans-bias=x\n\t\twhere x is in [1,infinity) representing the transition bias.\n\t\tDefault is T.\n \n\t-S=x\n\t\twhere x is the number of site-types. Default is A\n\n\t--Site-Types=x,y,z,...\n\t\tto customize site-type frequencies. Example, if a system has 5\n\t\tsite-types and their frequencies are 3, 4, 1, 7, 9 \n\t\trespectively, enter\n\t\t\t\tatinflate -S=5 --Site-Types=3,4,1,7,9\n\t\tDefault: all site-type frequencies are set to 1\n\n\t-T=x, --tRNAs=x \n\t\twhere x is a positive integer greater than 1 representing the\n\t\tnumber of tRNAs. Default is 4.\n\n\t--uniform-amino-acids\n\t\tSets amino acid physicochemical values to uniform across the\n\t\tunit interval, e.g. if A=5, then the amino acid\n\t\tphysicochemical values are 0, 0.25, 0.5, 0.75, and 1.\n\t\tDefault: random across unit interval.\n \n \n \nAUTHOR \n	Written by Andrea Collins-Hed \n \nREPORTING BUGS \n	Report atinflate bugs to <codebugs.retrieval223@passinbox.com> \n \nCOPYLEFT \n	Copyleft (2025) A.I. Collins-Hed All Wrongs Reversed.Please cite \n	Collins-Hed 2025 in published works using this software.\n        Model originally developed in the Ardell lab and extended from the\n        atINFLAT model and atinflat.py published in Collins-Hed and Ardell\n        2019. Versions 1 and 2 of atinflate and the atINFLATE model were\n        also developed in the Ardell lab from 2020-2023.\n \n	  	      	February-March 2025	      	 	ATINFLATE(1) \n'| less");
	if(!help)
	return 0;
      }

    if(str_arg.substr(0,7) == "--seed=")
      {
	common_variables->seed = stoul(str_arg.substr(7));
	common_variables->mersene_twister.seed(common_variables->seed);
	bl_seed=true;
      }

    if(str_arg.substr(0,10) == "--Nthread=")
      {
	common_variables->N_threads = stoi(str_arg.substr(10));
	if(common_variables->N_threads <= 0){
	  std::cout<<"\nThe number of threads must be greater than 0. Use -h or --help for more.\n\n";
	  return 0;
	}
      }
    
    if(str_arg.substr(0,3) == "-n=")
      {
	common_variables->N_int_interface = stoi(str_arg.substr(3));
	if(common_variables->N_int_interface <= 0 ){
	  std::cout<<"\nn must be a positive integer. Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_N_int_interface=true;
      }

    if(str_arg.substr(0,3) == "-S=")
      {
	common_variables->N_site_type = stoi(str_arg.substr(3));
	if(common_variables->N_site_type < 2 ){
	  std::cout<<"\nS must be a positive integer > 2. Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_N_site_type=true;
      }

    
    if(str_arg.substr(0,11) == "--num-traj=")
      {
	common_variables->N_trajectory = stoi(str_arg.substr(11));
	if(common_variables->N_trajectory < 1 ){
	  std::cout<<"\nnum-traj must be a positive integer. Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_N_trajectory=true;
      }
    

    if(str_arg.substr(0,13) == "--trans-bias=")
      {
	common_variables->transition_bias = convertToDouble(str_arg.substr(13,str_arg.length()-1));
	if(common_variables->transition_bias < 1){
	  std::cout<<"\ntrans-bias must be a positive number, at least 1. Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_transition_bias=true;
      }

    
    if(str_arg.substr(0,3) == "-N=" || str_arg.substr(0,10) == "--popsize=")
      {
	if(str_arg.substr(0,3) == "-N=")
	  common_variables->N_population = stoi(str_arg.substr(3));
	else
	  common_variables->N_population = stoi(str_arg.substr(10));
	if(common_variables->N_population <= 0){
	  std::cout<<"\nN must be a positve integer. Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_N_population=true;
      }

    if(str_arg.substr(0,3) == "-A=" || str_arg.substr(0,8) == "--aaRSs=")
      {
	if(str_arg.substr(0,3) == "-A=")
	  common_variables->N_aaRS = stoi(str_arg.substr(3));
	else
	  common_variables->N_aaRS = stoi(str_arg.substr(8));
	if(common_variables->N_aaRS < 2)
	  {
	    std::cout<<"\nThe number of aaRSs must be an integer greater than 1. Use -h or --help for more.\n\n";
	    return 0;
	  }
	else
	  bl_N_aaRS=true;
      }

    if(str_arg.substr(0,3) == "-T=" || str_arg.substr(0,8) == "--tRNAs=")
      {
	if(str_arg.substr(0,3) == "-T=")
	  common_variables->N_tRNA = stoi(str_arg.substr(3));
	else
	  common_variables->N_tRNA = stoi(str_arg.substr(8));
	if(common_variables->N_tRNA < 1)
	  {
	    std::cout<<"\nThe number of tRNAs must be an integer greater than 1. Use -h or --help for more.\n\n";
	    return 0;
	  }
	else
	  bl_N_tRNA=true;
	  }
    /*
    if(str_arg.substr(0,13) == "--starting-A=")
      {
	Acap = stoi(str_arg.substr(13));
	if(Acap < 1)
	  {
	    std::cout<<"\nThe number of beginning aaRSs must be an integer greater than 0. Use -h or --help for more.\n\n";
	    return 0;
	  }
	else
	  baarsstart=true;
      }
    
    if(str_arg.substr(0,13) == "--starting-T=")
      {
	Tcap = stoi(str_arg.substr(13));
	if(Tcap < 1)
	  {
	    std::cout<<"\nThe number of beginning tRNAs must be an integer greater than 0. Use -h or --help for more.\n\n";
	    return 0;
	  }
	else
	  btrnasstart=true;
      }
    */	
    if(str_arg.substr(0,3) == "-o=" || str_arg.substr(0,8) == "--ofile=")
      {
	if(str_arg.substr(0,3) == "-o=")
	  common_variables->output_filename = str_arg.substr(3,str_arg.length()-1);
	else
	  common_variables->output_filename = str_arg.substr(8,str_arg.length()-1);
	if(common_variables->output_filename.empty()){
	  std::cout<<"\nDid you forget to give your output file a name? Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_output_filename=true;
      }

    if(str_arg.substr(0,3) == "-i=" || str_arg.substr(0,8) == "--ifile=")
      {
	if(str_arg.substr(0,3) == "-i=")
	  common_variables->input_filename = str_arg.substr(3,str_arg.length()-1);
	else
	  common_variables->input_filename = str_arg.substr(8,str_arg.length()-1);
	if(common_variables->input_filename.empty()){
	  std::cout<<"\nDid you forget to give your input file a name? Use -h or --help for more.\n\n";
	  return 0;
	}
	common_variables->bl_input_filename=true;
      }

    if(str_arg.substr(0,21) == "--uniform-amino-acids")
      {
	bl_uniform_amino_acids=true;
      }    

    if(str_arg.substr(0,19) == "--no-double-mutants")
      {
	common_variables->bl_double_mutants=false;
      }    
    
    if(str_arg.substr(0,9) == "--no-rate")
      {
	//bl_rate=false;
	common_variables->rate=false;
      }

    if(str_arg.substr(0,9) == "--no-mask")
      {
	//bl_rate=false;
	common_variables->mask=false;
      }

    if(str_arg.substr(0,14) == "--proofreading")
      {
	//bl_proofreading=true;
	common_variables->proofreading=true;
      }
    
    if(str_arg.substr(0,2) == "-1")
      {
	bl_binom_p_1=true;
	common_variables->binom_p = 1;
	bl_binom_p = true;
      }

    if(str_arg.substr(0,2) == "-0")
      {
	bl_binom_p_0=true;
	common_variables->binom_p = 0;
	bl_binom_p = true;
      }
    
    if(str_arg.substr(0,13) == "--Site-Types=")
      {
	sitestring = str_arg.substr(13,str_arg.length()-1)+",";
	bl_stype_pcv=true;
      }

    if(str_arg.substr(0,15) == "--codon-space-2")
      {
	common_variables->codonspace2=true;
	common_variables->codon_ring_space=false;
	if(common_variables->codonspace4){
	  std::cout<<"Can only have either 2 bases for codons or 4 but not both. Use --help or -h for more.\n\n";
	  return 0;
	}
      }

        if(str_arg.substr(0,15) == "--codon-space-4")
      {
	common_variables->codonspace4=true;
	common_variables->codon_ring_space=false;
	if(common_variables->codonspace2){
	  std::cout<<"Can only have either 2 bases for codons or 4 but not both. Use --help or -h for more.\n\n";
	  return 0;
	}
      }
    
    if(str_arg.substr(0,5) == "--mu=")
      {
	common_variables->mu_id_feat=convertToDouble(str_arg.substr(5,str_arg.length()-1));
	if(common_variables->mu_id_feat <= 0 || common_variables->mu_id_feat >= 1){
	  std::cout<<"\nmu must be in (0,1). Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_mu_id_feat=true;
      }

    if(str_arg.substr(0,5) == "--bp=")
      {
	common_variables->binom_p=convertToDouble(str_arg.substr(5,str_arg.length()-1));
	if(common_variables->binom_p < 0 || common_variables->binom_p > 1){
	  std::cout<<"\np must be in [0,1]. Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_binom_p=true;
      }

    
    if(str_arg.substr(0,5) == "--Mu=")
      {
	common_variables->mu_per_codon=convertToDouble(str_arg.substr(5,str_arg.length()-1));
	if(common_variables->mu_per_codon <= 0 || common_variables->mu_per_codon >= 1){
	  std::cout<<"\nMu must be in (0,1). Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_mu_per_codon=true;
      }
    
    if(str_arg.substr(0,18) == "--halting-fitness=")
      {
	common_variables->halting_fitness=convertToDouble(str_arg.substr(18,str_arg.length()-1));
	if(common_variables->halting_fitness <= 0 || common_variables->halting_fitness > 1){
	  std::cout<<"\nfthalt must be in (0,1]. Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_halting_fitness=true;
      }

    if(str_arg.substr(0,19) == "--halting-fixation=")
      {
	common_variables->halting_fixation=stoi(str_arg.substr(19));
	if(common_variables->halting_fixation < 1){
	  std::cout<<"\nhalting-fixation must be in [1, infinity). Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_halting_fixation=true;
      }


    if(str_arg.substr(0,7) == "--kmax=")
      {
	common_variables->kmax=float(convertToDouble(str_arg.substr(7,str_arg.length()-1)));
	if(common_variables->kmax <= 0){
	  std::cout<<"\nkmax must be in (0,infinity). Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_kmax=true;
      }

        if(str_arg.substr(0,7) == "--kmin=")
      {
	common_variables->kmin=float(convertToDouble(str_arg.substr(7,str_arg.length()-1)));
	if(common_variables->kmin <= 0){
	  std::cout<<"\nkmin must be in (0,infinity). Use -h or --help for more.\n\n";
	  return 0;
	}
	bl_kmin=true;
      }
	if(str_arg.substr(0,3) == "-f=")
	  {
	    common_variables->rate_constant=convertToDouble(str_arg.substr(3,str_arg.length()-1));
	    if(common_variables->rate_constant <= 0 ){
	      std::cout<<"\nf must be in (0,infinity). Use -h or --help for more.\n\n";
	      return 0;
	    }
	    bl_rate_constant=true;
	  }
	
        if(str_arg.substr(0,6) == "--phi=")
	  {
	    common_variables->phi=convertToDouble(str_arg.substr(6,str_arg.length()-1));
	    if(common_variables->phi <= 0 || common_variables->phi >= 1){
	      std::cout<<"\nphi must be in (0,1). Use -h or --help for more.\n\n";
	      return 0;
	    }
	    bl_phi=true;
	  }
	if(str_arg.substr(0,7) != "--seed=" && str_arg.substr(0,6) != "--phi=" && str_arg.substr(0,7) != "--kmin=" && str_arg.substr(0,7) != "--kmax=" && str_arg.substr(0,5) != "--mu=" && str_arg.substr(0,9) != "--no-rate" && str_arg.substr(0,21) != "--uniform-amino-acids" && str_arg.substr(0,2) != "-m" && str_arg.substr(0,3) != "-o=" && str_arg.substr(0,8) != "--ofile=" && str_arg.substr(0,3) != "-T=" && str_arg.substr(0,8) != "--tRNAs=" && str_arg.substr(0,3) != "-A=" && str_arg.substr(0,8) != "--aaRSs=" && str_arg.substr(0,3) != "-N=" && str_arg.substr(0,10) != "--popsize=" && str_arg.substr(0,3) != "-n=" && str_arg.substr(0,2) != "-1" && str_arg.substr(0,2) != "-0" && str_arg.substr(0,11) != "--num-traj=" && str_arg.substr(0,5) != "--bp=" && str_arg.substr(0,5) != "--Mu=" && str_arg.substr(0,13) != "--Site-Types=" && str_arg.substr(0,18) != "--halting-fitness="&& str_arg.substr(0,19) != "--halting-fixation="&& str_arg.substr(0,3) != "-S="&&str_arg.substr(0,15) != "--codon-space-4"&& str_arg.substr(0,15) != "--codon-space-2" && str_arg.substr(0,13) != "--trans-bias=" && str_arg.substr(0,3) != "-i=" && str_arg.substr(0,8) != "--ifile=" && str_arg.substr(0,14) != "--proofreading" && str_arg.substr(0,3) != "-f=" && str_arg.substr(0,19) != "--no-double-mutants" && str_arg.substr(0,9) != "--no-mask" && str_arg.substr(0,10) != "--Nthread="){
	  std::cout<<std::endl<<str_arg<<" is not a recognized parameter. Use -h or --help for more.\n\n";
	  return 0;
	}
  }  

  unsigned int random_number_from_random_device = common_variables->rand_dev();

  if(common_variables->bl_input_filename){
    if(!Read_Input_File(common_variables,bl_halting_fixation))
      return 0;
  }
  else{
    if(!bl_seed)
      common_variables->mersene_twister.seed(random_number_from_random_device);
    if(!bl_N_int_interface)
      common_variables->N_int_interface=4;
    //if(!bns)
    //  k=n;
    if(!bl_N_tRNA)
      common_variables->N_tRNA = 4;
    /*
    if(!btrnasstart)
      Tcap = T;
    if(Tcap > T){
      std::cout<<"The beginning number tRNAs cannot exceed the ending number of tRNAs. Use -h or --help for more.\n\n";
      return 0;
    }
    */
    if(!bl_N_aaRS)
      common_variables->N_aaRS = 4;
    /*
    if(!baarsstart)
      Acap = A;
    if(Acap > A){
      std::cout<<"The beginning number aaRSs cannot exceed the ending number of aaRSs. Use -h or --help for more.\n\n";
      return 0;
    }
    */
    if(!bl_N_site_type)
      common_variables->N_site_type = common_variables->N_aaRS;
    if(common_variables->N_site_type < common_variables->N_aaRS){
      std::cout<<"There cannot have more aaRSs (A) than site-types (S). See --help or -h for more.\n\n";
      return 0;
    }

    common_variables->site_types.resize(common_variables->N_site_type);
    std::cout<<std::endl;
    
    if(!bl_N_population)
      common_variables->N_population=100;
    if(!bl_mu_id_feat)
      common_variables->mu_id_feat=1e-6;
    if(!bl_mu_per_codon)
      common_variables->mu_per_codon=1e-4;
    if(!bl_transition_bias)
      common_variables->transition_bias=1;
    if(!bl_phi)
      common_variables->phi=0.99;
    if(!bl_rate_constant){
      if(common_variables->proofreading)
	common_variables->rate_constant = ((double) 1/9680);
      else
	common_variables->rate_constant = ((double) 1/44);
    }
    if(!bl_halting_fitness)
      common_variables->halting_fitness = 0.001;
    if(!bl_halting_fixation)
      common_variables->halting_fixation = 1;
    if(!bl_binom_p)
      common_variables->binom_p=0.5;
    if(!bl_N_trajectory)
      common_variables->N_trajectory = 1;
    common_variables->tRNA_State_bits.resize(common_variables->N_trajectory,common_variables->N_tRNA);
    common_variables->aaRS_State_bits.resize(common_variables->N_trajectory,common_variables->N_aaRS);
    common_variables->tRNA_Mask_bits.resize(common_variables->N_trajectory,common_variables->N_tRNA);
    common_variables->aaRS_Mask_bits.resize(common_variables->N_trajectory,common_variables->N_aaRS);
    if(!bl_output_filename)
      common_variables->output_filename="run";
    if(!bl_kmax)
      common_variables->kmax = 10000;
    if(!bl_kmin)
      common_variables->kmin = 220;
    if(common_variables->kmax <= common_variables->kmin){
      std::cout<<"\n\nkmax must be larger than kmin. Use -h or --help for more.\n\n";
      return 0;
    }

    /*
    if(k > n){
      std::cout<<"\nk must be at most n.\n\n";
      return 0;
    }
    */
    if(common_variables->codonspace2){
      common_variables->N_base = 2;
      common_variables->codon_ring_space=false;
      if(common_variables->N_tRNA != 2 && common_variables->N_tRNA!= 4 && common_variables->N_tRNA!=8){
	std::cout<<"\nThe number of tRNAs must be 2, 4, or 8 for a codon space with 2 bases. Use --help or -h for more.\n\n";
	return 0;
      }
    }
    if(common_variables->codonspace4){
      common_variables->N_base = 4;
      common_variables->codon_ring_space=false;
      if(common_variables->N_tRNA != 4 && common_variables->N_tRNA!= 16 && common_variables->N_tRNA!=64){
	std::cout<<"\nThe number of tRNAs must be 4, 16, or 64 for a codon space with 4 bases. Use --help or -h for more.\n\n";
	return 0;
      }
      
    }
    
    common_variables->site_type_freqs.resize(common_variables->N_site_type);
    for(int i=0;i<common_variables->N_site_type;i++){
      common_variables->site_type_freqs(i) = 1;
    }
    if(bl_stype_pcv){
      for(int i=0, j=0, stp=0;i<((int)sitestring.size());i++){
	if(stp >= common_variables->N_site_type){
	  std::cout<<"There are more fequencies given than site-types. Use --help or -h for more.\n\n";
	  return 0;
	}
	if(sitestring.at(i)==','){
	  common_variables->site_type_freqs(stp) = stoi(sitestring.substr(j,i-j));
	  stp++;
	  j = i + 1;
	  if(i == ((int)sitestring.size())-1 && stp < common_variables->N_site_type){
	    std::cout<<"WARNING: there are fewer fequencies given than site-types! Use --help or -h for more.\n\n";
	  }
	}
      }
      
    }

    if(bl_uniform_amino_acids){
      for(int i=0;i<common_variables->N_site_type;i++)
	common_variables->site_types(i) = ((double) i)/((double) common_variables->N_site_type-1);
    }
    else{
      for(int i=0;i<common_variables->N_site_type;i++)
	common_variables->site_types(i) = aadist(common_variables->mersene_twister);
    }
    
    {
      Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> pm(common_variables->N_site_type);
      pm.setIdentity();
      shuffle(pm.indices().data(),pm.indices().data()+pm.indices().size(),common_variables->mersene_twister);
      common_variables->site_types = pm*common_variables->site_types;
      common_variables->site_type_freqs = pm*common_variables->site_type_freqs;
    }
    common_variables->amino_acids.resize(common_variables->N_aaRS);
    for(int i=0;i<common_variables->N_aaRS;i++){
      common_variables->amino_acids(i) = common_variables->site_types(i);
      common_variables->aa_to_st[common_variables->amino_acids(i)] = i;
    }

  }
  
  create_log_file(common_variables,bl_binom_p_0,bl_binom_p_1,bl_seed,random_number_from_random_device);
  
  //L = (common_variables->N_tRNA+common_variables->N_aaRS)*2*common_variables->N_int_interface;

  
  common_variables->selection_mat.resize(common_variables->N_site_type,common_variables->N_site_type);
  common_variables->mutation_mat.resize(common_variables->N_tRNA,common_variables->N_tRNA);
  
  /*
  n_d_wout = (T+A)*(2*(T+A)-1)*n*n;
  n_d_win = (n*(n-1))*(T+A);
  */

  common_variables->epsilon = (log(common_variables->kmax)-log(common_variables->kmin))/double(common_variables->N_int_interface);
  if(common_variables->proofreading && !common_variables->bl_input_filename){
    common_variables->epsilon = 2*(log(common_variables->kmax)-log(common_variables->kmin))/double(common_variables->N_int_interface);
    common_variables->kmax *= common_variables->kmax;
    common_variables->kmin *= common_variables->kmin;
  }

  
  for(int i=0;i<common_variables->N_site_type;i++){
    for(int j=0;j<common_variables->N_site_type;j++){
      common_variables->selection_mat(i,j) = pow(common_variables->phi,abs(common_variables->site_types(i)-common_variables->site_types(j)));
    }
  }
  
  if(common_variables->codon_ring_space){
    for(int i = 0;i<common_variables->N_tRNA;i++){
      for(int j = 0; j<common_variables->N_tRNA;j++){
	if(i == j)
	  common_variables->mutation_mat(i,j) = 1 - 2*common_variables->mu_per_codon;
	else{
	  if(j == i+1 || j == i-1 || (j == 0 && i == common_variables->N_tRNA - 1) || (j == common_variables->N_tRNA - 1 && i == 0))
	    common_variables->mutation_mat(i,j) = common_variables->mu_per_codon;
	  else
	    common_variables->mutation_mat(i,j) = 0;
	}
      }
    }
  }
  else{
    if(common_variables->N_base == 2){
      Eigen::MatrixXd Codon_Mutation(2,2);
      for(int i=0;i<2;i++){
	for(int j=0;j<2;j++){
	  if(i == j)
	    Codon_Mutation(i,j) = 1-common_variables->mu_per_codon;
	  else
	    Codon_Mutation(i,j) = common_variables->mu_per_codon;
	}
      }

      if(common_variables->N_tRNA == 2)
	common_variables->mutation_mat = Codon_Mutation;
      if(common_variables->N_tRNA == 4){
	common_variables->mutation_mat = kroneckerProduct(Codon_Mutation,Codon_Mutation);
      }
      if(common_variables->N_tRNA == 8){
	common_variables->mutation_mat = kroneckerProduct(Codon_Mutation,kroneckerProduct(Codon_Mutation,Codon_Mutation));
      }
    }
    if(common_variables->N_base == 4){
      Eigen::MatrixXd Codon_Mutation(4,4);
      
      for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	  if(i == j)
	    Codon_Mutation(i,j) = 1-common_variables->mu_per_codon;
	  else{
	    if((i == 0 && j == 1)||(i == 1 && j==0) ||(i == 2 && j == 3)||(i == 3 && j == 2))
	      Codon_Mutation(i,j) = ((double) common_variables->transition_bias) * common_variables->mu_per_codon/((double)common_variables->transition_bias+2);
	    else
	      Codon_Mutation(i,j) = common_variables->mu_per_codon/((double) common_variables->transition_bias + 2);
	  }
	}
      }
      
      if(common_variables->N_tRNA == 4)
	common_variables->mutation_mat = Codon_Mutation;
      else{
	
	if(common_variables->N_tRNA == 16){
	  
	  common_variables->mutation_mat = kroneckerProduct(Codon_Mutation,Codon_Mutation);
	  
	}
	if(common_variables->N_tRNA == 64)
	  common_variables->mutation_mat = kroneckerProduct(Codon_Mutation,kroneckerProduct(Codon_Mutation,Codon_Mutation));
      }      
    }
  }
  
  return 1;
}
