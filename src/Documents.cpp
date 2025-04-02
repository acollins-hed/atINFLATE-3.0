#include<iostream>
#include<fstream>

struct Documents{

  std::ofstream traj_file;
  std::ofstream prob_file;
  std::ofstream code_file;
  std::ofstream int_file;
  std::ofstream codon_file;
  std::ofstream ocheckpoint_file;

  void Open_files(std::string filename,bool input_file);
  void Close_files(void);  
  void Write_headers(void);
};


void Documents::Open_files(std::string filename,bool input_file){
  
  if(input_file){
    traj_file.open(filename+"_traj.dat",std::ios_base::app);
    code_file.open(filename+"_code.dat",std::ios_base::app);
    prob_file.open(filename+"_prob.dat",std::ios_base::app);
    int_file.open(filename+"_int.dat",std::ios_base::app);
    codon_file.open(filename+"_codon.dat",std::ios_base::app);
    ocheckpoint_file.open(filename+"_checkpoint.log");
  }
  else{
    traj_file.open(filename+"_traj.dat");
    code_file.open(filename+"_code.dat");
    prob_file.open(filename+"_prob.dat");
    int_file.open(filename+"_int.dat");
    codon_file.open(filename+"_codon.dat");
    ocheckpoint_file.open(filename+"_checkpoint.log");
  }
}

void Documents::Close_files(void){
  traj_file.close();
  prob_file.close();
  code_file.close();
  codon_file.close();
  int_file.close();
  ocheckpoint_file.close();
}
  
void Documents::Write_headers(void){
  traj_file<<"Trajectory Fixation Fitness Proportion_On Mutation_Type Mut.type.1 Mut.which.1 Mut.position.1 Mut.mask.1 Mut.type.2 Mut.which.2 Mut.position.2 Mut.mask.2\n";// P_Of_1Mut P_Of_2Mut\n";
  prob_file<<"Trajectory Fixation Site_Type Amino_Acid Probability\n";
  code_file<<"Trajectory Fixation tRNA aaRS Prob_Interaction Match kd\n";
  codon_file<<"Trajectory Fixation Site_Type Codon Codon_Frequency\n";
  int_file<<"Trajectory Fixation Molecule Number Type BValue DValue\n";
}
