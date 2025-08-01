#include "Common_Variables.hpp"

int Get_Mutant_Vector_Size(struct Common_Variables * common_variables){
  int single_mutants;
  int double_mutants;

//Explanation for the calculation of total mutants:
  /*
  the number of possible single mutants is 2*(N_tRNA+N_aaRS)*N_int_interface 
  is the total number of mutants 1 Hamming distance away. The factor of 2
  comes from the fact that there are state bits and mask bits.
  
  the number of possible double mutants is single mutants choose 2
  */

  if(common_variables->mask)
    single_mutants = 2*(common_variables->N_tRNA+common_variables->N_aaRS)*common_variables->N_int_interface;
  else{
    single_mutants = (common_variables->N_tRNA+common_variables->N_aaRS)*common_variables->N_int_interface;
  }
  double_mutants = (single_mutants*(single_mutants-1))/2;
  
  if(common_variables->bl_double_mutants)
    return (double_mutants + single_mutants) ;
  else
    return single_mutants;
}
