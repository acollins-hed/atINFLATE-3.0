#include "common_variables.hpp"

long double Transition_Probability(double current_fitness, double mutant_fitness,int hamming, struct Common_Variables * common_variables){
  long double Fixation_Probability;
  if(current_fitness == mutant_fitness)
    Fixation_Probability = 1/((long double) common_variables->N_population);
  else
    Fixation_Probability = (long double) (1 - (current_fitness)/(mutant_fitness))/(1 - pow((current_fitness)/(mutant_fitness),common_variables->N_population));
  return ((long double) common_variables->N_population*pow(common_variables->mu_id_feat,hamming)*Fixation_Probability); 
}
