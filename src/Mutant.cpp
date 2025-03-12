#include "Mutation.hpp"

struct Mutant{
  struct Mutation mutation_1;
  struct Mutation mutation_2;
  long double trans_prob;
  Mutant();
  Mutant(char tm, int which_one, int position, int mask);
  Mutant(char tm1, int wo1, int pos1, int mask1, char tm2, int wo2, int pos2, int mask2);
};

Mutant::Mutant(){
  mutation_1.kind_trans_machinery = '0';
  mutation_1.which_trans_machinery = -1;
  mutation_1.mutation_position = -1;
  mutation_1.mask = -1;
  mutation_2.kind_trans_machinery = '0';
  mutation_2.which_trans_machinery = -1;
  mutation_2.mutation_position = -1;
  mutation_2.mask = -1;
  trans_prob = 0;
}
Mutant::Mutant(char tm, int which_one, int position, int mask){
  mutation_1.kind_trans_machinery = tm;
  mutation_1.which_trans_machinery = which_one;
  mutation_1.mutation_position=position;
  mutation_1.mask = mask;
  mutation_2.kind_trans_machinery = '0';
  mutation_2.which_trans_machinery = -1;
  mutation_2.mutation_position = -1;
  mutation_2.mask = -1;
  trans_prob=0;
}
Mutant::Mutant(char tm1, int wo1, int pos1, int mask1, char tm2, int wo2, int pos2, int mask2){
  mutation_1.kind_trans_machinery = tm1;
  mutation_1.which_trans_machinery = wo1;
  mutation_1.mutation_position=pos1;
  mutation_1.mask = mask1;
  mutation_2.kind_trans_machinery = tm2;
  mutation_2.which_trans_machinery = wo2;
  mutation_2.mutation_position=pos2;
  mutation_2.mask = mask2;
  trans_prob=0;
}
