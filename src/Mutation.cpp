struct Mutation{
  int mutation_position;
  char kind_trans_machinery;
  int which_trans_machinery;
  int mask;
  Mutation();
};

Mutation::Mutation(){
  mutation_position=-1;
  kind_trans_machinery = '0';
  which_trans_machinery = -1;
  mask = -1;
}
