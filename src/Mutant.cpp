//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

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

//Mutant.cpp is part of atinflate_3
//atinflate_3 simulates the aaRS-tRNA Interaction Fitness LAndscape
//Topographer Express (atINFLATE) model, a dynamic extension of the
//atINFLAT model published in Collins-Hed and Ardell 2019
//https://www.sciencedirect.com/science/article/pii/S0040580918300789
//Copyright (C) 2026 Andrea Collins-Hed
//
//This file is part of atinflate_3
//
//atinflate_3 is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//atinflate_3 is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program. If not, see <https://www.gnu.org/licenses/>. 
