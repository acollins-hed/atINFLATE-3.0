#include<iostream>
#include<fstream>

#ifndef DOCUMENT_HPP
#define DOCUMENT_HPP

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

#endif
