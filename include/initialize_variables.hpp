#ifndef INITVAR_HPP
#define INITVAR_HPP

#include<iostream>
#include<fstream>
#include<random>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Read_Input_File.hpp"
#include "Convert_To_Double.hpp"
#include "Bad_Conversion.hpp"
#include "Common_Variables.hpp"
#include "create_log_file.hpp"
#include<cfloat>
#include<unsupported/Eigen/KroneckerProduct>

int initialize_variables(int argc, char* argv[], struct Common_Variables * common_variables);

#endif
