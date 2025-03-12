#ifndef READIFILE_HPP
#define READIFILE_HPP

#include <iostream>
#include <fstream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "BadConversion.hpp"
#include "convertToDouble.hpp"
#include "common_variables.hpp"


int Read_Input_File(struct Common_Variables * common_variables, bool bl_halting_fixation);

#endif
