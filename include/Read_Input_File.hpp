#ifndef READIFILE_HPP
#define READIFILE_HPP

#include <iostream>
#include <fstream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Bad_Conversion.hpp"
#include "Convert_To_Double.hpp"
#include "Common_Variables.hpp"


int Read_Checkpoint(struct Common_Variables * common_variables, bool bl_halting_fixation);

#endif
