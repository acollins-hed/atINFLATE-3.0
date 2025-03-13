#ifndef CREATE_LOG_FILE_HPP
#define CREATE_LOG_FILE_HPP

#include<iostream>
#include<fstream>
#include "common_variables.hpp"

void create_log_file(struct Common_Variables * common_variables, bool bl_binom_p_0, bool bl_binom_p_1, bool bl_seed, unsigned int random_number_from_random_device);

#endif
