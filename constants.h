//
// Created by User on 4/28/2019.
//


#ifndef MATRIX_BUILDER_CONSTANTS_H
#define MATRIX_BUILDER_CONSTANTS_H
#include <thread>
#include <unordered_map>
#include <vector>
#include <bitset>
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>

#define SS2_WIDTH 114
#define SS2_HEIGHT 4599
typedef std::unordered_map<int,std::vector<std::string>> Field;
typedef std::unordered_map<std::string,std::vector<std::string>> Reverse_Field;
typedef std::mt19937 MyRNG;
class file_does_not_exist : std::exception {} ;
class file_line_num_bad : std::exception {};

#endif //MATRIX_BUILDER_CONSTANTS_H
