//
// Created by User on 4/27/2019.
//

#ifndef MATRIX_BUILDER_FIELD_AUX_H
#define MATRIX_BUILDER_FIELD_AUX_H
#include <thread>
#include <unordered_map>
#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>

#include "constants.h"

/*
*Generates the field into the file field.csv and to the data structure field
*@param field a field data structure monitoring field elements
*@return zero if successful and not otherwise
*/
int field_generator(Field& field);

/*
*Updates the field into the file field2.csv using the data structure field
*@param field a field data structure monitoring field elements and their orbits
*@return zero if successful and not otherwise
*/
int field_update(Field& field);

/*
*Marks the element of the field whose degree is deg as a member of the orbit param
*@param field a field data structure monitoring field elements and their orbits
*@param deg the degree identifying the field element
*@param orbit the degree identifying the orbit representative
*/
void markOrbit(Field& field, int deg, int orbit);



/*
*Checks the data structure field to see if the field element deg is marked in some orbit
*@param field a field data structure monitoring field elements and their orbits
*@param deg the degree identifying the field element
*@return boolean value that is true iff the degree is marked into some orbit in field
*/
bool checkMark(Field& field, int deg);

/*
*Checks the data structure field to see if any of the diffs fall in the same orbit
*@param field a field data structure monitoring field elements and their orbits
*@param diffs a vector of diffs which are regarded as degrees of field elements
*@return boolean value that is true iff the field elements whose degrees are the diffs fall in DIFFERENT orbits
*/
bool checkDifferences(Field& field, int diffs[]);


std::vector<std::vector<int>> orbits_file_to_SS2_matrix(Field& field);

void field2_to_Field(Field& field, Reverse_Field& reverse_field);

bool check_field2_integrity();

bool check_SS2_integrity();



#endif //MATRIX_BUILDER_FIELD_AUX_H
