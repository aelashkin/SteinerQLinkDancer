//
// Created by User on 4/27/2019.
//

#ifndef MATRIX_BUILDER_VEC_19_H
#define MATRIX_BUILDER_VEC_19_H

#include <thread>
#include <unordered_map>
#include <map>
#include <array>
#include <sstream>
#include <vector>
#include <bitset>
#include <fstream>
#include <iostream>
#include "constants.h"

//#define S_SIZE 524287


template <int N=19>
class vec_n
{
private:

    int degree;
    std::bitset<N> bit_vector;
public:

    static const int S_SIZE;

    /*
    *uses the file generated in the first part of the program and returns the bit vector representation of the field
    *element whose degree is deg
    *@param deg the degree of a field element
    *@return the bit vector representation of the field element whose degree is deg
    */
    std::bitset<N> deg_to_vector(int deg){

        std::ifstream myfile("./field.csv");
        // Open an existing file
        int deg2 = 0;
        // Read the Data from the file
        // as String Vector
        std::vector<std::string> row;
        std::string line, word, temp;
        getline(myfile, line);//skip over first line of the zero element, @
        while (!myfile.eof()) {

            row.clear();

            // read an entire row and
            // store it in a string variable 'line'
            getline(myfile, line);

            // used for breaking words
            std::stringstream s(line);
            // read every column data of a row and
            // store it in a string variable, 'word'
            while (getline(s, word, ',')) {

                // add all the column data
                // of a row to a vector
                row.push_back(word);
            }
            // convert string to integer for comparision

            deg2 = stoi(row[0]);

            if (deg == deg2) {
                //std::cout << row[1] << std::endl;
                return std::bitset<N>(row[1]);
            }
        }
        std::cout << "error: the index " << deg << " is not found " << std::endl;
        return std::bitset<N>(0);
    }
    /*
    *uses the file generated in the first part of the program and returns the degree representation of the field
    *element whose bit vector is vec
    *@param vec the bit vector of a field element
    *@return the degree representation of vec
    */
    int vector_to_deg(std::bitset<N> vec){
        std::ifstream myfile("./field.csv");
        // Open an existing file
        std::bitset<N> vec2 = 0;
        // Read the Data from the file
        // as String Vector
        std::vector<std::string> row;
        std::string line, word, temp;
        getline(myfile, line);//skip over first line of the zero element, @
        while (!myfile.eof()) {

            row.clear();


            // read an entire row and
            // store it in a string variable 'line'
            getline(myfile, line);

            // used for breaking words
            std::stringstream s(line);
            // read every column data of a row and
            // store it in a string variable, 'word'
            while (getline(s, word, ',')) {

                // add all the column data
                // of a row to a vector
                row.push_back(word);
            }
            // convert string to integer for comparision

            vec2 = std::bitset<N>(row[1]);

            if (vec == vec2) {
                //std::cout << row[0] << std::endl;
                return stoi(row[0]);
            }
        }
        std::cout << "error: the vector " << vec << " is not found " << std::endl;
        return -1;
    }

    //basic constructor for vector
    vec_n(int deg) {
        degree = deg;
        bit_vector = deg_to_vector(deg);
    }

    //constructs an iteration of frobenius map a^(deg * 2^power)
    vec_n(int deg, int power) {

        long long int tmp_deg = deg;
        long long int tmp_power = myPow(2, power);
        degree = (int  ((tmp_deg*tmp_power) % S_SIZE));
        bit_vector = deg_to_vector(degree);
    }

    vec_n& operator+=(const vec_n& rhs)
    {
        bit_vector ^= rhs.bit_vector;
        degree = vector_to_deg(bit_vector);
        return *this; // return the result by reference
    }

    friend vec_n operator+(vec_n lhs,        // passing lhs by value helps optimize chained a+b+c
                            const vec_n& rhs)
    {
        lhs += rhs;
        return lhs; // return the result by value
    }

    vec_n& operator*=(const vec_n& rhs)
    {
        degree = ((degree + rhs.degree) % S_SIZE);
        bit_vector = deg_to_vector(degree);
        return *this; // return the result by reference
    }

    friend vec_n operator*(vec_n lhs,        // passing lhs by value helps optimize chained a+b+c
                            const vec_n& rhs)
    {
        lhs *= rhs;
        return lhs; // return the result by value
    }

    //returns new vector a^(-degree)
    vec_n minus()
    {
        return vec_n(cyclic_substract(0, degree));
    }

    int get_deg() {
        return degree;
    }

    std::bitset<N> get_vector(){
        return bit_vector;
    }

    /*
    *Performs a version of exponentiation
    *@param x the base
    *@param p the power
    *@return an integer representing x to the power of p
    */
    int myPow (int x, int p) {
        int i = 1;
        for (int j = 1; j <= p; j++)  i *= x;
        return i;
    }

    /*
    *receives 2 positive integers and subtracts one from another under cyclic group
    *@param a first integer
    *@param b second integer
    *@return an integer representing the cyclic subtraction
    */
    int cyclic_substract(int a, int b) {
        int c = a - b;
        if (c >= 0)
            return c;
        else
            return S_SIZE + c;
    }


    /*
    *receives 2 positive integers and return their sum under cyclic group
    *@param a first integer
    *@param b second integer
    *@return an integer representing the cyclic addition
    */
    int cyclic_add(int a, int b) {
        return (a+b) % S_SIZE;
    }

};


template <int N> const int vec_n<N>::S_SIZE = pow(2,N)-1;


#endif //MATRIX_BUILDER_VEC_19_H
