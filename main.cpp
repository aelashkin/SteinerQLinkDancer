//
//  main.cpp
//  Matrix builder
//
//  Created by Andrew Elashkin and Andy Berger on 1/1/19.
//  Copyright © 2019 NCE11 inc. & Davrosh LTD. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <bitset>
#include <vector>
#include <sstream>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>


#include <thread>
#include <functional>

#include "vec_n.h"
#include "field_aux.h"
#include "constants.h"
#include "dancing_links.h"
#include "matrix.h"


//code of GF(2^19)



/*
*Generates given number of 2 dimensional subspaces of a field under the normalizer of Singer subgroup
*@param orb_num number of orbits to generate, max value is 4599
*@return matrix of type vector<vector<int>> that holds all generated orbits
*/
std::vector<std::vector<int>> SS2_generator(Field& field,int orb_num)
{
    vec_n<> helper(0);
    int i = 1;
    std::vector<std::vector<int>> SS2_matrix;
    if (orb_num > 4599 || orb_num < 1) {
        throw "orb_num is out of borders";
    }
    
    
    while(SS2_matrix.size() < orb_num) {
        vec_n<> a_0(0);
        std::cout << "current i is " << i << "\n";
        std::cout << "current orbit is " << SS2_matrix.size() << "\n";
        vec_n<> a_i(i);
        if (checkMark(field,a_i.get_deg()) == true) {
            std::cout << i <<  " is already in some orbit" << "\n";
            i++;
            continue;
        }
        
        vec_n<> a_j = a_0 + a_i;         //calculating a_j and all 6 differences
        int diff_array[6] = {a_i.get_deg(), a_j.get_deg(),helper.cyclic_substract(a_i.get_deg(), a_j.get_deg()), helper.cyclic_substract(a_j.get_deg(), a_i.get_deg()), a_i.minus().get_deg(), a_j.minus().get_deg()};
        
        std::vector<int> SS2_orbit_i;
        for (int j = 0; j<6; j++) {
            long long int tmp_deg = diff_array[j];
            long long int tmp_power = 1;
            for(int k = 0; k<19; k++) {
                int degree = (int  ((tmp_deg*tmp_power) % vec_n<>::S_SIZE));
                //vec_n<> new_vec(diff_array[j], k);
                SS2_orbit_i.push_back(degree);
                markOrbit(field,degree, i);
                tmp_power = tmp_power * 2;
            }
        }
        SS2_matrix.push_back(SS2_orbit_i);
        i++;
    }
    return SS2_matrix;
}

/*
 *Generates given number of 3 dimentional subspaces of GF(19,2)
 *@param thread_num number of threads to run
 *@param ss_num number of 3D SS to generate
 *@return matrix of type int where i is number of orbit representing second vector and j is power of vector j
 */

int SS3_generator(Field& field, std::vector<std::vector<int>> SS2_matrix, int tread_num, int ss_num) {

    vec_n<> helper(0);
    vec_n<> a_0(0); //first element of any 3D SS
    srand((unsigned int) time(nullptr));
    std::vector<std::vector<int>> SS3_matrix;

    while (SS3_matrix.size() <= ss_num) {
        for (int i = 0; i < 4599; i++) { //Second element of 3D SS, representative of an orbit

            //do the threads here for each i, rotate them when the last i is done;

            vec_n<> a_i(SS2_matrix[i][0]); //vector i is one representative of the orbit i
            for(int j = 2*i; j < 10000; j++) { //j is at least twice bigger than i

                vec_n<> a_j(rand());

                int diff_array[7] = {a_i.get_deg(), a_j.get_deg(),helper.cyclic_substract(a_i.get_deg(), a_j.get_deg()), helper.cyclic_substract(a_j.get_deg(), a_i.get_deg()), a_i.minus().get_deg(), a_j.minus().get_deg(), helper.cyclic_substract(a_j.minus().get_deg(), a_i.get_deg())};

                //add func that checks if several vectors are in the same orbit and run it from diffs
                SS3_matrix[i][j] = a_j.get_deg();
                //add it to the algo x table

            }
        }
        break;
    }
    return 1;
}

/*
 *Generates orbits of 3D SS's and checks the number of such SS's that have all 7 of their 2 SS's are in different orbits.
 *@param field - a DS of the field as a map: from the degree to a vector containing the bit representation and the orbit degree (if it exists)
 *@param reverse_field - a DS of the field as a map: from the bit representation to a vector containing the degree and the orbit degree (if it exists)
 *@param SS2_matrix - a DS of the SS2 orbits as a matrix: each row i contains all 2D SS's of the i'th orbit (the first index is the orbit representative).
 *@return
 */
int SS3_generator_aux(Field& field,Reverse_Field& reverse_field, std::vector<std::vector<int>> SS2_matrix, int OrbNumber) {

    std::bitset<19> zero_bit = std::bitset<19>(field[0][0]);// bit rep. of a^0


    int i = OrbNumber;//pick some orbit.


    //Generate a file with an index to write the orbits to.
    std::string name = "SS3_orbits_";
    name+= std::to_string(i);
    name+= ".txt";


    std::fstream myfile(name, std::ios_base::binary|std::fstream::out);



    if(myfile.is_open() == 0)
        return 1;





    int deg = SS2_matrix[i][0];//deg is the orbit rep. of the i'th orbit: what we will refer to as a^i, though it's actually a^deg
    std::bitset<19> deg_bit = std::bitset<19>(field[deg][0]);//bit rep. of a^i
    int bad_j = std::stoi(reverse_field[(zero_bit^deg_bit).to_string()][0]);//a^j = a^0+a^i, this j is not independent of 0 and i


    for(int j = deg + 1; j < vec_n<>::S_SIZE; j++){


        if(j == bad_j){
            continue;
        }

        myfile << deg << "," << j << std::endl;//write row to file

    }
    myfile.close();
    std::cout << "Generated orbit n'" << i << std::endl;


    return 0;
}


int reduced_kramer_generator_aux(Field& field,Reverse_Field& reverse_field, std::vector<std::vector<int>> SS2_matrix, char cols[SS2_HEIGHT], int OrbNumber){
    int i = OrbNumber;//pick some orbit.
    vec_n<> helper(0);
    std::bitset<19> zero_bit = std::bitset<19>(field[0][0]);// bit rep. of a^0
    std::vector<std::vector<int>> SS3_matrix;
    int deg = SS2_matrix[i][0];//deg is the orbit rep. of the i'th orbit: what we will refer to as a^i, though it's actually a^deg
    std::bitset<19> deg_bit = std::bitset<19>(field[deg][0]);//bit rep. of a^i

    //Generate a file with an index to write the orbits to.
    std::ofstream kramer;
    std::string name = "kramer_mesner_";
    name+= std::to_string(i);
    name+= ".txt";
    kramer.open (name);
    if(kramer.is_open() == 0)
        return 1;


    for(int j = 2 * deg; j < vec_n<>::S_SIZE; j++){
    //for(int j = deg + 1; j < 2 * deg; j++){


        /*
         * The seven 2D SS's of a^0,a^i,a^j are:
         * (a^0,a^i) => i (actually deg) is the deg of the SS.
         * (a^0,a^j) => j
         * (a^i,a^j) == (a^0, a^(j-i)) => j-deg
         * (a^0,a^i+a^j) == (a^0,=a^k) => k
         * (a^i,a^0+a^j) == (a^i,a^u) == (a^0,a^(u-i)) => u-deg
         * (a^j,a^0+a^i) == (a^j,a^v) == (a^0,a^(v-j)) => v-j
         * (a^0+a^i,a^0+a^j) == (a^v,a^u) == (a^0,a^(u-v)) => u-v
         *
         *
         * */


        int bad_j = std::stoi(reverse_field[(zero_bit^deg_bit).to_string()][0]);//a^j = a^0+a^i, this j is not independent of 0 and i
        if(j == bad_j){
            continue;
        }

        std::bitset<19> j_bit = std::bitset<19>(field[j][0]);//bit rep. of a^j
        int k = std::stoi(reverse_field[(deg_bit^j_bit).to_string()][0]);//a^k = a^i+a^j, gets degree k
        int u = std::stoi(reverse_field[(zero_bit^j_bit).to_string()][0]);//a^u = a^0+a^j, gets degree u
        int v = std::stoi(reverse_field[(zero_bit^deg_bit).to_string()][0]);//a^v = a^0+a^i, gets degree v

        //calculate the 7 degrees.
        int diff_array[7] = {deg, j, helper.cyclic_substract(j, deg), k, helper.cyclic_substract(u,deg), helper.cyclic_substract(v, j), helper.cyclic_substract(u,v)};

        //if all 2D SS's fall in different orbits then count this 3D orbit.
        if(checkDifferences(field,diff_array)){


            int arr[SS2_HEIGHT];

            for(int orb = 0; orb < SS2_HEIGHT; orb++){
                arr[orb] = 0;
            }

            bool already_covered = false;
            for(int diff = 0; diff < 7; diff++){
                int index = std::stoi(field[diff_array[diff]][2]);
                if(cols[index] == 1){
                    already_covered = true;
                    break;
                }
                arr[index] = 1;

            }
            if(already_covered){
                continue;
            }


            kramer << deg << "." << j << ",";

            for(int orb = 0; orb < SS2_HEIGHT; orb++){
                if(cols[orb] == 0){
                    kramer << arr[orb] << ",";
                }
            }

            kramer << std::endl;



        }

    }
    kramer.close();
    std::cout << "Generated reduced kramer file n' " << i << std::endl;
    return 0;
}


int kramer_generator_aux(Field& field,Reverse_Field& reverse_field, std::vector<std::vector<int>> SS2_matrix, int OrbNumber) {
    vec_n<> helper(0);
    std::bitset<19> zero_bit = std::bitset<19>(field[0][0]);// bit rep. of a^0
    std::vector<std::vector<int>> SS3_matrix;
    int i = OrbNumber;//pick some orbit.
    int deg = SS2_matrix[i][0];//deg is the orbit rep. of the i'th orbit: what we will refer to as a^i, though it's actually a^deg
    std::bitset<19> deg_bit = std::bitset<19>(field[deg][0]);//bit rep. of a^i




    //Generate a partial Kramer Mesner matrix
    std::ofstream kramer;
    std::string kramer_name = "kramer_mesner_";
    kramer_name+= std::to_string(i);
    kramer_name+= ".txt";
    kramer.open (kramer_name);
    if(kramer.is_open() == 0)
        return 1;

    int rows_taken = 0;

    std::vector<int> taken;

    //TODO: turn back into 8
    while(rows_taken < 8){

        //generate random number
        std::random_device dev;
        MyRNG rng(dev());
        std::uniform_int_distribution<MyRNG::result_type> dist(2 * deg, vec_n<>::S_SIZE - 1);
        int j = dist(rng);


        //don't want to consider a good j more than once, all 8 rows must be unique.
        for (const auto &item : taken) {
            if(j == item){
                continue;
            }
        }



        /*
         * The seven 2D SS's of a^0,a^i,a^j are:
         * (a^0,a^i) => i (actually deg) is the deg of the SS.
         * (a^0,a^j) => j
         * (a^i,a^j) == (a^0, a^(j-i)) => j-deg
         * (a^0,a^i+a^j) == (a^0,=a^k) => k
         * (a^i,a^0+a^j) == (a^i,a^u) == (a^0,a^(u-i)) => u-deg
         * (a^j,a^0+a^i) == (a^j,a^v) == (a^0,a^(v-j)) => v-j
         * (a^0+a^i,a^0+a^j) == (a^v,a^u) == (a^0,a^(u-v)) => u-v
         *
         *
         * */

        std::bitset<19> j_bit = std::bitset<19>(field[j][0]);//bit rep. of a^j
        int k = std::stoi(reverse_field[(deg_bit^j_bit).to_string()][0]);//a^k = a^i+a^j, gets degree k
        int u = std::stoi(reverse_field[(zero_bit^j_bit).to_string()][0]);//a^u = a^0+a^j, gets degree u
        int v = std::stoi(reverse_field[(zero_bit^deg_bit).to_string()][0]);//a^v = a^0+a^i, gets degree v

        //calculate the 7 degrees.
        int diff_array[7] = {deg, j, helper.cyclic_substract(j, deg), k, helper.cyclic_substract(u,deg), helper.cyclic_substract(v, j), helper.cyclic_substract(u,v)};

        //if all 2D SS's fall in different orbits then count this 3D orbit.
        if(checkDifferences(field,diff_array)){


            int arr[SS2_HEIGHT];

            for(int orb = 0; orb < SS2_HEIGHT; orb++){
                arr[orb] = 0;
            }

            for(int diff = 0; diff < 7; diff++){
                arr[std::stoi(field[diff_array[diff]][2])] = 1;
            }

            kramer << deg << "." << j << ",";

            for(int orb = 0; orb < SS2_HEIGHT; orb++){
                kramer << arr[orb] << ",";
            }

            kramer << std::endl;

            rows_taken++;
            taken.push_back(j);

        }
    }
    kramer.close();

    std::cout << "Generated kramer file n'" << i << std::endl;

    return 0;
}


int matrix_generator_aux(Field& field,Reverse_Field& reverse_field, std::vector<std::vector<int>> SS2_matrix,matrix& m, int OrbNumber) {
    vec_n<> helper(0);
    std::bitset<19> zero_bit = std::bitset<19>(field[0][0]);// bit rep. of a^0
    std::vector<std::vector<int>> SS3_matrix;
    int i = OrbNumber;//pick some orbit.
    int deg = SS2_matrix[i][0];//deg is the orbit rep. of the i'th orbit: what we will refer to as a^i, though it's actually a^deg
    std::bitset<19> deg_bit = std::bitset<19>(field[deg][0]);//bit rep. of a^i


    int rows_taken = 0;

    std::vector<int> taken;

    //TODO: turn back into 8
    while(rows_taken < 8){

        //generate random number
        std::random_device dev;
        MyRNG rng(dev());
        std::uniform_int_distribution<MyRNG::result_type> dist(2 * deg, vec_n<>::S_SIZE - 1);
        int j = dist(rng);


        //don't want to consider a good j more than once, all 8 rows must be unique.
        for (const auto &item : taken) {
            if(j == item){
                continue;
            }
        }



        /*
         * The seven 2D SS's of a^0,a^i,a^j are:
         * (a^0,a^i) => i (actually deg) is the deg of the SS.
         * (a^0,a^j) => j
         * (a^i,a^j) == (a^0, a^(j-i)) => j-deg
         * (a^0,a^i+a^j) == (a^0,=a^k) => k
         * (a^i,a^0+a^j) == (a^i,a^u) == (a^0,a^(u-i)) => u-deg
         * (a^j,a^0+a^i) == (a^j,a^v) == (a^0,a^(v-j)) => v-j
         * (a^0+a^i,a^0+a^j) == (a^v,a^u) == (a^0,a^(u-v)) => u-v
         *
         *
         * */

        std::bitset<19> j_bit = std::bitset<19>(field[j][0]);//bit rep. of a^j
        int k = std::stoi(reverse_field[(deg_bit^j_bit).to_string()][0]);//a^k = a^i+a^j, gets degree k
        int u = std::stoi(reverse_field[(zero_bit^j_bit).to_string()][0]);//a^u = a^0+a^j, gets degree u
        int v = std::stoi(reverse_field[(zero_bit^deg_bit).to_string()][0]);//a^v = a^0+a^i, gets degree v

        //calculate the 7 degrees.
        int diff_array[7] = {deg, j, helper.cyclic_substract(j, deg), k, helper.cyclic_substract(u,deg), helper.cyclic_substract(v, j), helper.cyclic_substract(u,v)};

        //if all 2D SS's fall in different orbits then count this 3D orbit.
        if(checkDifferences(field,diff_array)){


            int arr[7];


            for(int col = 0; col < 7; col++){
                arr[col] = std::stoi(field[diff_array[col]][2]);
            }

            std::string indices = std::to_string(deg);
            indices += ".";
            indices += std::to_string(j);

            m.add_row(arr, indices);


            rows_taken++;
            taken.push_back(j);

        }
    }


    std::cout << "Taken rows from orbit n'" << i << std::endl;

    return 0;
}

class GeneratedData
{
private:
    Field field;
    Reverse_Field reverse_field;
    std::vector<std::vector<int>> orbits_matrix;
    char cols[SS2_HEIGHT];
public:
    GeneratedData(Field inp_field, Reverse_Field inp_reverse_field, std::vector<std::vector<int>> matrix, char inp_cols[SS2_HEIGHT] = nullptr) {
        field = inp_field;
        reverse_field = inp_reverse_field;
        orbits_matrix = matrix;
        if(inp_cols != nullptr){
            for(int i = 0; i < SS2_HEIGHT; i++){
                cols[i] = inp_cols[i];
            }
        }
    }
    Field& getField() {
        return field;
    }
    Reverse_Field& getReverse() {
        return reverse_field;
    }
    std::vector<std::vector<int>>& getOrbits() {
        return orbits_matrix;
    }

    char *getCols()  {
        return cols;
    }


};

void threadCalcOrbit(GeneratedData params, int currentThread, int totalThreads, std::function<int(Field&,Reverse_Field&, std::vector<std::vector<int>>, int)>& f) {
    int currentOrbit = currentThread;
    while (currentOrbit < SS2_HEIGHT) //SS2_HEIGHT = 4599
    {
        f(params.getField(), params.getReverse(), params.getOrbits(), currentOrbit);
        currentOrbit += totalThreads;
    }
}

void threadCalcReduced(GeneratedData params, int currentThread, int totalThreads){
    //, std::function<int(Field&,Reverse_Field&, std::vector<std::vector<int>>,char[SS2_HEIGHT], int)>& g
    int currentOrbit = currentThread;
    while (currentOrbit < SS2_HEIGHT) //SS2_HEIGHT = 4599
    {
//        g(params.getField(), params.getReverse(), params.getOrbits(),params.getCols(), currentOrbit);
        reduced_kramer_generator_aux(params.getField(), params.getReverse(), params.getOrbits(),params.getCols(), currentOrbit);
        currentOrbit += totalThreads;
    }
}



int part1(Field& field, Reverse_Field& reverse_field){

    std::cout << "You will now generate the GF(2,19) field and SS2 orbits. \n"
            "Warning: file field2.csv weighs roughly 16.76MB and SS2_orbits.txt roughly 3.56MB. \n"
            "Please make sure you have the necessary space on your machine, otherwise the program has undefined behavior! \n"
            "The generation of both files takes between 2 to 3 HOURS, just so you know :) \n"
            "If for any reason you wish to quit now input 'q',if you wish to start input 's':" << std::endl;

    std::string prompt;
    std::cin >> prompt;

    while(prompt != "s" && prompt != "q"){
        std::cout << "Please either input 's' to start or input 'q' to quit:" << std::endl;
        std::cin >> prompt;
    }
    if(prompt == "q"){
        std::cout << "Quitting, have a nice day!" << std::endl;
        return 0;
    }

    //starting...

    int res = field_generator(field);
    if (res != 0)
        std::cout << "Error!\n";

    //number of orbits to generate
    int num = 4599;
    std::vector<std::vector<int>> SS2_matrix = SS2_generator(field,num);
    std::cout << num << " orbits of SS2 Generated\n";

    std::ofstream outFile("SS2_orbits.txt");
    for (int i = 0 ; i <  SS2_matrix.size(); i++) {
        for (const auto &e : SS2_matrix[i]) outFile << e << ",";
        outFile << "\n";
    }
    std::cout << "SS2 is successfully saved\n";

    int res2 = field_update(field);
    if(res2 != 0){
        std::cout << "Error!\n";
        return 1;
    }

    return 0;
}

int merge_files(std::string filename, char cols[SS2_HEIGHT] = nullptr){


    std::string name = filename;
    name+= ".txt";
    std::ofstream myfile(name, std::ios_base::binary);
    if(myfile.is_open() == 0)
        return 1;

    if(filename == "kramer_mesner"){
        myfile << " ,";

        for(int i = 0; i < SS2_HEIGHT; i++){
            if(cols == nullptr) {
                myfile << i << ",";
            }
            else{
                if(cols[i] == 0){
                    myfile << i << ",";
                }
            }
        }

        myfile << std::endl;
    }


    for(int i=0; i< SS2_HEIGHT; i++){
        std::string curr_name = filename;
        curr_name += "_";
        curr_name += std::to_string(i);
        curr_name += ".txt";

        std::fstream curr (curr_name, std::ios_base::in|std::ios_base::out|std::ios_base::binary);

        if(!curr.is_open()){
            return 1;
        }

        if(curr.rdbuf()->in_avail()){
            myfile << curr.rdbuf();
        }


        curr.close();
        while(!remove(curr_name.c_str())){
            std::cout << "Trying to delete file n' " << i << ". If this message persists then the file cannot be deleted, make sure it isn't opened by any program." << std::endl;
//            std::cout << "Help! Get me out of here. Abort! Abort!" << std::endl;
        }
    }
    myfile.close();



    return 0;

}

int part2_checks(Field& field, Reverse_Field& reverse_field, std::vector<std::vector<int>>& file_generated_matrix){
    std::cout << "checking file field2.csv" << std::endl;

    if(!check_field2_integrity()){
        std::cout << "check of field2.csv failed" <<std::endl;
        return 1;
    }

    std::cout << "file field2.csv exists and good to go!" << std::endl;

    std::cout << "checking file SS2_orbits.txt" << std::endl;

    if(!check_SS2_integrity()){
        std::cout << "check of SS2_orbits.txt failed" <<std::endl;
        return 1;
    }

    std::cout << "file SS2_orbits.txt exists and good to go!" << std::endl;


    std::cout << "transferring field data from field2.csv" << std::endl;

    field2_to_Field(field, reverse_field);

    std::cout << "field transfer complete!" << std::endl;

    std::cout << "transferring SS2 orbit data from SS2_orbits.txt" << std::endl;

    file_generated_matrix = orbits_file_to_SS2_matrix(field);

    std::cout << "orbit transfer complete!" << std::endl;
    return 0;
}

int part2_threads(Field& field, Reverse_Field& reverse_field, std::vector<std::vector<int>>& file_generated_matrix, std::function<int(Field&,Reverse_Field&, std::vector<std::vector<int>>, int)>& f){

    //Generated data holds all auxiliary DB, such as field and orbits
    GeneratedData params(field, reverse_field, file_generated_matrix);

    std::vector<std::thread> threads;
    int num_threads;
    if(std::thread::hardware_concurrency() > 1){
        num_threads = std::thread::hardware_concurrency() - 1;
    }
    else{
        num_threads = 1;
    }

    for (int i = 0; i < num_threads; i++) {
        threads.push_back(std::thread (threadCalcOrbit, params, i, num_threads,std::ref(f)));
        std::cout << "Started thread n'" << i << "\n";
    }

    for (int i = 0; i < num_threads; i++) {
        threads[i].join();
        std::cout << "Joined thread n'" << i << "\n";
    }

    return 0;
}

int part2a_aux(Field& field, Reverse_Field& reverse_field){
    std::vector<std::vector<int>> file_generated_matrix;

    if(part2_checks(field,reverse_field,file_generated_matrix) != 0){
        return 1;
    }


    std::cout << "Running Kramer-Mesner generator:" << std::endl;
    std::function<int(Field&,Reverse_Field&, std::vector<std::vector<int>>, int)> f = kramer_generator_aux;

    part2_threads(field,reverse_field,file_generated_matrix,f);

    std::cout << "Thank you for waiting :) Merging kramer_mesner files to create kramer_mesner.txt. This could take between 1 to 5 minutes." << std::endl;
    merge_files("kramer_mesner");
    std::cout << "Done!" << std::endl;

    return 0;
}

int part2a(Field& field, Reverse_Field& reverse_field){
    std::cout << "You will now generate the Kramer-Mesner matrix \n"
                 "Warning: file kramer_mesner.txt weighs roughly 338.86MB, give or take 2MB due to the random nature of its generation! \n"
                 "Please make sure you have the necessary space on your machine, otherwise the program has undefined behavior! \n"
                 "The generation of the file takes between 30 seconds to 15 minutes, just so you know :) \n"
                 "If for any reason you wish to quit now input 'q',if you wish to start input 's':" << std::endl;

    std::string prompt;
    std::cin >> prompt;

    while(prompt != "s" && prompt != "q"){
        std::cout << "Please either input 's' to start or input 'q' to quit:" << std::endl;
        std::cin >> prompt;
    }
    if(prompt == "q"){
        std::cout << "Quitting, have a nice day!" << std::endl;
        return 0;
    }

    return part2a_aux(field,reverse_field);
    
}


int part2b(Field& field, Reverse_Field& reverse_field){
    std::cout << "You will now generate the SS3 orbits\n"
                 "Warning: file SS3_orbits.txt weighs roughly 28.06GB! \n"
                 "Please make sure you have the necessary space on your machine, otherwise the program has undefined behavior! \n"
                 "The generation of the file takes between 20 to 24 HOURS!!!, just so you know :) \n"
                 "Note: this function isn't meant to be run more than once, we all have to get out of the house every once in a while :)"
                 "If for any reason you wish to quit now input 'q',if you wish to start input 's':" << std::endl;

    std::string prompt;
    std::cin >> prompt;

    while(prompt != "s" && prompt != "q"){
        std::cout << "Please either input 's' to start or input 'q' to quit:" << std::endl;
        std::cin >> prompt;
    }
    if(prompt == "q"){
        std::cout << "Quitting, have a nice day!" << std::endl;
        return 0;
    }

    //start...

    std::vector<std::vector<int>> file_generated_matrix;

    if(part2_checks(field,reverse_field,file_generated_matrix) != 0){
        return 1;
    }


    std::cout << "Running SS3 generator:" << std::endl;
    std::function<int(Field&,Reverse_Field&, std::vector<std::vector<int>>, int)> f = SS3_generator_aux;

    part2_threads(field,reverse_field,file_generated_matrix,f);



    std::cout << "Thank you for waiting :) Merging SS3_orbits files to create SS3_orbits.txt. This could take between 20 and 35 minutes :) Do not attempt to open the files while this is running: this may result in the files not being able to be deleted and the program crashing. Not the thing you would want after 24 nerve-wrecking hours, innit?" << std::endl;
    merge_files("SS3_orbits");
    std::cout << "Done!" << std::endl;

    return 0;
}

int part2(Field& field, Reverse_Field& reverse_field){

    std::cout << "Action 2 has two independent sub-actions where only sub-action a is needed to later run action 3. \n"
                 "Please choose either 'a' for Kramer-Mesner matrix generation or 'b' for SS3 orbits generation or 'q' to quit" << std::endl;
    std::string abprompt;
    std::cin >> abprompt;

    while(abprompt != "a" && abprompt != "b" && abprompt != "q"){
        std::cout << "Please either input 'a' for Kramer-Mesner generation,'b' for SS3 orbits generation or input 'q' to quit:" << std::endl;
        std::cin >> abprompt;
    }
    if(abprompt == "q"){
        std::cout << "Quitting, have a nice day!" << std::endl;
        return 0;
    }
    else if(abprompt == "a"){
        return part2a(field,reverse_field);
    }

    //part 2b...
    return  part2b(field,reverse_field);

}

int part3(){
    std::cout << "Generating the table from file kramer_mesner.txt. This takes a few seconds..." << std::endl;
    dancing_links DLX("kramer_mesner_reduced");
    std::cout << "Done!" << std::endl;


    node* dummy = DLX.get_dummy();

    std::cout << "Running the Dancing Links solver. Generates solutions in all_solutions.txt" << std::endl;
    DLX.solve();

    std::cout << "Done!" << std::endl;

    return 0;
}

/// runs several iterations, in each: generates a random kramer file and tries to reduce it a selected number of times,
/// randomly, and get a partial solution, it then finds the largest partial solution and deletes all other files.
/// \param field
/// \param reverse_field
/// \return 0 if successful, 1 o.w.
int
partial_solution(Field &field, Reverse_Field &reverse_field, bool memory_based,bool use_existing_kramer_file, int kramer_generation_iterations_num,
                 int kramer_reduction_iterations_num, int max_sol_try, int &max_sol_size, int &max_sol_index_i,
                 int &max_sol_index_j) {



    std::vector<std::vector<int>> file_generated_matrix;

    //check that the SS2 and field files are correct and loads the data unto field, reverse_field
    // and file_generated_matrix
    if(part2_checks(field,reverse_field,file_generated_matrix) != 0){
        return 1;
    }

    std::function<int(Field&,Reverse_Field&, std::vector<std::vector<int>>, int)> f = kramer_generator_aux;

    std::cout << "Running Kramer-Mesner ultimate random partial solution generator:" << std::endl;


    //you get to choose your own number of iterations: the number of times to generate a kramer file and the number
    // of times to try to reduce it into a partial solution

    for(int i = 0; i < kramer_generation_iterations_num; i++){

        std::cout << "Running iteration n' " << i+1 << " out of " << kramer_generation_iterations_num << std::endl;

        matrix matrix_to_reduce;

        if(!memory_based) {

            if(!use_existing_kramer_file){
                std::cout << "Generating random kramer file..." << std::endl;

                //generate the kramer file
                part2_threads(field, reverse_field, file_generated_matrix, f);
                merge_files("kramer_mesner");
                std::cout << "File created" << std::endl;
            }
            


            //create a dancing links kramer object from the file
            //dancing_links matrix_to_reduce("kramer_mesner");
            matrix_to_reduce = matrix("kramer_mesner");
        }
        else{
            for(int orb = 0; orb < SS2_HEIGHT; orb++){
                matrix_generator_aux(field,reverse_field,file_generated_matrix, matrix_to_reduce, orb);
            }
        }

        //save the max solution size amongst all iterations INSIDE the reduce function and the index of the iteration
        // where the solution was found
        int curr_iter_max_sol_size;
        int curr_iter_max_sol_index_j;

        //reduce it said number of times and save the results onto the parameters. Tries to reach the max solution
        // number given.

        matrix_to_reduce.reduce(curr_iter_max_sol_size, curr_iter_max_sol_index_j, kramer_reduction_iterations_num, max_sol_try);

        //found better partial solution
        if(curr_iter_max_sol_size > max_sol_size){
            max_sol_index_i = i;
            max_sol_index_j = curr_iter_max_sol_index_j;
            max_sol_size = curr_iter_max_sol_size;
        }


        //this loop renames the files generated by the reduction file, from having index j to having index i_j
        for(int j = 0; j < kramer_reduction_iterations_num; j++){

            std::cout << "Running renaming n' " << i << "," << j << std::endl;


            std::string kramer_partial_before, kramer_partial_after, kramer_reduced_before, kramer_reduced_after;

            kramer_partial_before = kramer_partial_after = kramer_reduced_before = kramer_reduced_after = "kramer_mesner";

            kramer_partial_after += "_";
            kramer_partial_after += std::to_string(i);
            kramer_reduced_after += "_";
            kramer_reduced_after += std::to_string(i);

            kramer_partial_before += "_";
            kramer_partial_before += std::to_string(j);
            kramer_partial_after += "_";
            kramer_partial_after += std::to_string(j);
            kramer_reduced_before += "_";
            kramer_reduced_before += std::to_string(j);
            kramer_reduced_after += "_";
            kramer_reduced_after += std::to_string(j);

            kramer_partial_before += "_partial_solution.txt";
            kramer_partial_after += "_partial_solution.txt";

            kramer_reduced_before += "_reduced.txt";
            kramer_reduced_after += "_reduced.txt";

            while(rename(kramer_partial_before.c_str(), kramer_partial_after.c_str()) != 0){
                std::cout << "trying to rename partial solution " << i << "," << j << std::endl;
            }

            while(rename(kramer_reduced_before.c_str(), kramer_reduced_after.c_str()) != 0){
                std::cout << "trying to rename reduced kramer " << i << "," << j << std::endl;
            }


        }

    }


    //this double loop deletes all files except for the one with index i_j which is the index of the max partial sol.
    for(int i = 0; i < kramer_generation_iterations_num; i++){
        for(int j = 0; j < kramer_reduction_iterations_num; j++){
            if(!((i == max_sol_index_i) && (j == max_sol_index_j))){
                std::string kramer_random_filename = "kramer_mesner";
                kramer_random_filename += "_";
                kramer_random_filename += std::to_string(i);
                kramer_random_filename += "_";
                kramer_random_filename += std::to_string(j);
                std::string kramer_partial_filename = kramer_random_filename;
                std::string kramer_reduction_filename = kramer_random_filename;
                kramer_partial_filename += "_partial_solution.txt";
                kramer_reduction_filename += "_reduced.txt";

                while(remove(kramer_partial_filename.c_str()) != 0){
                    std::cout << "Removing partial solution file" << i << "," << j << std::endl;
                }

                while(remove(kramer_reduction_filename.c_str()) != 0){
                    std::cout << "Removing reduction file" << i << "," << j << std::endl;

                }
            }
        }
    }

    std::cout << "The max partial solution size is " << max_sol_size << " in partial solution file with indices " << max_sol_index_i
              << "," << max_sol_index_j << std::endl;




    return 0;
}



//this function sucks, 'nuff said...
int partial_increase(){


    dancing_links k_551("kramer_mesner_551_partial_solution");
    dancing_links k_547("kramer_mesner_547_partial_solution");
    dancing_links k_544("kramer_mesner_544_partial_solution");
    dancing_links k_538("kramer_mesner_538_partial_solution");
    dancing_links k_522("kramer_mesner_522_partial_solution");
    dancing_links k_507("kramer_mesner_507_partial_solution");
    dancing_links k_498("kramer_mesner_498_partial_solution");
    dancing_links k_491("kramer_mesner_491_partial_solution");
    dancing_links k_472("kramer_mesner_472_partial_solution");
    dancing_links k_461("kramer_mesner_461_partial_solution");
    dancing_links k_436("kramer_mesner_436_partial_solution");
    dancing_links k_419("kramer_mesner_419_partial_solution");


    k_551.portmanteau(k_547);
    k_551.portmanteau(k_544);
    k_551.portmanteau(k_538);
    k_551.portmanteau(k_522);
    k_551.portmanteau(k_507);
    k_551.portmanteau(k_498);
    k_551.portmanteau(k_491);
    k_551.portmanteau(k_472);
    k_551.portmanteau(k_461);
    k_551.portmanteau(k_436);
    k_551.portmanteau(k_419);


//    std::string file_a;
//    std::string file_b;
//
//    std::cin >> file_a;
//    std::cin >> file_b;
//
//    dancing_links bigger_kramer(file_a);
//    dancing_links smaller_kramer(file_b);
//
//    bigger_kramer.portmanteau(smaller_kramer);




    return 0;


}

//runs several iterations, for each: takes a set number of kramer rows from each orbit
// and creates a kramer file and tries to solve it with dancing links
bool overly_ambitious(Field& field, Reverse_Field& reverse_field){
    for (int i = 0; i < 1; ++i)
    {
        std::cout << "Running iteration n' " << (i+1) << std::endl;
        part2a_aux(field,reverse_field);
        dancing_links DLX("kramer_mesner");
        DLX.solve();
        if(DLX.get_num_of_solutions() > 0){
            return true;
        }
    }
    return false;
}

int build_reduced_kramer_from_partial_solution(Field& field, Reverse_Field& reverse_field, std::string& filename){

    std::vector<std::vector<int>> file_generated_matrix;

    //check that the SS2 and field files are correct and loads the data unto field, reverse_field
    // and file_generated_matrix
    if(part2_checks(field,reverse_field,file_generated_matrix) != 0){
        return 1;
    }

    //std::function<int(Field&,Reverse_Field&, std::vector<std::vector<int>>,char*, int)> g = reduced_kramer_generator_aux;

    char* cols = new char[SS2_HEIGHT];

    dancing_links partial_sol(filename);
    partial_sol.cols_arr(cols);


//Generated data holds all auxiliary DB, such as field and orbits
    GeneratedData params(field, reverse_field, file_generated_matrix,cols);

    std::vector<std::thread> threads;
    int num_threads;
    if(std::thread::hardware_concurrency() > 1){
        num_threads = std::thread::hardware_concurrency() - 1;
    }
    else{
        num_threads = 1;
    }

    for (int i = 0; i < num_threads; i++) {
        threads.push_back(std::thread (threadCalcReduced, params, i, num_threads));
        std::cout << "Started thread n'" << i << "\n";
    }

    for (int i = 0; i < num_threads; i++) {
        threads[i].join();
        std::cout << "Joined thread n'" << i << "\n";
    }

    merge_files("kramer_mesner",cols);

    delete[] cols;
    return 0;
}



int reduced_file_to_original(std::string filename_arg){
    std::string filename = filename_arg;
    filename += ".txt";

    std::ifstream myfile(filename);

    if(!myfile.is_open()){
        return 1;
    }


    std::ofstream original("orig.txt");

    if(!original.is_open()){
        return 1;
    }

    std::string line;
    std::string elem;



    char arr[SS2_HEIGHT];

    for(int i = 0; i < SS2_HEIGHT; i++){
        arr[i] = 0;
    }



    getline(myfile,line);

    std::istringstream stream(line);


    getline(stream,elem,',');

    while(getline(stream,elem,',')){
        int index = std::stoi(elem);
        arr[index] = 1;
    }




    original << " ,";
    for(int i = 0; i < SS2_HEIGHT; i++){
        original << std::to_string(i) << ",";
    }
    original << std::endl;

    while(getline(myfile,line)){
        std::istringstream stream(line);
        getline(stream,elem,',');
        original << elem << ",";

        for(int i = 0; i < SS2_HEIGHT; i++){

            if(arr[i] == 1){
                getline(stream,elem,',');
                original << elem << ",";
            }
            else{
                original << "0,";
            }


        }
        original << std::endl;

    }

    while(rename("orig.txt", filename.c_str()) != 0){
        std::cout << "trying to rename the file" << std::endl;
    }

    return 0;
}

int reduce_and_deduce(Field &field, Reverse_Field &reverse_field, bool start_from_partial_file) {

    int start = 500;
    int best_ever_sol_size = 593;
    //300 rows from each orbit, specified in matrix_generator_aux
    std::string curr_partial_filename = "kramer_mesner_593_partial_solution";


    for(int iter = 0; iter < 1 && best_ever_sol_size < 657; iter++) {

        if(!start_from_partial_file){
            int max_sol_size = 0;
            int max_sol_index_i = 0;
            int max_sol_index_j = 0;
            partial_solution(field, reverse_field, true, false, 1, 3, start, max_sol_size, max_sol_index_i,
                             max_sol_index_j);

            std::string r_filename = "kramer_mesner_";
            r_filename += std::to_string(max_sol_index_i);
            r_filename += "_";
            r_filename += std::to_string(max_sol_index_j);
            r_filename += "_reduced.txt";

            while (remove(r_filename.c_str()) != 0) {
                std::cout << "trying to remove the reduced file" << r_filename << std::endl;
            }

            std::string p_filename = "kramer_mesner_";
            p_filename += std::to_string(max_sol_index_i);
            p_filename += "_";
            p_filename += std::to_string(max_sol_index_j);
            p_filename += "_partial_solution.txt";

            std::string n_filename = "kramer_mesner_";
            n_filename += std::to_string(max_sol_size);
            n_filename += "_partial_solution.txt";

            while (rename(p_filename.c_str(), n_filename.c_str()) != 0) {
                std::cout << "trying to rename the max file" << p_filename << std::endl;
            }

            std::string par_sol_arg_filename = n_filename.substr(0, n_filename.length() - 4);

            build_reduced_kramer_from_partial_solution(field, reverse_field, par_sol_arg_filename);

            dancing_links reduced_kramer("kramer_mesner");
            reduced_kramer.solve();
            if (reduced_kramer.get_num_of_solutions() > 0) {
                std::cout << "Eureka!" << std::endl;
                return 0;
            }

            int max_red_sol_size = 0;
            int max_red_sol_index_i = 0;
            int max_red_sol_index_j = 0;

            partial_solution(field, reverse_field, false, true, 1, 10, (657 - max_sol_size), max_red_sol_size,
                             max_red_sol_index_i, max_red_sol_index_j);

            std::string r_r_filename = "kramer_mesner_";
            r_r_filename += std::to_string(max_red_sol_index_i);
            r_r_filename += "_";
            r_r_filename += std::to_string(max_red_sol_index_j);
            r_r_filename += "_reduced.txt";

            while (remove(r_r_filename.c_str()) != 0) {
                std::cout << "trying to remove the reduced file " << r_r_filename << std::endl;
            }

            std::string p_r_filename = "kramer_mesner_";
            p_r_filename += std::to_string(max_red_sol_index_i);
            p_r_filename += "_";
            p_r_filename += std::to_string(max_red_sol_index_j);
            p_r_filename += "_partial_solution.txt";

            std::string red_par_sol_arg_filename = p_r_filename.substr(0, p_r_filename.length() - 4);

            if (max_sol_size + max_red_sol_size > best_ever_sol_size) {
                best_ever_sol_size = max_sol_size + max_red_sol_size;

                reduced_file_to_original(red_par_sol_arg_filename);

                dancing_links smaller_kramer(red_par_sol_arg_filename);
                dancing_links bigger_kramer(par_sol_arg_filename);
                bigger_kramer.portmanteau(smaller_kramer);

                std::string new_best_filename = "kramer_mesner_";
                new_best_filename += std::to_string(best_ever_sol_size);
                new_best_filename += "_partial_solution.txt";

                while (rename("portmanteau.txt", new_best_filename.c_str()) != 0) {
                    std::cout << "trying to rename the new best file" << new_best_filename << std::endl;
                }
            }


            while (remove(p_r_filename.c_str()) != 0) {
                std::cout << "trying to remove the reduced file" << p_r_filename << std::endl;
            }

            while (remove(n_filename.c_str()) != 0) {
                std::cout << "trying to remove the reduced file" << n_filename << std::endl;
            }
        }
        else{
            build_reduced_kramer_from_partial_solution(field, reverse_field, curr_partial_filename);

            dancing_links reduced_kramer("kramer_mesner");
            reduced_kramer.solve();
            if (reduced_kramer.get_num_of_solutions() > 0) {
                std::cout << "Eureka!" << std::endl;
                return 0;
            }

            int max_red_sol_size = 0;
            int max_red_sol_index_i = 0;
            int max_red_sol_index_j = 0;

            partial_solution(field, reverse_field, false, true, 1, 10, (657 - best_ever_sol_size), max_red_sol_size,
                             max_red_sol_index_i, max_red_sol_index_j);

            std::string r_r_filename = "kramer_mesner_";
            r_r_filename += std::to_string(max_red_sol_index_i);
            r_r_filename += "_";
            r_r_filename += std::to_string(max_red_sol_index_j);
            r_r_filename += "_reduced.txt";

            while (remove(r_r_filename.c_str()) != 0) {
                std::cout << "trying to remove the reduced file " << r_r_filename << std::endl;
            }

            std::string p_r_filename = "kramer_mesner_";
            p_r_filename += std::to_string(max_red_sol_index_i);
            p_r_filename += "_";
            p_r_filename += std::to_string(max_red_sol_index_j);
            p_r_filename += "_partial_solution.txt";

            std::string red_par_sol_arg_filename = p_r_filename.substr(0, p_r_filename.length() - 4);

            if (best_ever_sol_size + max_red_sol_size > best_ever_sol_size) {

                best_ever_sol_size = best_ever_sol_size + max_red_sol_size;

                reduced_file_to_original(red_par_sol_arg_filename);

                dancing_links smaller_kramer(red_par_sol_arg_filename);
                dancing_links bigger_kramer(curr_partial_filename);
                bigger_kramer.portmanteau(smaller_kramer);

                std::string new_best_filename = "kramer_mesner_";
                new_best_filename += std::to_string(best_ever_sol_size);
                new_best_filename += "_partial_solution.txt";


                while (rename("portmanteau.txt", new_best_filename.c_str()) != 0) {
                    std::cout << "trying to rename the new best file" << new_best_filename << std::endl;
                }

                curr_partial_filename = new_best_filename.substr(0, new_best_filename.length() - 4);
            }

            while (remove(p_r_filename.c_str()) != 0) {
                std::cout << "trying to remove the reduced file" << p_r_filename << std::endl;
            }

        }


    }


    return 0;
}

int main(int argc, const char * argv[]) {

    Field field;
    Reverse_Field reverse_field;




    std::cout << "Hello and welcome to the q-analog solver for GF(2,19) \n"
            "Please select one of the following actions (1,2,3) or input 'q' to quit: \n"
            "1. Generate GF(2,19) field and the SS2 orbits. \n"
            "Output: \n"
            "field2.csv - a file containing all members of GF(2,19) with the pattern: degree,bit representation, orbit representative degree. \n"
            "SS2_orbits.txt - a file containing all SS2 orbits with each row i containing all degrees of members of the i'th orbit starting with the orbit representative. \n"
            "2. Generate SS3_orbits and Kramer-Mesner matrix used for the solving algorithm. \n"
            "Output: \n"
            "SS3_orbits.txt - a file containing all orbit representatives for the SS3 orbits with pattern: degree1,degree2 \n"
            "kramer_mesner.txt - a file of the Kramer-Mesner matrix, taking 8 rows from each of the 4599 orbits, each row has 7 1's in it\n"
            "3. Solve exact cover problem visualized by the Kramer-Mesner matrix using Knuth's DLX. \n"
            "Output: \n"
            "" << std::endl;

    std::string prompt;
    std::cin >> prompt;

    while(prompt != "1" && prompt != "2" && prompt != "3" && prompt != "q"){
        std::cout << "Please select one of the above actions (1,2,3) or input 'q' to quit:" << std::endl;
        std::cin >> prompt;
    }

    if(prompt == "q"){
        std::cout << "Quitting, have a nice day!" << std::endl;
        return 0;
    }
    else if(prompt == "1"){
        return part1(field,reverse_field);
    }
    else if(prompt == "2"){
        int res = part2(field,reverse_field);
        if(res != 0){
            std::cout << "Uh-oh! Looks like you've encountered a problem. \n"
                    " If there is a problem with one of the files either run the program again, this time SELECT 1 to generate the files correctly \n "
                    "or if you've received the files from a source that is not this program either check the authenticity of the source and if you trust it download again - if the problem persists run this program again and SELECT 1" << std::endl;
        }
        return res;
    }

    //part3...
    return part3();

}







