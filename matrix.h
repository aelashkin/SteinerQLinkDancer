//
// Created by andy on 8/23/19.
//

#ifndef MATRIX_BUILDER_MATRIX_H
#define MATRIX_BUILDER_MATRIX_H

#include "constants.h"
#include "dancing_links.h"


class row{
public:

    int arr[7];
    bool in_solution = false;
    std::string indices;

    row() = default;

    row(const row& r) {
        in_solution = r.in_solution;
        indices = r.indices;
        for(int i = 0; i < 7; i++){
            arr[i] = r.arr[i];
        }
    }
};

class matrix {

private:
    std::unordered_map<int,row> rows_arr;
    std::string problem_file;
    int n_cols = 0;
    int row_count = 0;
    bool memory_based = false;
    std::string first_row = "";

public:

    //copies the rows of the matrix object to the vector param.
    void rows_map_to_vector(std::vector<row>& v){
        for(int i = 0; i < row_count; i++){
            v.insert(v.begin() + i, rows_arr[i]);
        }


    }

    // Tries to find the best partial solution in the matrix object, greedily.
    // recieves vectors meant to hold the current partial sol. and best partial sol.
    // and an array of the currently covered cols
    void best_partial_solution(int k, int& best_par_sol_num,char*& cols,std::vector<int>& curr_par_sol,std::vector<int>& best_par_sol){

        std::cout << k << std::endl;

        for(int i = 0; i < row_count; i++){

            bool already_covered = false;
            for(int j = 0; j < 7; j++){
                if(cols[rows_arr[i].arr[j]] == 1){
                    already_covered = true;
                    break;
                }
            }
            if(already_covered){
                continue;
            }


            for(int j = 0; j < 7; j++){
                cols[rows_arr[i].arr[j]] = 1;
            }

            curr_par_sol.push_back(i);

            if(k + 1 > best_par_sol_num){
                best_par_sol_num = k + 1;
                best_par_sol = curr_par_sol;
            }

            best_partial_solution(k + 1,best_par_sol_num,cols,curr_par_sol,best_par_sol);

            curr_par_sol.pop_back();

            for(int j = 0; j < 7; j++){
                cols[rows_arr[i].arr[j]] = 0;
            }



        }


        std::cout << k << std::endl;

    }

    //tries to attach a number of rows to the current partial solution
    void attach_aux(int l, int& best_par_sol_num, char*& cols,std::vector<row>& par_sol,std::vector<row>& best_par_sol,int how_many_rows_to_attach, std::vector<row>& smaller_par_sol,bool& finished){

        std::cout << "Attach: " << l << std::endl;


        if(l >= how_many_rows_to_attach){

            if(par_sol.size() > best_par_sol_num) {
                best_par_sol_num = par_sol.size();
                best_par_sol = std::vector<row>(par_sol);

                finished = true;
            }

            return;
        }

        for(int i = 0; i < smaller_par_sol.size(); i++){

            row curr_row = smaller_par_sol[i];

            bool already_covered = false;
            for(int j = 0 ; j < 7; j++){
                if(cols[smaller_par_sol[i].arr[j]] == 1){
                    already_covered = true;
                    break;
                }
            }
            if(already_covered){
                continue;
            }

            for(int j = 0 ; j < 7; j++){
                cols[smaller_par_sol[i].arr[j]] = 1;
            }

            par_sol.push_back(curr_row);

            attach_aux(l + 1, best_par_sol_num,cols,par_sol,best_par_sol,how_many_rows_to_attach,smaller_par_sol,finished);

            par_sol.pop_back();

            for(int j = 0 ; j < 7; j++){
                cols[smaller_par_sol[i].arr[j]] = 0;
            }

            if(finished){
                break;
            }
        }

        std::cout << "Attach: " << l << std::endl;

    }

    //detaches some rows from the current partial sol. and then tries to attach some more rows from a vector of other rows
    void detach_some_and_attach_some_more(int k,int& best_par_sol_num,char*& cols,std::vector<row>& par_sol,std::vector<row>& best_par_sol, int how_many_rows_to_detach, int how_many_rows_to_attach, std::vector<row>& smaller_par_sol,bool& finished){
        std::cout << "Detach: " << k << std::endl;

        if(k >= how_many_rows_to_detach){

            attach_aux(0, best_par_sol_num,cols,par_sol,best_par_sol,how_many_rows_to_attach,smaller_par_sol,finished);


            return;
        }


        for(int i = 0; i < par_sol.size(); i++){

            row curr_row = par_sol[i];

            par_sol.erase(par_sol.begin() + i);

            for(int j = 0 ; j < 7; j++) {
                cols[curr_row.arr[j]] = 0;
            }

            detach_some_and_attach_some_more(k + 1, best_par_sol_num,cols,par_sol,best_par_sol,how_many_rows_to_detach,how_many_rows_to_attach,smaller_par_sol, finished);

            for(int j = 0 ; j < 7; j++) {
                cols[curr_row.arr[j]] = 1;
            }

            par_sol.insert(par_sol.begin() + i, curr_row);


            if(finished){
                break;
            }
        }


        std::cout << "Detach: " << k << std::endl;

    }

    //tries to remove a set number of rows from the partial solution and add one than that to it to increase the par. sol. by 1.
    static void detach_some_and_attach_some_more(int how_many_rows_to_detach,int& best_par_sol_num,std::vector<row>& best_par_sol){
        dancing_links a("kramer_mesner");
        char* cols = new char[SS2_HEIGHT];
        a.cols_arr(cols);
        std::vector<row> kramer_rows;
        matrix bigger_kramer("kramer_mesner");
        bigger_kramer.rows_map_to_vector(kramer_rows);

        std::vector<row> small_rows;
        matrix smaller_kramer("small_sol");
        smaller_kramer.rows_map_to_vector(small_rows);

        bool finished = false;
        bigger_kramer.detach_some_and_attach_some_more(0, best_par_sol_num,cols, kramer_rows,best_par_sol, how_many_rows_to_detach, how_many_rows_to_detach + 1, small_rows,finished);
        delete[] cols;
    }

    /// reduces the kramer matrix object to a partial solution and outputs two sets of files: the random partial
    /// solution of each iterations is one type of file and the other is the kramer matrix file after being reduced by
    /// complying with the columns of the partial solution
    /// \param max_partial_solution returns the number of rows of the chosen partial solution, the one with
    /// the largest amount of rows
    /// \param max_index returns the index of said chosen partial solution, amongst all iterations
    /// \param num_of_iterations number of random solutions to try and generate
    /// \param partial_solution_size the number of rows in the wanted partial solution
    /// \return 0 if successful, 1 o.w.
    int reduce(int& max_partial_solution, int& max_index, int num_of_iterations = 1,int partial_solution_size = 657){

        max_partial_solution = 0;
        max_index = 0;


        for(int iter = 0; iter < num_of_iterations; iter++) {

            std::cout << "Partial solution iteration n' " << (iter + 1) << " of " << num_of_iterations << std::endl;


            //to mark columns in the solution
            char bit_array[n_cols];

            for (int i = 0; i < n_cols; i++) {
                bit_array[i] = 0;
            }

            int rows_taken = 0;


            int num_iterations = 0;

            //we would like to take as many as many rows as we've prescribed in the parameter
            while (rows_taken < partial_solution_size) {

                //a break condition, empirically selected :)
                num_iterations++;
                if (num_iterations >= 10 * row_count) {
                    std::cout << "Couldn't find " << partial_solution_size << " rows this time , only " << rows_taken
                              << " :)" << std::endl;
                    break;
                }

                //generate random number
                std::random_device dev;
                MyRNG rng(dev());
                std::uniform_int_distribution<MyRNG::result_type> dist(0, row_count - 1);
                int j = dist(rng);

                row row = rows_arr[j];

                //if it is already in the solution
                if (row.in_solution) {
                    continue;
                }

                //check that the selected row isn't already covered by our current partial solution
                bool already_covered = false;
                for(int i = 0; i < 7; i++){
                    if(bit_array[row.arr[i]] == 1){
                        already_covered = true;
                        break;
                    }
                }
                if(already_covered){
                    continue;
                }


                //it is not covered, so we will add it to our partial solution and mark its columns
                for(int i = 0; i < 7; i++){
                    bit_array[row.arr[i]] = 1;
                }

                rows_arr[j].in_solution = true;

                rows_taken++;

                std::cout << "Taken row n' " << j << " to the partial solution. |Solution| = " << rows_taken
                          << std::endl;

            }

            //found better partial solution
            if(rows_taken > max_partial_solution){
                max_partial_solution = rows_taken;
                max_index = iter;
            }



            std::string reduced_filename;
            std::string partial_filename;
            std::string line;
            std::string elem;


            //open (possibly) three files: the original for reading and the partial solution and reduced files for writing
            reduced_filename = problem_file;
            reduced_filename += "_";
            reduced_filename += std::to_string(iter);
            reduced_filename += "_reduced.txt";
            std::ofstream kramer_reduced(reduced_filename);

            if (!kramer_reduced.is_open()) {
                return 1;
            }

            partial_filename = problem_file;
            partial_filename += "_";
            partial_filename += std::to_string(iter);
            partial_filename += "_partial_solution.txt";
            std::ofstream kramer_partial(partial_filename);

            if (!kramer_partial.is_open()) {
                return 1;
            }






            if(!memory_based) {
                std::string original_filename;

                original_filename = problem_file;
                original_filename += ".txt";
                std::ifstream kramer_original(original_filename);

                if (!kramer_original.is_open()) {
                    return 1;
                }

                ///


                getline(kramer_original, line);
                std::istringstream stream(line);
                getline(stream, elem, ',');

                //first row...
                kramer_reduced << " ,";
                kramer_partial << " ,";
                for (int i = 0; i < n_cols; i++) {
                    getline(stream,elem,',');

                    if (bit_array[i] == 0) {
                        kramer_reduced << elem << ",";
                    }

                    kramer_partial << elem << ",";
                }

                kramer_reduced << std::endl;

                kramer_partial << std::endl;




                //write info to files...
                for (int i = 0; i < row_count; i++) {

                    row row = rows_arr[i];

                    getline(kramer_original, line);


                    //if the row is not in the partial solution, we must first check that it complies with it
                    //before inserting it to the reduced file.
                    if (!row.in_solution) {
                        bool already_covered = false;
                        for(int j = 0; j < 7; j++){
                            if(bit_array[row.arr[j]] == 1){
                                already_covered = true;
                            }
                        }
                        if(already_covered){
                            continue;
                        }

                    }


                    //now, the row is either in the partial solution or complies with its columns, or both
                    //insert first element of row: the i.j indices of the orbit

                    std::istringstream stream(line);

                    getline(stream, elem, ',');
                    kramer_reduced << elem << ",";



                    if (row.in_solution) {


                        kramer_partial << elem << ",";

                    }

                    //add the other parts of the row. The reduced file contains only the uncovered columns
                    for (int j = 0; j < n_cols; j++) {

                        getline(stream, elem, ',');



                        if (bit_array[j] == 0) {
                            kramer_reduced << elem << ",";
                        }

                        if (row.in_solution) {
                            kramer_partial << elem << ",";
                        }


                    }

                    //add line break
                    kramer_reduced << std::endl;

                    if (row.in_solution) {
                        kramer_partial << std::endl;
                    }

                }


                ///
                kramer_original.close();
            }
            else{
                //@

                //first row...
                kramer_reduced << " ,";
                kramer_partial << " ,";

                if(!first_row.empty()) {
                    std::stringstream stream(first_row);
                    getline(stream, elem, ',');
                    for (int i = 0; i < n_cols; i++) {

                        getline(stream,elem,',');
                        if (bit_array[i] == 0) {
                            kramer_reduced << elem << ",";
                        }

                        kramer_partial << elem << ",";
                    }
                }
                else{
                    for(int i = 0; i < n_cols; i++){
                        elem = std::to_string(i);
                        if(bit_array[i] == 0){
                            kramer_reduced << elem << ",";
                        }

                        kramer_partial << elem << ",";

                    }
                }


                kramer_reduced << std::endl;

                kramer_partial << std::endl;




                //write info to files...
                for (int i = 0; i < row_count; i++) {

                    row row = rows_arr[i];

                    //if the row is not in the partial solution, we must first check that it complies with it
                    //before inserting it to the reduced file.
                    if (!row.in_solution) {
                        bool already_covered = false;
                        for(int j = 0; j < 7; j++){
                            if(bit_array[row.arr[j]] == 1){
                                already_covered = true;
                            }
                        }
                        if(already_covered){
                            continue;
                        }

                    }


                    //now, the row is either in the partial solution or complies with its columns, or both
                    //insert first element of row: the i.j indices of the orbit

                    kramer_reduced << rows_arr[i].indices << ",";


                    if (row.in_solution) {

                        kramer_partial << rows_arr[i].indices << ",";

                    }

                    //add the other parts of the row. The reduced file contains only the uncovered columns
                    for (int j = 0; j < n_cols; j++) {


                        elem = "0";

                        for(int k = 0; k < 7; k++){
                            if(rows_arr[i].arr[k] == j){
                                elem = "1";
                            }
                        }



                        if (bit_array[j] == 0) {
                            kramer_reduced << elem << ",";
                        }

                        if (row.in_solution) {
                            kramer_partial << elem << ",";
                        }


                    }

                    //add line break
                    kramer_reduced << std::endl;

                    if (row.in_solution) {
                        kramer_partial << std::endl;
                    }

                }
                //@
            }





            //fix rows_arr
            for (int i = 0; i < row_count; i++) {
                rows_arr[i].in_solution = false;
            }


            //close files
            kramer_partial.close();
            kramer_reduced.close();

        }

        return 0;
    }

    matrix(std::string file = ""){

        if(!file.empty()) {
        	problem_file = file;

            file += ".txt";


            try {
                make(file);
            }
            catch (std::bad_alloc &e) {
                destroy();
            }
        }
        else{
            memory_based = true;
            n_cols = SS2_HEIGHT;
            problem_file = "kramer_mesner";
        }
    }

    void add_row(int arr[7], std::string str){
        if(memory_based) {
            for (int i = 0; i < 7; i++) {
                rows_arr[row_count].arr[i] = arr[i];
            }
            rows_arr[row_count].indices = str;
            row_count++;
        }
    }

    void add_first_row(std::string f){
        first_row = f;
    }

    void make(std::string filename){


        std::ifstream myfile (filename);
        std::string line;

        getline(myfile,line);

        std::istringstream first_stream(line);
        std::string col;



        getline(first_stream,col,',');

        int n = 0;

        //build headers
        while( getline(first_stream, col, ',') ){
            if(!col.empty()){
                n++;
            }
        }

        n_cols = n;



        int row_index = 0;
        while(getline(myfile,line)){

            std::istringstream stream(line);

            std::string elem;

            getline(stream, elem, ',');

            rows_arr[row_index].indices = elem;


            int col_index = 0;
            int place_in_row = 0;

            while( getline(stream, elem, ',') )
            {

                if(!elem.empty()){


                    int bit = std::stoi(elem);

                    if(bit == 1){

                        rows_arr[row_index].arr[place_in_row] = col_index;
                        place_in_row++;



                    }

                }


                col_index++;

            }


            row_index++;
        }

        myfile.close();

        row_count = row_index;
    }

    void destroy(){

    }

    ~matrix(){
        destroy();
    }



};


#endif //MATRIX_BUILDER_MATRIX_H
