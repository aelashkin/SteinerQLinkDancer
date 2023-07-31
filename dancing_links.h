//
// Created by User on 8/4/2019.
//

#ifndef MATRIX_BUILDER_DANCING_LINKS_H
#define MATRIX_BUILDER_DANCING_LINKS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include "constants.h"
#include "vec_n.h"


class node{
public:
    node* left,*right,*up,*down,*header;
    int row,col, count;
};



class dancing_links {
private:
    node* dummy;
    node** header_arr;
    int n_cols;
    std::vector<int> solutions;
    int num_of_solutions;
    std::unordered_map<int,node*> rows_arr;
    //std::unordered_map<int,std::string> rows_orbit_degrees;//haven't used it yet
    int row_count;
    std::string problem_file;

public:


    //Treats arr is a bit array, sets the bits in arr in the indices where the cols are covered, nulls bits where they are not.
    void cols_arr(char*& arr){
        for(int i = 0; i < n_cols; i++){
            int cnt = header_arr[i]->count;
            if(cnt > 0){
                arr[i] = 1;
            }
            else{
                arr[i] = 0;
            }
        }
    }

    ///
    /// \return if all the rows of the matrix make up a partial solution - each column has at most one 1.
    bool is_matrix_partial_solution(){

        int cnt_array[n_cols];


        for (int i = 0; i < n_cols; i++) {
            cnt_array[i] = 0;
        }

        for(int i = 0; i < row_count; i++){
            node *row = rows_arr[i];
            node* curr = row;

            do {

                cnt_array[curr->header->col]++;
                curr = curr->right;
            } while (curr != row);

        }



        bool res = true;

        for (int i = 0; i < n_cols; i++) {
            if(!((cnt_array[i] == 0) || (cnt_array[i] == 1))){
                res = false;
                std::cout << i << "," << cnt_array[i] << std::endl;
            }
        }



        return res;

    }


    /// A function that is to be invoked by an object signifying a partial solution that tries to add rows from another
    /// partial solution matrix that don't intersect with the object's rows
    /// \param smaller_kramer the object of the other matrix whose rows we are trying to add to our partial solution matrix
    /// \return 0 if successful, 1 o.w.
    int portmanteau(dancing_links& smaller_kramer){


        std::string filename_a = problem_file;
        filename_a += ".txt";
        std::ifstream file_a(filename_a);

        if (!file_a.is_open()) {
            return 1;
        }

        std::string filename_b = smaller_kramer.problem_file;
        filename_b += ".txt";
        std::ifstream file_b(filename_b);

        if (!file_b.is_open()) {
            return 1;
        }

        std::string filename_c = "portmanteau.txt";
        std::ofstream file_c(filename_c);

        if (!file_c.is_open()) {
            return 1;
        }



        std::cout << "The current partial solution is of size: " << row_count << std::endl;



        int bit_array[n_cols];

        for (int i = 0; i < n_cols; i++) {
            bit_array[i] = 0;
        }



        for(int i = 0; i < row_count; i++){
            node *row = rows_arr[i];
            node* curr = row;

            do {

                bit_array[curr->header->col] = 1;
                curr = curr->right;
            } while (curr != row);

        }

        file_c << file_a.rdbuf();

        file_a.close();


        int rows_taken = 0;

        std::string line;

        getline(file_b,line);

        for(int i = 0; i < smaller_kramer.row_count; i++){


            /////check that the row is consistent with our partial solution
            node *row = smaller_kramer.rows_arr[i];

            node *curr = row;
            bool already_covered = false;
            do {
                if (bit_array[curr->header->col] == 1) {
                    already_covered = true;
                }
                curr = curr->right;
            } while (curr != row);


            if (already_covered) {
                getline(file_b,line);
                continue;
            }

            /////add it to the partial solution

            curr = row;

            do{

                bit_array[curr->header->col] = 1;
                curr = curr->right;
            } while(curr != row);

            getline(file_b,line);

            file_c << line << std::endl;
            rows_taken++;


        }

        file_b.close();
        file_c.close();

        std::cout << "The partial solution in file portmanteau.txt is of size: " << (row_count + rows_taken) << std::endl;

        return 0;


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

            //to restore table object info at the end of the iteration
            std::vector<int> row_indices(partial_solution_size);
            std::vector<node *> row_pointers(partial_solution_size);


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

                node *row = rows_arr[j];

                //if it is already in the solution
                if (row == nullptr) {
                    continue;
                }

                //check that the selected row isn't already covered by our current partial solution
                node *curr = row;
                bool already_covered = false;
                do {
                    if (bit_array[curr->header->col] == 1) {
                        already_covered = true;
                    }
                    curr = curr->right;
                } while (curr != row);

                if (already_covered) {
                    continue;
                }


                //it is not covered, so we will add it to our partial solution and mark its columns
                curr = row;
                do {
                    bit_array[curr->header->col] = 1;
                    curr = curr->right;
                } while (curr != row);

                rows_arr[j] = nullptr;


                //save its data (since the row array pointer now points to null) and restore it in the iteration's end
                row_indices[rows_taken] = j;
                row_pointers[rows_taken] = row;

                rows_taken++;

                std::cout << "Taken row n' " << j << " to the partial solution. |Solution| = " << rows_taken
                          << std::endl;

            }

            //found better partial solution
            if(rows_taken > max_partial_solution){
                max_partial_solution = rows_taken;
                max_index = iter;
            }



            //open three files: the original for reading and the partial solution and reduced files for writing
            std::string reduced_filename = problem_file;
            reduced_filename += "_";
            reduced_filename += std::to_string(iter);
            reduced_filename += "_reduced.txt";
            std::ofstream kramer_reduced(reduced_filename);

            if (!kramer_reduced.is_open()) {
                return 1;
            }

            std::string original_filename = problem_file;
            original_filename += ".txt";
            std::ifstream kramer_original(original_filename);

            if (!kramer_original.is_open()) {
                return 1;
            }

            std::string partial_filename = problem_file;
            partial_filename += "_";
            partial_filename += std::to_string(iter);
            partial_filename += "_partial_solution.txt";
            std::ofstream kramer_partial(partial_filename);

            if (!kramer_partial.is_open()) {
                return 1;
            }

            std::string line;

            getline(kramer_original, line);


            std::istringstream stream(line);
            std::string elem;
            getline(stream,elem,',');

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

                node *row = rows_arr[i];

                getline(kramer_original, line);

                //if the row is not in the partial solution, we must first check that it complies with it
                //before inserting it to the reduced file.
                if (row != nullptr) {
                    node *curr = row;
                    bool already_covered = false;
                    do {
                        if (bit_array[curr->header->col] == 1) {
                            already_covered = true;
                        }
                        curr = curr->right;
                    } while (curr != row);

                    if (already_covered) {
                        continue;
                    }


                }


                //now, the row is either in the partial solution or complies with its columns, or both
                //insert first element of row: the i.j indices of the orbit
                std::istringstream stream(line);
                getline(stream, elem, ',');
                kramer_reduced << elem << ",";

                if (row == nullptr) {

                    kramer_partial << elem << ",";
                }

                //add the other parts of the row. The reduced file contains only the uncovered columns
                for (int j = 0; j < n_cols; j++) {
                    getline(stream, elem, ',');
                    if (bit_array[j] == 0) {
                        kramer_reduced << elem << ",";
                    }

                    if (row == nullptr) {
                        kramer_partial << elem << ",";
                    }
                }

                //add line break
                kramer_reduced << std::endl;

                if (row == nullptr) {
                    kramer_partial << std::endl;
                }

            }


            //fix rows_arr
            for (int i = 0; i < rows_taken; i++) {
                rows_arr[row_indices[i]] = row_pointers[i];
            }

            //close files
            kramer_partial.close();
            kramer_reduced.close();
            kramer_original.close();

        }

        return 0;
    }

    int choose_col(){
        int min = dummy->right->count;
        int min_index = dummy->right->col;

        for(node* curr = dummy->right; curr != dummy; curr = curr->right){

            int curr_count = curr->count;
            if(curr_count < min){
                min = curr_count;
                min_index = curr->col;
            }

        }
        return min_index;


    }

    void cover(int i){
        node* c = header_arr[i];
        c->right->left = c->left;
        c->left->right = c->right;
        for(node* row = c->down; row != c; row = row->down){
            for(node* col = row->right; col != row; col = col->right){
                col->down->up = col->up;
                col->up->down = col->down;
                col->header->count--;
            }
        }
    }

    void uncover(int i){
        node* c = header_arr[i];
        for(node* row = c->up; row != c; row = row->up){
            for(node* col = row->left; col != row; col = col->left){
                col->header->count++;
                col->down->up = col;
                col->up->down = col;

            }
        }
        c->right->left = c;
        c->left->right = c;
    }

    int get_num_of_solutions(){
        return num_of_solutions;
    }


    void solve(){
        std::ofstream myfile ("all_solutions.txt");

        if(!myfile.is_open()){
            return;
        }

        num_of_solutions = 0;
        solve(0,myfile);

        myfile.close();

    }


    void solve(int k, std::ofstream& myfile){

        std::cout << k << std::endl;


        if(dummy->right == dummy){


            for(auto iter = solutions.begin(); iter != solutions.end(); iter++){
                myfile << *iter << ",";
            }
            myfile << std::endl;

            num_of_solutions++;
            std::cout << "Found solution n' " << num_of_solutions << " at depth " << k << std::endl;

            return;
        }

        int i = choose_col();
        node* c = header_arr[i];

        cover(i);

        for(node* row = c->down; row != c; row = row->down){

            solutions.push_back(row->row);


            for(node* col = row->right; col != row; col = col->right){

                cover(col->header->col);

            }

            solve(k+1, myfile);

            solutions.pop_back();

            c = row->header;

            for(node* col = row->left; col != row; col = col->left){
                uncover(col->header->col);

            }

        }

        uncover(c->header->col);

    }




    dancing_links(std::string file = "kramer_mesner"){
        problem_file = file;
        file += ".txt";

        try {
            make(file);
        }
        catch(std::bad_alloc& e){
            destroy();
        }
    }

    node*& get_dummy(){
        return dummy;
    }


    void make(std::string filename){


        std::ifstream myfile (filename);
        std::string line;

        getline(myfile,line);

        std::istringstream first_stream(line);
        std::string col;

        dummy = new node();

        int n = 0;

        node* curr = dummy;


        getline(first_stream,col,',');


        //build headers
        while( getline(first_stream, col, ',') ){
            if(!col.empty()){
                node* next = new node();
                curr->right = next;
                next->left = curr;
                next->count = 0;

                //next->col = std::stoi(col);
                next->col = n;

                next->header = next;

                next->up = next;
                next->down = next;

                curr = next;
                n++;
            }
        }

        curr->right = dummy;
        dummy->left = curr;

        header_arr = new node*[n];
        n_cols = n;


        int i = 0;
        for(node* curr = dummy->right; curr != dummy; curr = curr->right){
            header_arr[i] = curr;
            i++;
        }




        int row_index = 0;
        while(getline(myfile,line)){

            std::istringstream stream(line);

            std::string elem;

            getline(stream, elem, ',');

            //rows_orbit_degrees[row_index] = elem;

            int col_index = 0;
            bool first_in_row = true;
            node* curr,*first;

            while( getline(stream, elem, ',') )
            {

                if(!elem.empty()){


                    int bit = std::stoi(elem);

                    if(bit == 1){
                        node* next = new node();
                        next->row = row_index;
                        node* col_head = header_arr[col_index];
                        next->header = col_head;

                        node* bottom_one = col_head->up;

                        bottom_one->down = next;
                        next->up = bottom_one;
                        next->down = col_head;
                        col_head->up = next;
                        
                        if(first_in_row){
                            first_in_row = false;
                            first = next;
                            curr = next;

                            rows_arr[row_index] = next;

                        }
                        else{

                            curr->right = next;
                            next->left = curr;
                            curr = next;
                        }

                        col_head->count++;

                    }

                }


                col_index++;

            }


            if(!first_in_row){
                curr->right = first;
                first->left = curr;
            }


            row_index++;
        }

        myfile.close();

        row_count = row_index;
    }

    void destroy(){
        node* col = dummy->right;
        delete dummy;

        while(col != dummy){

            node* col_elm = col->down;

            while(col_elm != col){

                node* temp_elm = col_elm;
                col_elm = col_elm->down;
                delete temp_elm;
            }

            node* temp_col = col;
            col = col->right;
            delete temp_col;

        }

        delete[] header_arr;
    }

    ~dancing_links(){
        destroy();
    }



};


#endif //MATRIX_BUILDER_DANCING_LINKS_H
