//
// Created by User on 4/27/2019.
//


#include <sstream>
#include "field_aux.h"
#include <sys/stat.h>

void markOrbit(Field& field, int deg, int orbit){
    field[deg].push_back(std::to_string(orbit));
}

bool checkMark(Field& field, int deg){
    return field[deg].size() > 1;
}

bool checkDifferences(Field& field, int diffs[]){
    //TODO: Assumes diffs is of size 7. Consider passing size as parameter.
    for(int i=0; i < 7;i++){
        for(int j=i+1; j < 7;j++){
            if(field[diffs[i]][1] == field[diffs[j]][1]){
                return false;
            }
        }
    }
    return true;
}

//Generates GF(2, 19) to external file
int field_generator(Field& field) {



    std::cout << "Field generation started\n";
    std::ofstream myfile;
    myfile.open ("field.csv");
    if(myfile.is_open() == 0)
        return 1;
    myfile << "@," << std::bitset<19>(0) << std::endl;

    int elem = 0b1;
    //using polynomial x^19 + x^5 + x^2 + x^1 + 1
    for(int i=0;i<524287;i++){
        myfile << i << "," << std::bitset<19>(elem) << std::endl;
        field[i].push_back(std::bitset<19>(elem).to_string());
        if((elem&0b1000000000000000000) == 0b1000000000000000000){
            elem = (((elem&(0b1000000000000000000 - 1))<<1)^0b100111);
        }
        else{
            elem = (elem<<1);
        }
    }
    myfile.close();
    std::cout << "Field successfully generated\n";
    return 0;
}

int field_update(Field& field){
    std::cout << "Field update started\n";
    std::ofstream myfile;
    myfile.open ("field2.csv");
    if(myfile.is_open() == 0)
        return 1;
    myfile << "@," << std::bitset<19>(0) << std::endl;

    int elem = 0b1;
    //using polynomial x^19 + x^5 + x^2 + x^1 + 1
    for(int i=0;i < field.size();i++){
        if(checkMark(field,i)){
            myfile << i << "," << std::bitset<19>(field[i][0]) << "," << field[i][1] << std::endl;
        }
        else{
            myfile << i << "," << std::bitset<19>(field[i][0]) << std::endl;
        }


    }
    myfile.close();

    std::cout << "Field successfully updated\n";
    return 0;
}

bool check_SS2_integrity(){

    std::string filename = "SS2_orbits.txt";
    struct stat buffer;
    if((stat (filename.c_str(), &buffer) != 0)){
        std::cout << "file SS2_orbits.txt does not exist" << std::endl;
        return false;
    }


    std::vector<std::vector<int>> output;
    std::ifstream myfile ("SS2_orbits.txt");
    std::string line;
    int line_count = 0;
    while(getline(myfile,line)){
        line_count++;
        std::istringstream stream(line);
        std::string elem;
        int comma_count = 0;
        while( getline(stream, elem, ',') )
        {
            comma_count++;
            if(!elem.empty()){
                int deg;
                try {
                    deg = std::stoi(elem);
                }
                catch(std::logic_error& e){
                    std::cout << "one of the degrees in file SS2_orbits.txt is not an acceptable integer" << std::endl;
                    return false;
                }

                if(deg <= 0 || deg >= pow(2,19)-1){
                    std::cout << "one of the degrees in file SS2_orbits.txt is not in the range [1,2^19-2]" << std::endl;
                    return false;
                }
            }
        }

        if(comma_count != 114){
            std::cout << "file SS2_orbits.txt should have 114 degrees in each orbit, each followed by a comma" << std::endl;
            return false;
        }
    }

    if(line_count != 4599){
        std::cout << "file SS2_orbits.txt should have 4599 consecutive non-empty lines of orbits." << std::endl;
        return false;
    }

    return true;

}


std::vector<std::vector<int>> orbits_file_to_SS2_matrix(Field& field){

    std::vector<std::vector<int>> output;
    std::ifstream myfile ("SS2_orbits.txt");
    std::string line;
    while(getline(myfile,line)){
        std::istringstream stream(line);
        std::vector<int> row;
        std::string elem;

        while( getline(stream, elem, ',') )
        {
            if(!elem.empty()){
                row.push_back(std::stoi(elem));

                field[std::stoi(elem)].push_back(std::to_string(output.size()));
            }

        }
        output.push_back(row);
    }
    return output;
}


bool check_field2_integrity(){
    std::string filename = "field2.csv";

    struct stat buffer;
    if((stat (filename.c_str(), &buffer) != 0)){
        std::cout << "file field2.csv does not exist" << std::endl;
        return false;
    }

    std::ifstream myfile ("field2.csv");
    std::string line;
    getline(myfile,line);// the line with @, no need for it and it screws up the code below since @ is not an int
    int line_count = 1;
    bool flag_that_turns_on_when_getting_to_a_1 = false;
    while(getline(myfile,line)){
        line_count++;
        std::istringstream stream(line);

        std::string elem;

        getline(stream, elem, ',');
        int deg;
        try {
            deg = std::stoi(elem);
        }
        catch(std::logic_error& e){
            std::cout << "one of the degrees in file field2.csv is not an acceptable integer" << std::endl;
            return false;
        }

        if(deg < 0 || deg >= pow(2,19)-1){
            std::cout << "one of the degrees in file field2.csv is not in the range [0,2^19-2]" << std::endl;
            return false;
        }


        getline(stream, elem, ',');
        std::string bit_rep = elem;

        if(bit_rep.size() != 19){
            std::cout << "a bit representation of one of the orbit degrees in not of size 19" << std::endl;
            return false;
        }
        std::bitset<19> bit_set;
        try{
            bit_set = std::bitset<19>(bit_rep);
        }
        catch(std::invalid_argument& e){
            std::cout << "a bit representation of one of the orbit degrees does not consist purely of ones and zeros" << std::endl;
            return false;
        }

        if(getline(stream, elem, ',') && flag_that_turns_on_when_getting_to_a_1){
            std::string orbit_rep = elem;
            int deg_orbit;
            try {
                deg_orbit = std::stoi(elem);
            }
            catch(std::logic_error& e){
                std::cout << "one of the orbit representatives in file field2.csv is not an acceptable integer" << std::endl;
                return false;
            }

            if(deg_orbit < 0 || deg_orbit >= pow(2,19)-1){
                std::cout << "one of the orbit representatives in file field2.csv is not in the range [0,2^19-2]" << std::endl;
                return false;
            }

        }
        flag_that_turns_on_when_getting_to_a_1 = true;
    }

    if(line_count != pow(2,19)){
        std::cout << "file field2.csv should have" << pow(2,19) << "consecutive non-empty lines starting from the second line" << std::endl;
        return false;
    }

    return true;

}

void field2_to_Field(Field& field, Reverse_Field& reverse_field){


    std::ifstream myfile ("field2.csv");
    std::string line;
    getline(myfile,line);// the line with @, no need for it and it screws up the code below since @ is not an int
    while(getline(myfile,line)){
        std::istringstream stream(line);
        std::vector<std::string> row;
        std::vector<std::string> other_row;
        std::string elem;

        getline(stream, elem, ',');
        int deg = std::stoi(elem);
        other_row.push_back(elem);
        getline(stream, elem, ',');
        std::string bit_rep = elem;
        row.push_back(bit_rep);
        if(getline(stream, elem, ',')){
            std::string orbit_rep = elem;
            row.push_back(orbit_rep);
            other_row.push_back(orbit_rep);
        }
        field[deg] = row;
        reverse_field[bit_rep] = other_row;
    }


}


