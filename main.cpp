#include <iostream>
#include <vector>
#include <map>
#include <random>
#include <list>


#include "myClass.h"
#include "myFunctions.h"
#include "GillespieForHouseholds.h"


int main() {


    std::string inputpath = "../Input/Input_Gillespie_Household.txt";
    std::string outputpath = "../Output/gillespie_Household";


    int tot_simulations = 100;

    parameter par;
    read_Parameters_From_File(inputpath, par);



    //par.initialize_generator(seed=0);

    if (par.beta1 != par.beta2) {
        outputpath = "../Output/gillespie_Household_lockdown";
    }



    std::random_device myRandomDevice;
    //unsigned seed = myRandomDevice();
    unsigned seed = 10;
    par.initialize_generator(seed);



    // Gillespie algorithm.
    // SEIR is a matrix that contains the data relative to the number of infected recovred etc..
    // tmp is the vector of the time (each entry is the time at which an event occurred)
    for (int i = 0; i < tot_simulations; i++) {
        std::vector<double> tempo;
        std::vector<double> time_lockdown;
        std::vector<std::vector<int> > SEIR = gillespie_for_Households(par, tempo,
                                                                       time_lockdown);

        if (par.beta1 != par.beta2) {
            write_lock_down_files(outputpath + std::to_string(i) + "lock_down_time" + ".txt", time_lockdown);
        }
        write_the_csv_file(outputpath + std::to_string(i) + ".csv", SEIR, tempo);
        std::cout << i << "\n";
    }
    return 0;

}