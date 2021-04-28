#include <iostream>
#include <vector>
#include <map>
#include "myFunctions.h"
#include "GillespieForHouseholds.h"




int main() {

    parameter parameters;

    std::string inputpath = "../Input/Input_Gillespie_Household.txt";
    std::string outputpath = "../Output/gillespie_Household";




    int tot_simulations = 100;


    read_Parameters_From_File(inputpath, parameters);
    parameters.nh_max=parameters.nh;

    if (parameters.beta1 != parameters.beta2) {
        outputpath = "../Output/gillespie_Household_lockdown";
    }

    parameters.N = parameters.number_of_Households * parameters.nh;


    // Gillespie algorithm.
    // SEIR is a matrix that contains the data relative to the number of infected recovred etc..
    // tmp is the vector of the time (each entry is the time at which an event occurred)
    for (int i = 0; i < tot_simulations; i++) {
        std::vector<double> tempo;
        std::vector<double> time_lockdown;
        std::vector<std::vector<int> > SEIR = gillespie_for_Households(parameters, tempo,
                                                                       time_lockdown);

        if (parameters.beta1 != parameters.beta2) {
            write_lock_down_files(outputpath + std::to_string(i) + "lock_down_time" + ".txt", time_lockdown);
        }
        write_the_csv_file(outputpath + std::to_string(i) + ".csv", SEIR, tempo);
        std::cout << i << "\n";
    }
    return 0;

}