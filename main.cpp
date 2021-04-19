#include <iostream>
#include <vector>
#include "myFunctions.h"
#include "GillespieForHouseholds.h"

int main() {


    std::string inputpath = "../Input/Input_Gillespie_Household.txt";
    std::string outputpath = "../Output/gillespie_Household";


    //Number of steps
    int nSteps;

    // number of Households
    int number_of_Households;

    // number of people in one Household
    int number_of_people_in_one_Household;

    //S->E the initial beta
    double beta1;

    //S->E after the lockdown (automatically activated when a certain % of the population is recovered)
    double beta2;

    //the % above which the general contacts happens at ratio beta2
    double threshold_above_which_one_to_two;

    //the % under which the general contacts happens at ratio beta1
    double threshold_under_which_two_to_one;

    //S->E in the household
    double betaH;

    // E-> I
    double ny;

    // I->R
    double gamma;


    int number_of_infected_compartments;

    int number_of_exposed_compartments;

    int tot_simulations = 100;


    read_Parameters_From_File(inputpath, nSteps, number_of_Households, number_of_people_in_one_Household, beta1, beta2,
                              threshold_above_which_one_to_two, threshold_under_which_two_to_one,
                              betaH, ny, gamma, number_of_infected_compartments,
                              number_of_exposed_compartments);

    if (beta1 != beta2) {
        outputpath = "../Output/gillespie_Household_lockdown";
    }

    int N = number_of_Households * number_of_people_in_one_Household;


    // Gillespie algorithm.
    // SEIR is a matrix that contains the data relative to the number of infected recovred etc..
    // tmp is the vector of the time (each entry is the time at which an event occurred)
    for (int i = 0; i < tot_simulations; i++) {
        std::vector<double> tempo;
        std::vector<double> time_lockdown;
        std::vector<std::vector<int> > SEIR = gillespie_for_Households(nSteps, N, number_of_Households,
                                                                       number_of_people_in_one_Household, beta1, beta2,
                                                                       threshold_above_which_one_to_two,
                                                                       threshold_under_which_two_to_one, betaH, ny,
                                                                       gamma, number_of_infected_compartments,
                                                                       number_of_exposed_compartments, tempo,
                                                                       time_lockdown);

        if (beta1 != beta2) {
            write_lock_down_files(outputpath + std::to_string(i) + "lock_down_time" + ".txt", time_lockdown);
        }
        write_the_csv_file(outputpath + std::to_string(i) + ".csv", SEIR, tempo);
        std::cout << i << "\n";
    }
    return 0;

}