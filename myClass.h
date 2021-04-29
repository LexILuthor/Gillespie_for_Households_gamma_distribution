//
// Created by Cecilia Meggio on 29/04/21.
//

#ifndef GILLESPIE_FOR_HOUSEHOLDS_GAMMA_DISTRIBUTION_MYCLASS_H
#define GILLESPIE_FOR_HOUSEHOLDS_GAMMA_DISTRIBUTION_MYCLASS_H

#endif //GILLESPIE_FOR_HOUSEHOLDS_GAMMA_DISTRIBUTION_MYCLASS_H

#include <list>

class parameter {

public:
    //Number of steps
    int nSteps;

    // number of Households
    int number_of_Households;

    // number of people in one Household
    int nh_mean;

    // number of people in one Household
    int nh_max;

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

    int N;

    double beta;

    std::default_random_engine generator;

    void initialize_generator(unsigned seed) {
        std::default_random_engine tmp(seed);
        generator = tmp;
    }
};

class house {
public:
    std::vector<int> state;

    //nh
    int dimension;

    house(int nh, int s, int e, int i, parameter par);
    house(std::vector<int> &state_to_copy,int nh);
};

class state_to_household_map {
public:
    std::vector<std::vector<std::vector<std::list<house >>>> matrix;

    state_to_household_map(parameter par);


    std::vector<std::vector<std::vector<double >>> sumsHiH_state;

    std::list<house> *get_households_list(int s, int e, int i);


    void add_household(house household, int s, int e, int i);

    void add_household(int nh, int s, int e, int i, parameter par);

    void add_household(house household, parameter par);

    house return_random_household_in_state(int s, int e, int i, parameter &par);

    house select_household_in_state_based_on_infectivity(int s, int e, int i, parameter &par);


};