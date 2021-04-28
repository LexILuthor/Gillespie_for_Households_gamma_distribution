//
// Created by popcorn on 01/01/2021.
//

#ifndef GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H
#define GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H

#endif //GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H


class house {
public:
    std::vector<int> state;
    int dimension;
};




class parameter {

public:
    //Number of steps
    int nSteps;

    // number of Households
    int number_of_Households;

    // number of people in one Household
    int nh;

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
};

class state_to_household_map {
public:
    std::vector<std::vector<std::vector<std::list<house>>>> matrix;
    state_to_household_map(parameter par){
        std::vector<std::vector<std::vector<std::list<house>>>> tmp(par.nh + 1, std::vector<std::vector<std::list<house>>>(
                par.nh + 1, std::vector<std::list<house>>(par.nh + 1, std::list<house>())));
        matrix= tmp;
    }



    std::list<house> get_households_list(int s, int e, int i){
        return matrix[s][e][i];
    }

};

void new_Exposed_outside_the_household(std::vector<std::vector<int>> &SEIR,
                                       std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                                       int &sumsHiH,
                                       std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                                       std::vector<std::vector<int> > &households, parameter &par, int &j);

void new_exposed_inside_the_household(std::vector<std::vector<int>> &SEIR,
                                      std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                                      int &sumsHiH,
                                      std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                                      std::vector<std::vector<int> > &households, parameter &par, int &j);

void new_Infected(std::vector<std::vector<int> > &SEIR,
                  std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                  int &sumsHiH, std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                  std::vector<std::vector<int> > &households, parameter &par, int &j);

void new_Recovered(std::vector<std::vector<int>> &SEIR,
                   std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                   int &sumsHiH, std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                   std::vector<std::vector<int> > &households, parameter &par, int &j);

void initializeSEIRandTemp(std::vector<std::vector<int>> &SEIR, std::vector<double> &temp, int &N);

void read_Parameters_From_File(std::string inputpath, parameter &parameters);

void write_the_csv_file(std::string outputpath, std::vector<std::vector<int> > &SEIR, std::vector<double> &temp);

void write_lock_down_files(std::string outputpath, std::vector<double> &time_lockdown);

void initialize_household_with_Susceptible_Infected_Exposed(
        std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
        int number_of_Households, int nh);

void initialize_Households(std::vector<house> &households, parameter &par,
                           std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households);


double generateUnif_from_zeroExcluded_to(double to);