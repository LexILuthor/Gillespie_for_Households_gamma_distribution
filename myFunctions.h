//
// Created by popcorn on 01/01/2021.
//

#ifndef GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H
#define GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H

#endif //GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H


class parameter {

public:
    //Number of steps
    int nSteps;

    // number of Households
    int number_of_Households;

    // number of people in one Household
    int nh;

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
};

class house {
public:
    std::vector<int> state;

    //nh
    int dimension;

    house(int nh, int s, int e, int i, parameter par) {
        if (s + e + i > nh) {
            std::cout << "Error! household has dimension only " << dimension << std::endl;
        }
        dimension = nh;
        std::vector<int> tmp(1 + par.number_of_exposed_compartments + par.number_of_infected_compartments, 0);
        state = tmp;
        state[0] = s;
        state[1] = e;
        state[par.number_of_exposed_compartments + 1] = i;
    }
};

class state_to_household_map {
public:
    std::vector<std::vector<std::vector<std::list<house >>>> matrix;

    state_to_household_map(int max_household_dim) {
        std::vector<std::vector<std::vector<std::list<house>>>> tmp(max_household_dim + 1,
                                                                    std::vector<std::vector<std::list<house>>>(
                                                                            max_household_dim + 1,
                                                                            std::vector<std::list<house>>(
                                                                                    max_household_dim + 1,
                                                                                    std::list<house>())));
        matrix = tmp;
    }


    std::list<house> *get_households_list(int s, int e, int i) {
        //define a pointer to a list:
        // list<int> *p;
        //p->pop_front() to user the pointer
        return &matrix[s][e][i];
    }


    void add_household(house household, int s, int e, int i) {
        matrix[s][e][i].push_front(household);
    }

    void add_household(int nh, int s, int e, int i, parameter par) {
        house household_tmp(nh, s, e, i, par);
        matrix[s][e][i].push_front(household_tmp);
    }

    void add_household(house household, parameter par) {
        int s = household.state[0];
        int e = 0;
        int i = 0;
        for (int j = 1; j < par.number_of_exposed_compartments + 1; j++) {
            e = e + household.state[j];
        }
        for (int j = par.number_of_exposed_compartments + 1; j < household.state.size(); j++) {
            i = i + household.state[j];
        }
        matrix[s][e][i].push_front(household);
    }


};

void new_Exposed_outside_the_household(std::vector<std::vector<int> > &SEIR,
                                       state_to_household_map &states_to_household,
                                       int &sumsHiH_nh, parameter &par, int &j);

void new_exposed_inside_the_household(std::vector<std::vector<int> > &SEIR,
                                      state_to_household_map &states_to_household,
                                      int &sumsHiH_nh, parameter &par, int &j);

void new_Infected(std::vector<std::vector<int> > &SEIR,
                  state_to_household_map &states_to_household,
                  int &sumsHiH_nh, parameter &par, int &j);

void new_Recovered(std::vector<std::vector<int> > &SEIR,
                   state_to_household_map &states_to_household,
                   int &sumsHiH_nh, parameter &par, int &j);

void initializeSEIRandTemp(std::vector<std::vector<int>> &SEIR, std::vector<double> &temp, int &N);

void read_Parameters_From_File(std::string inputpath, parameter &parameters);

void write_the_csv_file(std::string outputpath, std::vector<std::vector<int> > &SEIR, std::vector<double> &temp);

void write_lock_down_files(std::string outputpath, std::vector<double> &time_lockdown);

void initialize_household_with_Susceptible_Infected_Exposed(
        std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
        int number_of_Households, int nh);

int initialize_Households(parameter &par, state_to_household_map &states_to_household);


double generateUnif_from_zeroExcluded_to(double to);