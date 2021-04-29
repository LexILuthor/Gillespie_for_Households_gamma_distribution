//
// Created by Cecilia Meggio on 29/04/21.
//


#include <vector>
#include <list>
#include <iostream>
#include <random>

#include "myClass.h"

//----------------------------------------------------------------------------------------------------------------------


//house
house::house(int nh, int s, int e, int i, parameter par) {
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


house::house(std::vector<int> &state_to_copy,int nh){
    dimension=nh;
    //state.assign(state_to_copy.begin(),state_to_copy.end());
    state=state_to_copy;

}
//----------------------------------------------------------------------------------------------------------------------

//state_to_household_map


state_to_household_map::state_to_household_map(parameter par) {
    std::vector<std::vector<std::vector<std::list<house>>>> tmp(par.nh_max + 1,
                                                                std::vector<std::vector<std::list<house>>>(
                                                                        par.nh_max + 1,
                                                                        std::vector<std::list<house>>(
                                                                                par.nh_max + 1,
                                                                                std::list<house>())));
    matrix = tmp;
    std::vector<std::vector<std::vector<double >>> sumsHiH_state_tmp(par.nh_max + 1,
                                                                 std::vector<std::vector<double >>(
                                                                         par.nh_max + 1,
                                                                         std::vector<double>(
                                                                                 par.nh_max + 1, 0
                                                                         )
                                                                 ));
    sumsHiH_state=sumsHiH_state_tmp;

}


std::list<house> *state_to_household_map::get_households_list(int s, int e, int i) {
    //define a pointer to a list:
    // list<int> *p;
    //p->pop_front() to user the pointer
    return &matrix[s][e][i];
}


void state_to_household_map::add_household(house household, int s, int e, int i) {
    matrix[s][e][i].push_front(household);
    sumsHiH_state[s][e][i] = sumsHiH_state[s][e][i] + (double) s * i / household.dimension;
}

void state_to_household_map::add_household(int nh, int s, int e, int i, parameter par) {
    house household_tmp(nh, s, e, i, par);
    matrix[s][e][i].push_front(household_tmp);
    sumsHiH_state[s][e][i] = sumsHiH_state[s][e][i] + (double) s * i / household_tmp.dimension;
}


void state_to_household_map::add_household(house household, parameter par) {
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
    sumsHiH_state[s][e][i] = sumsHiH_state[s][e][i] + (double) s * i / household.dimension;
}


house state_to_household_map::return_random_household_in_state(int s, int e, int i, parameter &par) {

    // estremi included
    std::uniform_int_distribution<int> random_int(0, matrix[s][e][i].size() - 1);

    int position = random_int(par.generator);
    std::list<house>::iterator it;
    it = matrix[s][e][i].begin();
    advance(it, position);

    house selected_household (it->state,it->dimension);//= *it;
    matrix[s][e][i].erase(it);
    sumsHiH_state[s][e][i] = sumsHiH_state[s][e][i] - (double) s * i / selected_household.dimension;
    return selected_household;

}


house state_to_household_map::select_household_in_state_based_on_infectivity(int s, int e, int i, parameter &par) {

    // estremi included
    std::uniform_real_distribution<double> random(0, sumsHiH_state[s][e][i]);

    double position = random(par.generator);
    std::list<house>::iterator it;
    it = matrix[s][e][i].begin();
    double sum = (double) s * i / it->dimension;
    while (position > sum) {

        advance(it, 1);
        sum = sum + (double) s * i / it->dimension;

    }
    house selected_household = *it;
    matrix[s][e][i].erase(it);
    sumsHiH_state[s][e][i] = sumsHiH_state[s][e][i] - (double) s * i / selected_household.dimension;
    return selected_household;

}



//----------------------------------------------------------------------------------------------------------------------
