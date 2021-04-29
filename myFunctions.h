//
// Created by popcorn on 01/01/2021.
//

#ifndef GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H
#define GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H

#endif //GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H


void new_Exposed_outside_the_household(std::vector<std::vector<int> > &SEIR,
                                       state_to_household_map &states_to_household,
                                       double &sumsHiH_nh, parameter &par, int &j);

void new_exposed_inside_the_household(std::vector<std::vector<int> > &SEIR,
                                      state_to_household_map &states_to_household,
                                      double &sumsHiH_nh, parameter &par, int &j);

void new_Infected(std::vector<std::vector<int> > &SEIR,
                  state_to_household_map &states_to_household,
                  double &sumsHiH_nh, parameter &par, int &j);

void new_Recovered(std::vector<std::vector<int> > &SEIR,
                   state_to_household_map &states_to_household,
                   double &sumsHiH_nh, parameter &par, int &j);

void initializeSEIRandTemp(std::vector<std::vector<int>> &SEIR, std::vector<double> &temp, int &N);

void read_Parameters_From_File(std::string inputpath, parameter &parameters);

void write_the_csv_file(std::string outputpath, std::vector<std::vector<int> > &SEIR, std::vector<double> &temp);

void write_lock_down_files(std::string outputpath, std::vector<double> &time_lockdown);

void initialize_household_with_Susceptible_Infected_Exposed(
        std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
        int number_of_Households, int nh);

double initialize_Households(parameter &par, state_to_household_map &states_to_household);


double generateUnif_from_zeroExcluded_to(double to, parameter &par);