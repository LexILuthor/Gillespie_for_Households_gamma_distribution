//
// Created by popcorn on 01/01/2021.
//

#ifndef GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H
#define GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H

#endif //GILLESPIE_FOR_HOUSEHOLDS_MYFUNCTIONS_H


void new_Exposed_outside_the_household(std::vector<std::vector<int>> &SEIR,
                                       std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                                       int &sumsHiH,
                                       std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                                       std::vector<std::vector<int> > &households, int number_of_infected_compartments,
                                       int number_of_exposed_compartments, int &j);

void new_exposed_inside_the_household(std::vector<std::vector<int>> &SEIR,
                                      std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                                      int &sumsHiH,
                                      std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                                      std::vector<std::vector<int> > &households, int number_of_infected_compartments,
                                      int number_of_exposed_compartments, int &j);

void new_Infected(std::vector<std::vector<int> > &SEIR,
                  std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                  int &sumsHiH, std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                  std::vector<std::vector<int> > &households, int number_of_infected_compartments,
                  int number_of_exposed_compartments, int &j);

void new_Recovered(std::vector<std::vector<int>> &SEIR,
                   std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                   int &sumsHiH, std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                   std::vector<std::vector<int> > &households, int number_of_infected_compartments,
                   int number_of_exposed_compartments, int &j);

void initializeSEIRandTemp(std::vector<std::vector<int>> &SEIR, std::vector<double> &temp, int &N);

void read_Parameters_From_File(std::string inputpath, int &nSteps, int &number_of_Households,
                               int &number_of_people_in_one_Household, double &beta1, double &beta2,
                               double &threshold_above_which_one_to_two, double &threshold_under_which_two_to_one,
                               double &betaH, double &ny, double &gamma, int &number_of_infected_compartments,
                               int &number_of_exposed_compartments);

void write_the_csv_file(std::string outputpath, std::vector<std::vector<int> > &SEIR, std::vector<double> &temp);

void write_lock_down_files(std::string outputpath, std::vector<double> &time_lockdown);

void initialize_household_with_Susceptible_Infected_Exposed(
        std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
        int number_of_Households, int number_of_people_in_one_Household);

void initialize_Households(std::vector<std::vector<int> > &households, int number_of_people_in_one_Household,
                           int number_of_exposed_compartments,
                           std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households);


double generateUnif_from_zeroExcluded_to(double to);