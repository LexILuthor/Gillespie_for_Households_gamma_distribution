//
// Created by popcorn on 01/01/2021.
//

#ifndef GILLESPIE_FOR_HOUSEHOLDS_GILLESPIEFORHOUSEHOLDS_H
#define GILLESPIE_FOR_HOUSEHOLDS_GILLESPIEFORHOUSEHOLDS_H

#endif //GILLESPIE_FOR_HOUSEHOLDS_GILLESPIEFORHOUSEHOLDS_H

std::vector<std::vector<int> >
gillespie_for_Households(int nSteps, int N, int number_of_Households, int number_of_people_in_one_Household,
                         double beta1, double beta2, double threshold_above_which_one_to_two,
                         double threshold_under_which_two_to_one, double betaH, double ny, double gamma,
                         int number_of_infected_compartments, int number_of_exposed_compartments,
                         std::vector<double> &temp, std::vector<double> &time_lockdown);