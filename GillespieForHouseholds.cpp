//
// Created by popcorn on 18/04/2021.
//


#include <iostream>
#include <vector>
#include<cstdlib>
#include <random>
#include <time.h>
#include <map>
#include<tuple>
#include "myFunctions.h"
#include "GillespieForHouseholds.h"

std::vector<std::vector<int> >
gillespie_for_Households(int nSteps, int N, int number_of_Households, int number_of_people_in_one_Household,
                         double beta1, double beta2, double threshold_above_which_one_to_two,
                         double threshold_under_which_two_to_one, double betaH, double ny, double gamma,
                         int number_of_infected_compartments, int number_of_exposed_compartments,
                         std::vector<double> &temp, std::vector<double> &time_lockdown) {
    //Here you can change the seed of the generator
    std::default_random_engine generator(time(0));
    //srand(time(0));
    //std::default_random_engine generator;

    // vector containing al the households, each cell will be of the shape (s,e1,e2,e2,...,i1,i2,i3...)
    std::vector<std::vector<int> > households(number_of_Households, std::vector<int>(
            1 + number_of_exposed_compartments + number_of_infected_compartments, 0));

    // map that connects each possible state (s,e,i,r) with the households that are in that state
    std::map<std::tuple<int, int, int>, std::vector<int> > states_to_households;

    initialize_Households(households, number_of_people_in_one_Household, number_of_exposed_compartments,
                          states_to_households);

    std::vector<std::vector<int> > SEIR(4, std::vector<int>(1, 0));
    double move = (double) 1 / N;

    //setting the initial conditions with N-1 susceptible 1 infected and zero exposed and recovered
    initializeSEIRandTemp(SEIR, temp, N);

    //the sum of (number of susceptible)*(number of infected) over all household
    int sumsHiH = number_of_people_in_one_Household - 1;


    //the matrix "household_with_Susceptible_Infected_Exposed" will contain in the position [s,i,e] the number of households with i susceptible people, i infected and e exposed
    std::vector<std::vector<std::vector<int>>> household_with_Susceptible_Infected_Exposed(
            number_of_people_in_one_Household + 1, std::vector<std::vector<int>>(number_of_people_in_one_Household + 1,
                                                                                 std::vector<int>(
                                                                                         number_of_people_in_one_Household +
                                                                                         1, 0)));

    //setting the initial condition of the matrix "household_with_Susceptible_Infected_Exposed"
    initialize_household_with_Susceptible_Infected_Exposed(household_with_Susceptible_Infected_Exposed,
                                                           number_of_Households, number_of_people_in_one_Household);

    double beta = beta1;


    // here we simulate the process
    int j = 1;
    while (j < nSteps) {
        //number of Susceptible
        int s = SEIR[0][j - 1];

        //number of Exposed
        int e = SEIR[1][j - 1];

        //Number of Infected
        int i = SEIR[2][j - 1];

        //number of Recovered
        int r = SEIR[3][j - 1];


        if (r == N) {
            //everyone recovered
            return SEIR;
        }

        //change beta when we have 10% of the population recovered

        if (i >= (N / 100) * threshold_above_which_one_to_two && beta != beta2) {
            beta = beta2;
            std::cout << "beta decrease at time t= " << temp.back() << "\n";
            time_lockdown.push_back(temp.back());
        } else if (i < (N / 100) * threshold_under_which_two_to_one && beta != beta1) {
            beta = beta1;
            std::cout << "beta increase at time t= " << temp.back() << "\n";
            time_lockdown.push_back(temp.back());
        }


        // compute the parameter lambda of the exponential and the probabilities of
        // S->E, E->I, I->R

        double se = beta * s * i * move;
        double seH = betaH * sumsHiH / number_of_people_in_one_Household;
        double ei = ny * e;
        double ir = gamma * i;
        double lambda = (se + seH + ei + ir);
        if (lambda == 0) {
            break;
        }
        se = se / lambda;
        seH = seH / lambda;
        ei = ei / lambda;
        ir = ir / lambda;


        //generate the time of the next event with an exponential with parameter lambda
        std::exponential_distribution<double> exp_distribution(lambda);
        double event = exp_distribution(generator);
        temp.push_back(temp.back() + event);


        //Randomly decide which event happened
        double tmp = rand() / ((double) RAND_MAX + 1);
        if (tmp < se) {
            //new Exposed from a contact outside the household
            new_Exposed_outside_the_household(SEIR, household_with_Susceptible_Infected_Exposed, sumsHiH, j);
        } else if (tmp < (se + seH)) {
            //new Exposed from a contact within the household
            new_exposed_inside_the_household(SEIR, household_with_Susceptible_Infected_Exposed, sumsHiH, j);
        } else if (tmp < (se + seH + ei)) {
            //new infected
            new_Infected(SEIR, household_with_Susceptible_Infected_Exposed, sumsHiH, states_to_households, households,
                         number_of_infected_compartments, number_of_exposed_compartments, j);
        } else {
            //new Recovered
            new_Recovered(SEIR, household_with_Susceptible_Infected_Exposed, sumsHiH, states_to_households, households,
                          number_of_infected_compartments, number_of_exposed_compartments, j);
        }
        j++;
    }

    return SEIR;

}