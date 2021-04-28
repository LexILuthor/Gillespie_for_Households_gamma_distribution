//
// Created by popcorn on 18/04/2021.
//


#include <iostream>
#include <vector>
#include<cstdlib>
#include <random>
#include <map>
#include<tuple>
#include <time.h>
#include "myFunctions.h"
#include "GillespieForHouseholds.h"


std::vector<std::vector<int> >
gillespie_for_Households(parameter &par, std::vector<double> &temp, std::vector<double> &time_lockdown) {
    //Here you can change the seed of the generator
    //std::default_random_engine generator(time(0));
    //std::default_random_engine generator(0);
    std::default_random_engine generator;
    //srand(0);
    //srand(time(0));


    // vector containing al the households, each cell will be of the shape (s,e1,e2,e2,...,i1,i2,i3...)
    std::vector<house> households(par.number_of_Households);

    // map that connects each possible state (s,e,i,r) with the households that are in that state

    state_to_household_map states_to_household(par);


    initialize_Households(households, par, states_to_household);

    std::vector<std::vector<int> > SEIR(4, std::vector<int>(1, 0));
    double move = (double) 1 / par.N;

    //setting the initial conditions with N-1 susceptible 1 infected and zero exposed and recovered
    initializeSEIRandTemp(SEIR, temp, par.N);

    //the sum of (number of susceptible)*(number of infected) over all household
    int sumsHiH = par.nh - 1;


    //the matrix "household_with_Susceptible_Infected_Exposed" will contain in the position [s,i,e] the number of households with i susceptible people, i infected and e exposed
    std::vector<std::vector<std::vector<int>>> household_with_Susceptible_Infected_Exposed(
            par.nh + 1, std::vector<std::vector<int>>(par.nh + 1, std::vector<int>(par.nh + 1, 0)));

    //setting the initial condition of the matrix "household_with_Susceptible_Infected_Exposed"
    initialize_household_with_Susceptible_Infected_Exposed(household_with_Susceptible_Infected_Exposed,
                                                           par.number_of_Households, par.nh);

    par.beta = par.beta1;

    std::exponential_distribution<double> exp_distribution(1);
    std::uniform_real_distribution<double> uniform_Real_Distribution(0.0, 1.0);


    // here we simulate the process
    int j = 1;
    while (j < par.nSteps) {
        //number of Susceptible
        int s = SEIR[0][j - 1];

        //number of Exposed
        int e = SEIR[1][j - 1];

        //Number of Infected
        int i = SEIR[2][j - 1];

        //number of Recovered
        int r = SEIR[3][j - 1];


        if (r == par.N) {
            //everyone recovered
            return SEIR;
        }

        //change beta when we have 10% of the population recovered

        if (i >= (par.N / 100) * par.threshold_above_which_one_to_two && par.beta != par.beta2) {
            par.beta = par.beta2;
            std::cout << "beta decrease at time t= " << temp.back() << "\n";
            time_lockdown.push_back(temp.back());
        } else if (i < (par.N / 100) * par.threshold_under_which_two_to_one && par.beta != par.beta1) {
            par.beta = par.beta1;
            std::cout << "beta increase at time t= " << temp.back() << "\n";
            time_lockdown.push_back(temp.back());
        }


        // compute the parameter lambda of the exponential and the probabilities of
        // S->E, E->I, I->R

        double se = par.beta * s * i * move;
        double seH = par.betaH * sumsHiH / par.nh;
        double ei = par.ny * e;
        double ir = par.gamma * i;
        double lambda = (se + seH + ei + ir);
        if (lambda == 0) {
            break;
        }
        se = se / lambda;
        seH = seH / lambda;
        ei = ei / lambda;
        ir = ir / lambda;


        //generate the time of the next event with an exponential with parameter lambda
        double event = exp_distribution(generator);


        event = event / lambda;
        temp.push_back(temp.back() + event);


        //Randomly decide which event happened
        //double tmp = rand() / ((double) RAND_MAX);
        double tmp = uniform_Real_Distribution(generator);

        if (tmp < se) {
            //new Exposed from a contact outside the household
            new_Exposed_outside_the_household(SEIR, household_with_Susceptible_Infected_Exposed, sumsHiH,
                                              states_to_households, households,
                                              par, j);
        } else if (tmp < (se + seH)) {
            //new Exposed from a contact within the household
            new_exposed_inside_the_household(SEIR, household_with_Susceptible_Infected_Exposed, sumsHiH,
                                             states_to_households, households, par, j);
        } else if (tmp < (se + seH + ei)) {
            //new infected
            new_Infected(SEIR, household_with_Susceptible_Infected_Exposed, sumsHiH, states_to_households, households,
                         par, j);
        } else {
            //new Recovered
            new_Recovered(SEIR, household_with_Susceptible_Infected_Exposed, sumsHiH, states_to_households, households,
                          par, j);
        }
        j++;
    }

    return SEIR;

}