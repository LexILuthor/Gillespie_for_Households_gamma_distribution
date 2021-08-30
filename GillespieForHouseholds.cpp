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

#include "myClass.h"
#include "myFunctions.h"
#include "GillespieForHouseholds.h"


std::vector<std::vector<int> >
gillespie_for_Households(parameter &par, std::vector<double> &temp, std::vector<double> &time_lockdown) {
    //Here you can change the seed of the generator
    //std::default_random_engine generator(time(0));
    //std::default_random_engine generator(0);
    //srand(0);
    //srand(time(0));



    // map that connects each possible state (s,e,i,r) with the households that are in that state
    state_to_household_map states_to_household(par);

    //the sum of (number of susceptible)*(number of infected)/nh over all household
    double sumsHiH_nh = initialize_Households(par, states_to_household);




    //setting the initial conditions with N-1 susceptible 1 infected and zero exposed and recovered
    std::vector<std::vector<int> > SEIR(4, std::vector<int>(1, 0));
    initializeSEIRandTemp(SEIR, temp, par.N);


    par.beta = par.beta1;

    std::exponential_distribution<double> exp_distribution(1);
    std::uniform_real_distribution<double> uniform_Real_Distribution(0.0, 1.0);


    double move = (double) 1 / par.N;

    // here we simulate the process
    int daily_infected = 0;
    int epidemic_day = 0;
    bool sync_day_reaced = false;
    int lockdown_day = 9999999;


    int j = 1;
    while (j < par.nSteps) {

        if (j % 1000 == 0) {
            std::cout << j << std::endl;
        }

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

        if (temp[j - 1] > epidemic_day) {
            epidemic_day = floor(temp[j - 1]) + 1;
            //std::cout << "daily infected: " << daily_infected << "\n";
            daily_infected = 0;
        }





        //activate lockdown par.lockdown_delay days after we reach par.daily_infected_sync new daily infected
        if (!sync_day_reaced && daily_infected > par.daily_infected_sync) {
            sync_day_reaced = true;
            std::cout << "sync day reached\n";
            lockdown_day = int(temp[j - 1] + par.lockdown_delay);
        }
        if (sync_day_reaced && temp[j - 1] > lockdown_day) {
            //std::cout << "lockdown starts\n";
            par.beta = par.beta2;
        }


        //activate lockdown at time par.time_activate_lockdown
        //if (temp[j - 1] > par.time_activate_lockdown) {
        //    par.beta = par.beta2;
        //}
        //if (temp[j - 1] > par.time_end_lockdown) {
        //    par.beta = par.beta3;
        //}


        //change beta when we have threshold_above_which_one_to_two % of the population infected
        // or les than threshold_under_which_two_to_one % is infected
        /*
        if (e >= ((double) par.N / 100) * par.threshold_above_which_one_to_two && par.beta != par.beta2) {
            par.beta = par.beta2;
            std::cout << "beta decrease at time t= " << temp.back() << "\n";
            time_lockdown.push_back(temp.back());
        } else if (e < ((double) par.N / 100) * par.threshold_under_which_two_to_one && par.beta != par.beta1) {
            par.beta = par.beta1;
            std::cout << "beta increase at time t= " << temp.back() << "\n";
            time_lockdown.push_back(temp.back());
        }
        */

        //sumsHiH_nh gives some problems of approximation.
        // since it is a double it may happen that it does not go to zero
        // here i check if it should be zero
        if (sumsHiH_nh < (0.1 / par.nh_max) || (i == 0 && s == 0)) {
            sumsHiH_nh = 0;
        }

        // compute the parameter lambda of the exponential and the probabilities of
        // S->E, E->I, I->R
        double se = par.beta * s * i * move;
        double seH = par.betaH * sumsHiH_nh;
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
        double event = exp_distribution(par.generator);
        event = event / lambda;
        temp.push_back(temp.back() + event);


        //Randomly decide which event happened
        //double tmp = rand() / ((double) RAND_MAX);
        double tmp = uniform_Real_Distribution(par.generator);

        if (tmp < se) {
            //new Exposed from a contact outside the household
            new_Exposed_outside_the_household(SEIR,
                                              states_to_household,
                                              sumsHiH_nh,
                                              par, j);


        } else if (tmp < (se + seH)) {
            //new Exposed from a contact within the household
            new_exposed_inside_the_household(SEIR,
                                             states_to_household,
                                             sumsHiH_nh,
                                             par, j);
        } else if (tmp < (se + seH + ei)) {
            //new infected


            double tmp_SumsHiH;
            new_Infected(SEIR,
                         states_to_household,
                         sumsHiH_nh,
                         par, j);
            if (tmp_SumsHiH != sumsHiH_nh) {
                //if sumsHiH has changes it means we actually have a new infected
                daily_infected++;
            }

        } else {
            //new Recovered
            new_Recovered(SEIR,
                          states_to_household,
                          sumsHiH_nh,
                          par, j);
        }
        j++;
    }

    return SEIR;

}