//
// Created by popcorn on 18/04/2021.
//


#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <map>
#include<tuple>
#include <list>
#include <time.h>
#include <math.h>

#include "myClass.h"
#include "myFunctions.h"


void new_Exposed_outside_the_household(std::vector<std::vector<int> > &SEIR,
                                       state_to_household_map &states_to_household,
                                       double &sumsHiH_nh, parameter &par, int &j) {
    SEIR[0].push_back(SEIR[0][j - 1] - 1);
    SEIR[1].push_back(SEIR[1][j - 1] + 1);
    SEIR[2].push_back(SEIR[2][j - 1]);
    SEIR[3].push_back(SEIR[3][j - 1]);

    //update households with susceptible based only on how many susceptible an house has
    //update also sumsHiH
    // generate a random number and decide which household will change



    std::uniform_int_distribution<int> random_int(1, SEIR[0][j - 1]);
    int random_i = random_int(par.generator);

    //int size = par.nh_max;
    int cumulativeSum = 0;

    for (int e = 0; e <= par.nh_max; e++) {

        for (int i = 0; i <= par.nh_max - e; i++) {

            for (int s = 0; s <= par.nh_max - (e + i); s++) {
                cumulativeSum = cumulativeSum + ((int) states_to_household.matrix[s][e][i].size() * s);
                if (random_i <= cumulativeSum) {
                    //allora abbiamo estratto il numero (s,e,i)

                    //we have a new exposed

                    house selected_household = states_to_household.return_random_household_in_state(s, e, i, par);
                    selected_household.state[0]--;
                    selected_household.state[1]++;

                    states_to_household.add_household(selected_household, s - 1, e + 1, i);

                    // this is the rewrite of:
                    // sumsHiH = sumsHiH - (s * i)/nh + ((s - 1) * i )/nh

                    sumsHiH_nh = sumsHiH_nh - (double) i / selected_household.dimension;



                    //--------------------------------------------------------------------------------------------------


                    goto skip;

                }
            }
        }
    }
    skip:;
}


void new_exposed_inside_the_household(std::vector<std::vector<int> > &SEIR,
                                      state_to_household_map &states_to_household,
                                      double &sumsHiH_nh, parameter &par, int &j) {
    SEIR[0].push_back(SEIR[0][j - 1] - 1);
    SEIR[1].push_back(SEIR[1][j - 1] + 1);
    SEIR[2].push_back(SEIR[2][j - 1]);
    SEIR[3].push_back(SEIR[3][j - 1]);

    //update households with susceptible based on how the quantity susceptible * infected/nh
    //update also sumsHiH
    // generate a random number and decide which household will change
    //double randomUnif = generateUnif_from_zeroExcluded_to(sumsHiH_nh, par);

    std::uniform_real_distribution<double> random_real(0, sumsHiH_nh);
    double randomUnif = random_real(par.generator);


    double cumulativeSum = 0;

    for (int e = 0; e <= par.nh_max; e++) {

        for (int i = 0; i <= par.nh_max - e; i++) {

            for (int s = 0; s <= par.nh_max - (e + i); s++) {
                cumulativeSum = cumulativeSum + (states_to_household.sumsHiH_state[s][e][i]);
                if (randomUnif <= cumulativeSum) {

                    //allora abbiamo estratto il numero (s,e,i)

                    house selected_household = states_to_household.select_household_in_state_based_on_infectivity(
                            s, e, i, par);

                    selected_household.state[0]--;
                    selected_household.state[1]++;

                    states_to_household.add_household(selected_household, s - 1, e + 1, i);

                    // this is the rewrite of:
                    // sumsHiH = sumsHiH - (s * i)/nh + ((s - 1) * i )/nh

                    sumsHiH_nh = sumsHiH_nh - (double) i / selected_household.dimension;



                    //-------------------------------------------------------------------------------------------------


                    goto skip;


                }

            }
        }
    }
    skip:;

}

void new_Infected(std::vector<std::vector<int> > &SEIR,
                  state_to_household_map &states_to_household,
                  double &sumsHiH_nh, parameter &par, int &j) {


    //update households with susceptible based only on how many exposed an house has
    //Note it may happen that thi is only a passage rome an exposed compartment to another and no new infected are created
    //update also sumsHiH




    // generate a random number and decide which household will change

    std::uniform_int_distribution<int> random_int(1, SEIR[1][j - 1]);
    int random_i = random_int(par.generator);

    int cumulativeSum = 0;

    for (int e = 0; e <= par.nh_max; e++) {

        for (int i = 0; i <= par.nh_max - e; i++) {

            for (int s = 0; s <= par.nh_max - (e + i); s++) {
                cumulativeSum = cumulativeSum + ((int) states_to_household.matrix[s][e][i].size() * e);
                if (random_i <= cumulativeSum) {
                    //allora abbiamo estratto il numero (s, e, i)

                    //move exposed

                    house selected_household = states_to_household.return_random_household_in_state(s, e, i, par);

                    // decide in which compartment happens the change based on how many exposed each compartment has

                    std::uniform_int_distribution<int> random_int_compartment(1, e);
                    int random_compartment_e = random_int_compartment(par.generator);

                    int compartment = 1;
                    int cum_sum = selected_household.state[compartment];
                    while (cum_sum < random_compartment_e) {
                        compartment++;
                        cum_sum = cum_sum + selected_household.state[compartment];
                    }

                    selected_household.state[compartment]--;
                    selected_household.state[compartment + 1]++;

                    if (compartment == par.number_of_exposed_compartments) {
                        //abbiamo un nuovo infetto

                        SEIR[0].push_back(SEIR[0][j - 1]);
                        SEIR[1].push_back(SEIR[1][j - 1] - 1);
                        SEIR[2].push_back(SEIR[2][j - 1] + 1);
                        SEIR[3].push_back(SEIR[3][j - 1]);

                        // this is the rewrite of:
                        // sumsHiH = sumsHiH - (s * i)/nh + (s * (i + 1) )/nh

                        sumsHiH_nh = sumsHiH_nh + (double) s / selected_household.dimension;

                        states_to_household.add_household(selected_household, s, e - 1, i + 1);


                    } else {
                        // abbiamo solo uno spostamento da un exposed compartment al sucessivo
                        SEIR[0].push_back(SEIR[0][j - 1]);
                        SEIR[1].push_back(SEIR[1][j - 1]);
                        SEIR[2].push_back(SEIR[2][j - 1]);
                        SEIR[3].push_back(SEIR[3][j - 1]);

                        states_to_household.add_household(selected_household, s, e, i);


                    }

                    //--------------------------------------------------------------------------------------------------

                    goto skip;

                }
            }
        }
    }
    skip:;
}


void new_Recovered(std::vector<std::vector<int> > &SEIR,
                   state_to_household_map &states_to_household,
                   double &sumsHiH_nh, parameter &par, int &j) {


    //update households with susceptible based only on how many infected an house has
    //Note it may happen that thi is only a passage from an Infected compartment to another and no new Recovered are created
    //update also sumsHiH



    // generate a random number and decide which household will change


    std::uniform_int_distribution<int> random_int(1, SEIR[2][j - 1]);
    int random_in = random_int(par.generator);

    int cumulativeSum = 0;
    for (int e = 0; e <= par.nh_max; e++) {

        for (int i = 0; i <= par.nh_max - e; i++) {

            for (int s = 0; s <= par.nh_max - (e + i); s++) {
                cumulativeSum = cumulativeSum + ((int) states_to_household.matrix[s][e][i].size() * i);
                if (random_in <= cumulativeSum) {
                    //allora abbiamo estratto il numero (s,i,e)


                    //move infected

                    house selected_household = states_to_household.return_random_household_in_state(s, e, i, par);

                    // decide in which compartment happens the change based on how many infected each compartment has

                    std::uniform_int_distribution<int> random_int_compartment(1, i);
                    int random_compartment_i = random_int_compartment(par.generator);

                    int compartment = 1 + par.number_of_exposed_compartments;
                    int cum_sum = selected_household.state[compartment];
                    while (cum_sum < random_compartment_i) {
                        compartment++;
                        cum_sum = cum_sum + selected_household.state[compartment];
                    }


                    selected_household.state[compartment]--;

                    // note we make ++1 only if compartment is not the last one (i.e. there is no recovered)
                    //selected_household.state[compartment + 1]++;

                    if (compartment == par.number_of_exposed_compartments + par.number_of_infected_compartments) {
                        //abbiamo un nuovo Recovered
                        SEIR[0].push_back(SEIR[0][j - 1]);
                        SEIR[1].push_back(SEIR[1][j - 1]);
                        SEIR[2].push_back(SEIR[2][j - 1] - 1);
                        SEIR[3].push_back(SEIR[3][j - 1] + 1);

                        // this is the rewrite of:
                        // sumsHiH = sumsHiH - (s * i)/nh + (s * (i -1) )/nh

                        sumsHiH_nh = sumsHiH_nh - (double) s / selected_household.dimension;


                        states_to_household.add_household(selected_household, s, e, i - 1);


                    } else {
                        // abbiamo solo uno spostamento da un infected compartment al sucessivo
                        SEIR[0].push_back(SEIR[0][j - 1]);
                        SEIR[1].push_back(SEIR[1][j - 1]);
                        SEIR[2].push_back(SEIR[2][j - 1]);
                        SEIR[3].push_back(SEIR[3][j - 1]);


                        selected_household.state[compartment + 1]++;
                        states_to_household.add_household(selected_household, s, e, i);

                    }


                    goto skip;

                    //--------------------------------------------------------------------------------------------------





                    goto skip;

                }

            }
        }
    }
    skip:;

}


void new_vacinated(std::vector<std::vector<int> > &SEIR,
                   state_to_household_map &states_to_household,
                   int &sumsHiH_nh, parameter &par,
                   int &j) {
    //move some people randomly to recovered



    for (int l = 0; l < par.number_of_vaccinated; l++) {
        SEIR[0].push_back(SEIR[0][j] - 1);
        SEIR[1].push_back(SEIR[1][j]);
        SEIR[2].push_back(SEIR[2][j]);
        SEIR[3].push_back(SEIR[3][j] + 1);

        //update households with susceptible based only on how many susceptible an house has
        //update also sumsHiH
        // generate a random number and decide which household will change



        std::uniform_int_distribution<int> random_int(1, SEIR[0][j - 1]);
        int random_i = random_int(par.generator);

        //int size = par.nh_max;
        int cumulativeSum = 0;

        for (int e = 0; e <= par.nh_max; e++) {

            for (int i = 0; i <= par.nh_max - e; i++) {

                for (int s = 0; s <= par.nh_max - (e + i); s++) {
                    cumulativeSum = cumulativeSum + ((int) states_to_household.matrix[s][e][i].size() * s);
                    if (random_i <= cumulativeSum) {
                        //allora abbiamo estratto il numero (s,e,i)

                        //we have a new recovered

                        house selected_household = states_to_household.return_random_household_in_state(s, e, i, par);
                        selected_household.state[0]--;

                        states_to_household.add_household(selected_household, s - 1, e, i);

                        // this is the rewrite of:
                        // sumsHiH = sumsHiH - (s * i)/nh + ((s - 1) * i )/nh

                        sumsHiH_nh = sumsHiH_nh - (double) i / selected_household.dimension;



                        //--------------------------------------------------------------------------------------------------


                        goto skip;

                    }
                }
            }
        }
        skip:;
    }

}


void initializeSEIRandTemp(std::vector<std::vector<int> > &SEIR, std::vector<double> &temp, int &N) {
    SEIR[0][0] = N - 1;
    SEIR[1][0] = 0;
    SEIR[2][0] = 1;
    SEIR[3][0] = 0;
    temp.push_back(0);

}

double initialize_Households(parameter &par, state_to_household_map &states_to_household) {
    double sumsHiH_nh = 0;
    int total_population = 0;

    int n = 50;
    double p = (double) par.nh_mean / n;

    //std::binomial_distribution<int> bin_distribution(n, p);
    std::poisson_distribution<int> poiss_distribution(par.nh_mean);
    int nh = poiss_distribution(par.generator);

    //make sure that nh is in (0,par.nh_max]
    while (nh > par.nh_max || nh <= 0) {
        nh = poiss_distribution(par.generator);
    }

    //------------------------------------------------------------------------------------------------------------------

    nh = par.nh_mean;

    //------------------------------------------------------------------------------------------------------------------


    //initialize first household with one infected
    int s = nh - 1;
    int e = 0;
    int i = 1;
    house household_tmp(nh, s, e, i, par);
    states_to_household.add_household(household_tmp, s, e, i);
    sumsHiH_nh = sumsHiH_nh + ((double) s * i / nh);

    total_population = total_population + nh;


    //initialize remaining households with zero infected

    for (int z = 1; z < par.number_of_Households; z++) {

        nh = poiss_distribution(par.generator);
        //make sure that nh is in (0,par.nh_max]
        while (nh > par.nh_max || nh <= 0) {
            nh = poiss_distribution(par.generator);
        }

        //------------------------------------------------------------------------------------------------------------------

        nh = par.nh_mean;

        //------------------------------------------------------------------------------------------------------------------

        s = nh;
        e = 0;
        i = 0;

        house household_tmp(nh, s, e, i, par);
        states_to_household.add_household(household_tmp, s, e, i);
        sumsHiH_nh = sumsHiH_nh + (s * i / nh);

        total_population = total_population + nh;
    }

    par.N = total_population;
    return sumsHiH_nh;
}


void read_Parameters_From_File(std::string inputpath, parameter &parameters) {

    std::string line;
    std::ifstream infile(inputpath);
    if (infile.is_open()) {

        //number of steps
        getline(infile, line, ':');
        getline(infile, line);
        parameters.nSteps = std::stoi(line);

        //number of Households
        getline(infile, line, ':');
        getline(infile, line);
        parameters.number_of_Households = std::stoi(line);

        // maximum number of people in one Household
        getline(infile, line, ':');
        getline(infile, line);
        parameters.nh_max = std::stoi(line);

        // mean number of people in one Household
        getline(infile, line, ':');
        getline(infile, line);
        parameters.nh_mean = std::stoi(line);

        //beta1 is the initial beta
        getline(infile, line, ':');
        getline(infile, line);
        parameters.beta1 = std::stod(line);

        //beta2 is the second beta after the start of the lockdown
        getline(infile, line, ':');
        getline(infile, line);
        parameters.beta2 = std::stod(line);

        //beta3 is the third beta after the end of the lockdown
        getline(infile, line, ':');
        getline(infile, line);
        parameters.beta3 = std::stod(line);

        //
        getline(infile, line, ':');
        getline(infile, line);
        parameters.threshold_above_which_one_to_two = std::stod(line);

        // the
        getline(infile, line, ':');
        getline(infile, line);
        parameters.threshold_under_which_two_to_one = std::stod(line);

        //betaH
        getline(infile, line, ':');
        getline(infile, line);
        parameters.betaH = std::stod(line);


        //ny
        getline(infile, line, ':');
        getline(infile, line);
        parameters.ny = std::stod(line);


        //gamma
        getline(infile, line, ':');
        getline(infile, line);
        parameters.gamma = std::stod(line);

        //number_of_infected_compartments
        getline(infile, line, ':');
        getline(infile, line);
        parameters.number_of_infected_compartments = std::stod(line);


        //number_of_exposed_compartments
        getline(infile, line, ':');
        getline(infile, line);
        parameters.number_of_exposed_compartments = std::stod(line);


        //time at which we activate the lockdown
        getline(infile, line, ':');
        getline(infile, line);
        parameters.time_activate_lockdown = std::stod(line);


        //time at which we end the lockdown
        getline(infile, line, ':');
        getline(infile, line);
        parameters.time_end_lockdown = std::stod(line);

        infile.close();
    } else std::cout << "Unable to open file";
}

void write_the_csv_file(std::string outputpath, std::vector<std::vector<int> > &SEIR, std::vector<double> &temp) {
    std::ofstream outfile(outputpath);
    if (!outfile.is_open()) {
        std::cout << "Unable to open file";
    } else {
        double print_time = 0;
        for (int i = 0; i < temp.size(); i++) {
            //write only every one unit of time
            if (temp[i] >= print_time) {
                print_time = print_time + 1;
                outfile << SEIR[0][i] << ",\t" << SEIR[1][i] << ",\t" << SEIR[2][i] << ",\t" << SEIR[3][i] << ",\t"
                        << temp[i] << "\n";
            }

        }
        outfile.close();
    }
}


void write_lock_down_files(std::string outputpath, std::vector<double> &time_lockdown) {
    std::ofstream outfile(outputpath);
    if (!outfile.is_open()) {
        std::cout << "Unable to open file";
    } else {
        for (int i = 0; i < time_lockdown.size(); i++) {
            if (i % 2 == 1) {
                outfile << '\n';
            }
            outfile << time_lockdown[i] << " ";
        }
        outfile.close();
    }

}


double generateUnif_from_zeroExcluded_to(double to, parameter &par) {
    //std::default_random_engine generator(time(0));
    std::uniform_real_distribution<double> uniform_distribution(0.0, to);

    double randomUnif = uniform_distribution(par.generator);
    while (randomUnif == 0) {
        uniform_distribution.reset();
        randomUnif = uniform_distribution(par.generator);
    }
    return randomUnif;
}


