//
// Created by popcorn on 18/04/2021.
//


#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <map>
#include<tuple>
#include <math.h>
#include "myFunctions.h"


void new_Exposed_outside_the_household(std::vector<std::vector<int> > &SEIR,
                                       std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                                       int &sumsHiH,
                                       std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                                       std::vector<std::vector<int> > &households, int number_of_infected_compartments,
                                       int number_of_exposed_compartments, int &j) {
    SEIR[0].push_back(SEIR[0][j - 1] - 1);
    SEIR[1].push_back(SEIR[1][j - 1] + 1);
    SEIR[2].push_back(SEIR[2][j - 1]);
    SEIR[3].push_back(SEIR[3][j - 1]);

    //update households with susceptible based only on how many susceptible an house has
    //update also sumsHiH
    // generate a random number and decide which household will change

    double randomUnif = generateUnif_from_zeroExcluded_to(SEIR[0][j - 1]);
    int size = household_with_Susceptible_Infected_Exposed[0].size();
    int cumulativeSum = 0;

    for (int e = 0; e < size; e++) {

        for (int i = 0; i < size - e; i++) {

            for (int s = 0; s < size - (e + i); s++) {
                cumulativeSum = cumulativeSum + (household_with_Susceptible_Infected_Exposed[s][i][e] * s);
                if (randomUnif <= cumulativeSum) {
                    //allora abbiamo estratto il numero (s,i,e)
                    household_with_Susceptible_Infected_Exposed[s][i][e]--;
                    household_with_Susceptible_Infected_Exposed[s - 1][i][e + 1]++;

                    // this is the rewrite of:
                    // sumsHiH = sumsHiH - (s * i) + ((s - 1) * i )
                    sumsHiH = sumsHiH - i;

                    //--------------------------------------------------------------------------------------------------


                    //we have a new exposed
                    std::tuple<int, int, int> extracted_tuple = std::make_tuple(s, e, i);
                    std::vector<int> possible_households = states_to_households[extracted_tuple];
                    int random_int = rand() % possible_households.size();
                    int choosed_household = possible_households[random_int];
                    households[choosed_household][0]--;
                    households[choosed_household][1]++;

                    // we move the chosen_household
                    //check if the key (s, e - 1, i + 1) exist
                    std::map<std::tuple<int, int, int>, std::vector<int> >::iterator it;
                    it = states_to_households.find(std::make_tuple(s - 1, e + 1, i));
                    if (it != states_to_households.end())
                        states_to_households[std::make_tuple(s - 1, e + 1, i)].push_back(
                                choosed_household);
                    else {
                        states_to_households[std::make_tuple(s - 1, e + 1, i)] = std::vector<int>(0);
                        states_to_households[std::make_tuple(s - 1, e + 1, i)].push_back(
                                choosed_household);



                    }
                    states_to_households[extracted_tuple].erase(
                            states_to_households[extracted_tuple].begin() + random_int);

                    //--------------------------------------------------------------------------------------------------


                    goto skip;

                }
            }
        }
    }
    skip:;
}


void new_exposed_inside_the_household(std::vector<std::vector<int>> &SEIR,
                                      std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                                      int &sumsHiH,
                                      std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                                      std::vector<std::vector<int> > &households, int number_of_infected_compartments,
                                      int number_of_exposed_compartments, int &j) {
    SEIR[0].push_back(SEIR[0][j - 1] - 1);
    SEIR[1].push_back(SEIR[1][j - 1] + 1);
    SEIR[2].push_back(SEIR[2][j - 1]);
    SEIR[3].push_back(SEIR[3][j - 1]);

    //update households with susceptible based on how many susceptible and infected an house has
    //update also sumsHiH
    // generate a random number and decide which household will change
    double randomUnif = generateUnif_from_zeroExcluded_to(sumsHiH);
    int size = household_with_Susceptible_Infected_Exposed[0].size();
    int cumulativeSum = 0;

    for (int e = 0; e < size; e++) {

        for (int i = 0; i < size - e; i++) {

            for (int s = 0; s < size - (e + i); s++) {
                cumulativeSum = cumulativeSum + (household_with_Susceptible_Infected_Exposed[s][i][e] * s * i);
                if (randomUnif <= cumulativeSum) {

                    //allora abbiamo estratto il numero (s,i,e)
                    household_with_Susceptible_Infected_Exposed[s][i][e]--;
                    household_with_Susceptible_Infected_Exposed[s - 1][i][e + 1]++;

                    // this is the rewrite of:
                    // sumsHiH = sumsHiH - (s * i) + ((s - 1) * i )
                    sumsHiH = sumsHiH - i;


                    //--------------------------------------------------------------------------------------------------


                    //we have a new exposed
                    std::tuple<int, int, int> extracted_tuple = std::make_tuple(s, e, i);
                    std::vector<int> possible_households = states_to_households[extracted_tuple];
                    int random_int = rand() % possible_households.size();
                    int choosed_household = possible_households[random_int];
                    households[choosed_household][0]--;
                    households[choosed_household][1]++;

                    // we move the chosen_household
                    //check if the key (s, e - 1, i + 1) exist
                    std::map<std::tuple<int, int, int>, std::vector<int> >::iterator it;
                    it = states_to_households.find(std::make_tuple(s - 1, e + 1, i));
                    if (it != states_to_households.end())
                        states_to_households[std::make_tuple(s - 1, e + 1, i)].push_back(
                                choosed_household);
                    else {
                        states_to_households[std::make_tuple(s - 1, e + 1, i)] = std::vector<int>(0);
                        states_to_households[std::make_tuple(s - 1, e + 1, i)].push_back(
                                choosed_household);





                    }
                    states_to_households[extracted_tuple].erase(
                            states_to_households[extracted_tuple].begin() + random_int);


                    //--------------------------------------------------------------------------------------------------

                    goto skip;


                }

            }
        }
    }
    skip:;

}

void new_Infected(std::vector<std::vector<int> > &SEIR,
                  std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                  int &sumsHiH, std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                  std::vector<std::vector<int> > &households, int number_of_infected_compartments,
                  int number_of_exposed_compartments, int &j) {


    //update households with susceptible based only on how many exposed an house has
    //update also sumsHiH
    // generate a random number and decide which household will change

    double randomUnif = generateUnif_from_zeroExcluded_to(SEIR[1][j - 1]);

    int size = household_with_Susceptible_Infected_Exposed[0].size();
    int cumulativeSum = 0;

    for (int e = 0; e < size; e++) {

        for (int i = 0; i < size - e; i++) {

            for (int s = 0; s < size - (e + i); s++) {
                cumulativeSum = cumulativeSum + (household_with_Susceptible_Infected_Exposed[s][i][e] * e);
                if (randomUnif <= cumulativeSum) {
                    //allora abbiamo estratto il numero (s, e, i)

                    //move exposed
                    std::tuple<int, int, int> extracted_tuple = std::make_tuple(s, e, i);
                    std::vector<int> possible_households = states_to_households[extracted_tuple];

                    // we know there are e*household_with_Susceptible_Infected_Exposed[s][e][i] exposed individuals
                    int random_int = rand() % possible_households.size();
                    int choosed_household = possible_households[random_int];


                    int cumsum = 0;
                    double random_compartment = generateUnif_from_zeroExcluded_to(e);
                    for (int z = 1; z < 1 + number_of_exposed_compartments; z++) {
                        cumsum = cumsum + households[choosed_household][z];
                        if (cumsum >= random_compartment) {
                            //the compartment in which the change will happen is z
                            households[choosed_household][z]--;
                            households[choosed_household][z + 1]++;
                            if (z < number_of_exposed_compartments) {
                                // we just have a move from one exposed compartment to another
                                SEIR[0].push_back(SEIR[0][j - 1]);
                                SEIR[1].push_back(SEIR[1][j - 1]);
                                SEIR[2].push_back(SEIR[2][j - 1]);
                                SEIR[3].push_back(SEIR[3][j - 1]);

                            } else {
                                //we have a new infected

                                SEIR[0].push_back(SEIR[0][j - 1]);
                                SEIR[1].push_back(SEIR[1][j - 1] - 1);
                                SEIR[2].push_back(SEIR[2][j - 1] + 1);
                                SEIR[3].push_back(SEIR[3][j - 1]);

                                household_with_Susceptible_Infected_Exposed[s][i][e]--;
                                household_with_Susceptible_Infected_Exposed[s][i + 1][e - 1]++;

                                // this is the rewrite of:
                                // sumsHiH = sumsHiH - (s * i) + (s * (i+1) )
                                sumsHiH = sumsHiH + s;

                                // we move the chosen_household
                                //check if the key (s, e - 1, i + 1) exist
                                std::map<std::tuple<int, int, int>, std::vector<int> >::iterator it;
                                it = states_to_households.find(std::make_tuple(s, e - 1, i + 1));
                                if (it != states_to_households.end())
                                    states_to_households[std::make_tuple(s, e - 1, i + 1)].push_back(
                                            choosed_household);
                                else {
                                    states_to_households[std::make_tuple(s, e - 1, i + 1)] = std::vector<int>(0);
                                    states_to_households[std::make_tuple(s, e - 1, i + 1)].push_back(
                                            choosed_household);


                                }
                                states_to_households[extracted_tuple].erase(
                                        states_to_households[extracted_tuple].begin() + random_int);
                            }
                            break;
                        }

                    }


                    goto skip;

                }
            }
        }
    }
    skip:;
}


void new_Recovered(std::vector<std::vector<int> > &SEIR,
                   std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
                   int &sumsHiH, std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households,
                   std::vector<std::vector<int> > &households, int number_of_infected_compartments,
                   int number_of_exposed_compartments, int &j) {
    SEIR[0].push_back(SEIR[0][j - 1]);
    SEIR[1].push_back(SEIR[1][j - 1]);
    SEIR[2].push_back(SEIR[2][j - 1] - 1);
    SEIR[3].push_back(SEIR[3][j - 1] + 1);

    //update households with susceptible based only on how many infected an house has
    //update also sumsHiH
    // generate a random number and decide which household will change

    double randomUnif = generateUnif_from_zeroExcluded_to(SEIR[2][j - 1]);
    int size = household_with_Susceptible_Infected_Exposed[0].size();
    int cumulativeSum = 0;

    for (int e = 0; e < size; e++) {

        for (int i = 0; i < size - e; i++) {

            for (int s = 0; s < size - (e + i); s++) {
                cumulativeSum = cumulativeSum + (household_with_Susceptible_Infected_Exposed[s][i][e] * i);
                if (randomUnif <= cumulativeSum) {
                    //allora abbiamo estratto il numero (s,i,e)
                    household_with_Susceptible_Infected_Exposed[s][i][e]--;
                    household_with_Susceptible_Infected_Exposed[s][i - 1][e]++;

                    // this is the rewrite of:
                    // sumsHiH = sumsHiH - (prec[0] * prec[1]) + (prec[0] * (prec[1]-1) )
                    sumsHiH = sumsHiH - s;


                    //--------------------------------------------------------------------------------------------------


                    //we have a new exposed
                    std::tuple<int, int, int> extracted_tuple = std::make_tuple(s, e, i);
                    std::vector<int> possible_households = states_to_households[extracted_tuple];
                    int random_int = rand() % possible_households.size();
                    int choosed_household = possible_households[random_int];
                    households[choosed_household][number_of_exposed_compartments + 1]--;
                    //households[choosed_household][-1]++;


                    // we move the chosen_household
                    //check if the key (s, e , i - 1) exist
                    std::map<std::tuple<int, int, int>, std::vector<int> >::iterator it;
                    it = states_to_households.find(std::make_tuple(s, e, i - 1));
                    if (it != states_to_households.end())
                        states_to_households[std::make_tuple(s, e, i - 1)].push_back(
                                choosed_household);
                    else {
                        states_to_households[std::make_tuple(s, e, i - 1)] = std::vector<int>(0);
                        states_to_households[std::make_tuple(s, e, i - 1)].push_back(
                                choosed_household);
                    }
                    states_to_households[extracted_tuple].erase(
                            states_to_households[extracted_tuple].begin() + random_int);

                    //--------------------------------------------------------------------------------------------------





                    goto skip;

                }

            }
        }
    }
    skip:;

}


void initializeSEIRandTemp(std::vector<std::vector<int> > &SEIR, std::vector<double> &temp, int &N) {
    SEIR[0][0] = N - 1;
    SEIR[1][0] = 0;
    SEIR[2][0] = 1;
    SEIR[3][0] = 0;
    temp.push_back(0);

}

void initialize_Households(std::vector<std::vector<int> > &households, int nh, int number_of_exposed_compartments,
                           std::map<std::tuple<int, int, int>, std::vector<int> > &states_to_households) {
    //first household will have one infected,the others none
    households[0][0] = nh - 1;
    states_to_households[std::make_tuple(nh - 1, 0, 1)] = std::vector<int>(1, 0);

    states_to_households[std::make_tuple(nh, 0, 0)] = std::vector<int>(0);
    //alternative way
    //states_to_households.insert(std::pair<std::tuple<int, int, int>,std::vector<int>>(std::make_tuple(nh, 0, 0), std::vector<int>(0)));


    households[0][number_of_exposed_compartments + 1] = 1;
    for (int i = 1; i < households.size(); i++) {
        households[i][0] = nh;
        states_to_households[std::make_tuple(nh, 0, 0)].push_back(i);
    }


}


void read_Parameters_From_File(std::string inputpath, int &nSteps, int &number_of_Households,
                               int &number_of_people_in_one_Household, double &beta1, double &beta2,
                               double &threshold_above_which_one_to_two, double &threshold_under_which_two_to_one,
                               double &betaH, double &ny, double &gamma, int &number_of_infected_compartments,
                               int &number_of_exposed_compartments) {

    std::string line;
    std::ifstream infile(inputpath);
    if (infile.is_open()) {

        //number of steps
        getline(infile, line, ':');
        getline(infile, line);
        nSteps = std::stoi(line);

        //number of Households
        getline(infile, line, ':');
        getline(infile, line);
        number_of_Households = std::stoi(line);

        // number of people in one Household
        getline(infile, line, ':');
        getline(infile, line);
        number_of_people_in_one_Household = std::stoi(line);


        //beta1 is the initial beta
        getline(infile, line, ':');
        getline(infile, line);
        beta1 = std::stod(line);

        //beta2 is the second beta after the start of the lockdown
        getline(infile, line, ':');
        getline(infile, line);
        beta2 = std::stod(line);

        //
        getline(infile, line, ':');
        getline(infile, line);
        threshold_above_which_one_to_two = std::stod(line);

        // the
        getline(infile, line, ':');
        getline(infile, line);
        threshold_under_which_two_to_one = std::stod(line);

        //betaH
        getline(infile, line, ':');
        getline(infile, line);
        betaH = std::stod(line);


        //ny
        getline(infile, line, ':');
        getline(infile, line);
        ny = std::stod(line);


        //gamma
        getline(infile, line, ':');
        getline(infile, line);
        gamma = std::stod(line);

        //number_of_infected_compartments
        getline(infile, line, ':');
        getline(infile, line);
        number_of_infected_compartments = std::stod(line);


        //number_of_exposed_compartments
        getline(infile, line, ':');
        getline(infile, line);
        number_of_exposed_compartments = std::stod(line);

        infile.close();
    } else std::cout << "Unable to open file";
}

void write_the_csv_file(std::string outputpath, std::vector<std::vector<int> > &SEIR, std::vector<double> &temp) {
    std::ofstream outfile(outputpath);
    if (!outfile.is_open()) {
        std::cout << "Unable to open file";
    } else {
        for (int i = 0; i < temp.size(); i++) {
            //write only every 10
            if (i % 1 == 0) {
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
        for (int i = 1; i < time_lockdown.size(); i++) {
            if (i % 2 == 1) {
                outfile << '\n';
            }
            outfile << time_lockdown[i] << " ";
        }
        outfile.close();
    }

}


void initialize_household_with_Susceptible_Infected_Exposed(
        std::vector<std::vector<std::vector<int>>> &household_with_Susceptible_Infected_Exposed,
        int number_of_Households, int number_of_people_in_one_Household) {
    household_with_Susceptible_Infected_Exposed[number_of_people_in_one_Household][0][0] = number_of_Households - 1;
    household_with_Susceptible_Infected_Exposed[number_of_people_in_one_Household - 1][1][0] = 1;
}


double generateUnif_from_zeroExcluded_to(double to) {
    std::default_random_engine generator(0);
    std::uniform_real_distribution<double> uniform_distribution(0.0, to);

    double randomUnif = uniform_distribution(generator);
    while (randomUnif == 0) {
        uniform_distribution.reset();
        randomUnif = uniform_distribution(generator);
    }
    return randomUnif;
}
