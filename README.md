# Gillespie_for_Households_gamma_distribution
Gillespie algorithm with Gamma-Distributed times

An implementation of Gillespie algorithm for epidemics.

In the input file it is possible to specify the parameters of the epidemy:

NumberOfSteps: -
number_of_Households: -
maximum_number_of_people_in_one_Household: -
mean_number_of_people_in_one_household: -
Beta_1: Gobal transmission rate of an infect
Beta_2: Gobal transmission rate of an infect during lockdown
threshold_above_which_we_change_from_beta1_to_beta2: percentage of exposed above which we activate the lockdown policy
threshold_under_which_we_change_from_beta2_to_beta1: percentage of exposed under which we ease the lockdown policy
BetaH: Trsmisison rate within the household
Ny: Exposed time parameter
Gamma: Recover rate parameter
number_of_infected_compartments: to be used to transform the infected time from exponential to a gamma
number_of_exposed_compartments: to be used to transform the exposed time from exponential to a gamma

The household size distribution follows a binomial distribution centered around the mean an then truncated at "maximum_number_of_people_in_one_Household".

The output is a csv file of five columnse S, E, I, R, time:
-S suceptibles at time t
-E exposed at time t
-I infected at time t
-R recoveres at time t
-time t

In the case of lockdown it is created also a .txt file containig the times at which the lock down poicy is activated/eased.
