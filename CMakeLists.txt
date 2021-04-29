cmake_minimum_required(VERSION 3.17)
project(Gillespie_for_Households_gamma_distribution)

set(CMAKE_CXX_STANDARD 20)

add_executable(Gillespie_for_Households_gamma_distribution main.cpp myFunctions.h myFunctions.cpp GillespieForHouseholds.h GillespieForHouseholds.cpp myClass.h myClass.cpp)