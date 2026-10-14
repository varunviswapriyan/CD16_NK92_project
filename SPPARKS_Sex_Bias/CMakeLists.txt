cmake_minimum_required(VERSION 3.22)
project(SPPARKS)

set(CMAKE_CXX_STANDARD 14)

add_executable(SPPARKS main.cpp
        app_lattice.cpp
        app_rxn_diff_custom.cpp
        app_rxn_diff.cpp)
