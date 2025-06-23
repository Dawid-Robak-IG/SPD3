cmake_minimum_required(VERSION 3.30)
project(spd3)

set(CMAKE_CXX_STANDARD 20)

add_executable(spd3 main.cpp
        Task.cpp
        Task.h
        Problem.cpp
        Problem.h)
