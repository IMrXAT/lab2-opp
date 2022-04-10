cmake_minimum_required(VERSION 3.21)

project(lab2 C)

set(CMAKE_C_STANDARD 99)
add_executable(lab2 main.c)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 ")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fopenmp")
target_link_libraries(lab2 m)
