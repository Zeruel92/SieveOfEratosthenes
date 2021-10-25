cmake_minimum_required(VERSION 3.20)
project(Eratosthenes C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-O3 -Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG  "-g -O2 -Wall -Wextra")
set(CMAKE_C_COMPILER mpiicc)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(sieve sieve.c)
target_link_libraries(sieve ${MPI_LIBRARIES})
add_executable(sieve_noeven sieve_noeven.c)
target_link_libraries(sieve_noeven ${MPI_LIBRARIES})
