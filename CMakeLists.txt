cmake_minimum_required(VERSION 3.25)
project(Project)

include(scripts/FindUMFPACK.cmake)
include_directories(${UMFPACK_INCLUDE_DIRS})

add_executable(solve_matrix_market solve_matrix_market.cpp)
target_link_libraries(solve_matrix_market  PUBLIC ${UMFPACK_LIBS})
