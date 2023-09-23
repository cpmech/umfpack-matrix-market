cmake_minimum_required(VERSION 3.22)
project(Project)

add_definitions(-Wall)

include(scripts/FindUMFPACK.cmake)
if(UMFPACK_FOUND)
    include_directories(${UMFPACK_INCLUDE})
else()
    message(FATAL_ERROR "cannot compile the code without UMFPACK")
endif()

add_executable(solve_matrix_market solve_matrix_market.cpp)
target_link_libraries(solve_matrix_market PUBLIC ${UMFPACK_LIBS})

set(SAMPLES
    bfwb62.mtx
    pressurized-cylinder-linear-elastic-symmetric-1.mtx
    pressurized-cylinder-linear-elastic-symmetric-2.mtx
)

foreach(sample ${SAMPLES})
    configure_file(samples/${sample} ${sample} @ONLY)
endforeach()
