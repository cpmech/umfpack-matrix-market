set(UMFPACK_INCLUDE_SEARCH_PATH
  /usr/include/suitesparse
  /usr/local/include/umfpack)

set(UMFPACK_LIBRARY_SEARCH_PATH
  /usr/lib/x86_64-linux-gnu
  /usr/local/lib/umfpack
  /usr/lib)

find_path(UMFPACK_AMD_H      amd.h      ${UMFPACK_INCLUDE_SEARCH_PATH})
find_path(UMFPACK_UMFPACK_H  umfpack.h  ${UMFPACK_INCLUDE_SEARCH_PATH})

find_library(UMFPACK_AMD     NAMES amd     PATHS ${UMFPACK_LIBRARY_SEARCH_PATH})
find_library(UMFPACK_UMFPACK NAMES umfpack PATHS ${UMFPACK_LIBRARY_SEARCH_PATH})

set(UMFPACK_FOUND 1)

foreach(var UMFPACK_UMFPACK_H UMFPACK_AMD_H UMFPACK_UMFPACK UMFPACK_AMD)
    if(NOT ${var})
        set(UMFPACK_FOUND 0)
    endif()
endforeach()

message(STATUS "UMFPACK_FOUND = ${UMFPACK_FOUND}")

if(UMFPACK_FOUND)
    set(UMFPACK_INCLUDE ${UMFPACK_AMD_H} ${UMFPACK_UMFPACK_H})
    set(UMFPACK_LIBS    ${UMFPACK_UMFPACK} ${UMFPACK_AMD})
endif()

message(STATUS "UMFPACK_INCLUDE_DIRS = ${UMFPACK_INCLUDE_DIRS}")
message(STATUS "UMFPACK_LIBRARIES = ${UMFPACK_LIBRARIES}")

if(UMFPACK_FOUND)
    message(STATUS "UMFPACK: Found Successfully!!!")
endif()
