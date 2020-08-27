cmake_minimum_required(VERSION 3.10)
#find_package(Python REQUIRED)
#find_package(HDF5 REQUIRED)
#
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(SCS_LIB NAMES scsdir libscsdir
  REQUIRED)
find_path(SCS_PATH NAMES scsdir libscsdir scs REQUIRED)

message(STATUS "${SCS_LIB}")
message(STATUS "${SCS_PATH}")
if(NOT SCS_LIB)
  message(FATAL_ERROR "SCSDIR library not found")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CARGO_CMD cargo build)
  set(TARGET_DIR "debug")
else ()
  set(CARGO_CMD cargo build --release)
  set(TARGET_DIR "release")
endif ()

set(SENT_LIB "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_DIR}/libs1_noisefloor.so")
set(SENT_BIN "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_DIR}/denoise_s1")
set(SENT_RLIB "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_DIR}/libs1_noisefloor_engine.rlib")

set(C_SOLVE_DUR "${CMAKE_CURRENT_BINARY_DIR}/engine/src/c_solve")

add_library(lp_solve STATIC
  ${CMAKE_CURRENT_BINARY_DIR}/engine/src/c_solve/scs_solve.c
  ${CMAKE_CURRENT_BINARY_DIR}/engine/src/c_solve/scs_solve.h)
set_target_properties(lp_solve PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/engine/src/c_solve/)

add_custom_target(sentlib ALL
    COMMENT "Compiling Sentinel library"
    COMMAND CARGO_TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR} ${CARGO_CMD}
    # copy result from library and bin.
    COMMAND cp ${SENT_LIB} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND cp ${SENT_RLIB} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND cp ${SENT_BIN} ${CMAKE_CURRENT_BINARY_DIR}
    # link library so it is accessible from python
    COMMAND ln -s "${CMAKE_CURRENT_BINARY_DIR}/libs1_noisefloor.so" "${CMAKE_CURRENT_BINARY_DIR}/s1_noisefloor.so"
    
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(sentlib PROPERTIES LOCATION ${CMAKE_CURRENT_BINARY_DIR})
