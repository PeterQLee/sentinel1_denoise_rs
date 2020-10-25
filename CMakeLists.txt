
find_library(SCS_LIB NAMES scsdir  REQUIRED)
find_path(SCS_PATH NAMES scs REQUIRED)

message(STATUS "${SCS_LIB}")
message(STATUS "${SCS_PATH}")

message( STATUS "Found SCS_LIB=${SCS_LIB} and SCS_PATH=${SCS_PATH}")
if(NOT SCS_LIB)
  message(FATAL_ERROR "SCS library not found")
endif()
if(NOT SCS_LIB)
  message(FATAL_ERROR "SCS include files not found")
endif()

get_filename_component(SCS_LIB_DIR ${SCS_LIB} DIRECTORY)
message(STATUS "${SCS_LIB_DIR}")


set(SENT_LIB "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_DIR}/libs1_noisefloor${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(SENT_BIN "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_DIR}/denoise_s1${CMAKE_EXECUTABLE_SUFFIX}")
set(SENT_RLIB "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_DIR}/libs1_noisefloor_engine.rlib")
set(SENT_INCLUDE_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/interface.h")

set(C_SOLVE_DIR "${CMAKE_CURRENT_BINARY_DIR}/engine/src/c_solve")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# supporting lp solve c program
add_library(lp_solve STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/engine/src/c_solve/scs_solve.c
  ${CMAKE_CURRENT_SOURCE_DIR}/engine/src/c_solve/scs_solve.h)
include_directories(${SCS_PATH})
set_target_properties(lp_solve PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/engine/src/c_solve/)

#
# message( STATUS "Found SCS_LIB_DIR=${SCS_LIB_DIR}")
# add_custom_target(sentlib ALL

#     COMMENT "Compiling Sentinel library"
#     COMMAND ${CMAKE_COMMAND} -E env SCS_LIB=${SCS_LIB_DIR} env CARGO_TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR} ${CARGO_CMD}
#     # copy result from library and bin.
#     COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SENT_LIB} ${CMAKE_CURRENT_BINARY_DIR}
#     COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SENT_RLIB} ${CMAKE_CURRENT_BINARY_DIR}
#     COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SENT_BIN} ${CMAKE_CURRENT_BINARY_DIR}
    
#     # link library so it is accessible from python
#     COMMAND ${CMAKE_COMMAND} -E create_symlink "${SENT_LIB}" "${CMAKE_CURRENT_BINARY_DIR}/s1_noisefloor${CMAKE_SHARED_LIBRARY_SUFFIX}"

    
#     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# add_dependencies(sentlib lp_solve)

# execute_process(
#   COMMAND "${PYTHON_EXECUTABLE}" -c "
# from distutils import sysconfig as sc
# print(sc.get_python_lib(plat_specific=True))"
#   OUTPUT_VARIABLE PYTHON_SITE
#   OUTPUT_STRIP_TRAILING_WHITESPACE)


# message (STATUS "python sitelib ${PYTHON_SITE}")

# set(SENT_LIB_INSTALL "${CMAKE_INSTALL_LIBDIR}/libs1_noisefloor${CMAKE_SHARED_LIBRARY_SUFFIX}")


# # Installation entries
# include(GNUInstallDirs)
# install(FILES ${SENT_INCLUDE_HEADER}
#   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/s1_noisefloor/
#     PERMISSIONS OWNER_WRITE OWNER_READ
#     GROUP_READ
#     WORLD_READ
#   )

# install(FILES ${SENT_LIB}
#   DESTINATION ${CMAKE_INSTALL_LIBDIR}
#   PERMISSIONS OWNER_EXECUTE OWNER_READ
#   GROUP_EXECUTE GROUP_READ
#   WORLD_READ WORLD_EXECUTE
#   )


# install(FILES ${SENT_BIN}
#   DESTINATION ${CMAKE_INSTALL_BINDIR}
#   PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
#   GROUP_EXECUTE GROUP_READ
#   WORLD_READ WORLD_EXECUTE
#   )
# install(CODE "execute_process( \
#     COMMAND ${CMAKE_COMMAND} -E create_symlink \
#     ${SENT_LIB_INSTALL} \
#     ${PYTHON_SITE}/s1_noisefloor   \
#     )")
