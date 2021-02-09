cmake_minimum_required( VERSION 2.8 )
project(decifer)

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

set( mergestatetrees_src/cpp
  mergestatetreesmain.cpp
  statetree.cpp
  stategraph.cpp
  utils.cpp
  readmatrix.cpp
  basematrix.cpp
)

set( mergestatetrees_hdr
  statetree.h
  stategraph.h
  utils.h
  readmatrix.h
  basematrix.h
)

set( generatestatetrees_src/cpp
  generatestatetreesmain.cpp
  statetree.cpp
  stategraph.cpp
  utils.cpp
  readmatrix.cpp
  basematrix.cpp
  solver.cpp
)

set( generatestatetrees_hdr
  statetree.h
  stategraph.h
  utils.h
  readmatrix.h
  basematrix.h
  solver.h
)

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#set( Boost_USE_STATIC_LIBS ON )

find_package( Boost REQUIRED COMPONENTS thread system filesystem )

MESSAGE( STATUS "found boost:"
  "\nINCLUDE: ${Boost_INCLUDE_DIRS}"
  "\nLIB: ${Boost_LIBRARIES}"
)

include_directories( "${LIBLEMON_ROOT}/include" ${Boost_INCLUDE_DIRS} )
link_directories( ${Boost_LIBRARY_DIRS} "${LIBLEMON_ROOT}/lib" ${GUROBI_LIBRARIES} ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} )

set( CommonLibs
  m
  emon
  pthread
  ${Boost_LIBRARIES}
)

set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" )

add_executable( mergestatetrees ${mergestatetrees_src/cpp} ${mergestatetrees_hdr} )
target_link_libraries( mergestatetrees ${CommonLibs} )

add_executable( generatestatetrees ${generatestatetrees_src/cpp} ${generatestatetrees_hdr} )
target_link_libraries( generatestatetrees ${CommonLibs} )

