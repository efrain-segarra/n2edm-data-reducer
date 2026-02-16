cmake_minimum_required(VERSION 3.11)

project(data_reduction)
# set Makefile verbosity on
#set( CMAKE_VERBOSE_MAKEFILE on )

# configure compiler version
#set(CMAKE_CXX_COMPILER "g++-7.3")
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 7.3)
        message("Found g++ >= 7.3. Continuing")
else()
        message( FATAL_ERROR "Need g++ >= 7.3. Exiting" )
endif()
# configure compiler flags and echo them during cmake
set(CMAKE_CXX_FLAGS "-ggdb -O0 -Wall -std=c++20  -fdiagnostics-color=always")





find_package(ROOT REQUIRED)

# Find n2dataread headers and ROOT
include_directories( 
	${CMAKE_SOURCE_DIR}/../n2dataread/src 
	${ROOT_INCLUDE_DIRS}
	)

add_executable( data_reducer data_reducer.cpp )
target_link_libraries( data_reducer 
	LINK_PUBLIC 
		${CMAKE_SOURCE_DIR}/../n2dataread/install/libN2readData_dyn.dylib 
		${ROOT_LIBRARIES} 
	)
