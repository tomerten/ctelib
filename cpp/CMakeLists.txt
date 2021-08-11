cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

SET(GCC_OPENCL_COMPILE_FLAGS "-I/usr/include -lOpenCL")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_OPENCL_COMPILE_FLAGS}")

project(CTELib VERSION 0.0.0)

# include dir
include_directories(/usr/local/include)

# Src
AUX_SOURCE_DIRECTORY(src SRC_FILES)

# Headers
set(PROJECT_SOURCE_DIR "src")
set(PROJECT_INCLUDE_DIR "include/cte_bits")

# Source files
set(SOURCE_FILES
${PROJECT_INCLUDE_DIR}/random.hpp
${PROJECT_INCLUDE_DIR}/radiation.hpp
${PROJECT_INCLUDE_DIR}/synch.hpp
${PROJECT_INCLUDE_DIR}/hamiltonian.hpp
${PROJECT_INCLUDE_DIR}/distributions.hpp
${PROJECT_INCLUDE_DIR}/physics.hpp
${PROJECT_INCLUDE_DIR}/utils.hpp
${PROJECT_INCLUDE_DIR}/bunch.hpp
${PROJECT_SOURCE_DIR}/random.cpp
${PROJECT_SOURCE_DIR}/radiation.cpp
${PROJECT_SOURCE_DIR}/hamiltonian.cpp
${PROJECT_SOURCE_DIR}/distributions.cpp
${PROJECT_SOURCE_DIR}/synch.cpp
${PROJECT_SOURCE_DIR}/physics.cpp
${PROJECT_SOURCE_DIR}/utils.cpp
${PROJECT_SOURCE_DIR}/bunch.cpp
)

# Set up such that XCode organizes the files correctly
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES})

# Add Library
find_library(CL CL HINTS /usr/include)
#find_library(IBSLIB_LIB IBSlib  HINTS /usr/local/include)
find_library(IBSLIB_LIB IBSLib HINTS /usr/local/lib/)

message(STATUS "${IBSLIB_LIB}")

add_library(${PROJECT_NAME}  SHARED ${SOURCE_FILES})
#add_library(IBSLib  SHARED ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC ${IBSLIB_LIB} )

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE include/)

# Install
install(TARGETS ${PROJECT_NAME} DESTINATION lib)

# Install the headers
install(FILES include/cte DESTINATION include)

# create base directory
install(DIRECTORY include/ DESTINATION include)