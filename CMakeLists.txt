cmake_minimum_required(VERSION 3.16)
project(MonteCarloPricer LANGUAGES CXX)

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenMP (optional but used if available)
find_package(OpenMP)

# Define the executable
add_executable(pricer
    src/main.cpp
)

# Include our headers
target_include_directories(pricer PRIVATE include)

# Link OpenMP if found
if(OpenMP_CXX_FOUND)
    target_link_libraries(pricer PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found — parallel build enabled")
else()
    message(WARNING "OpenMP not found — running single-threaded")
endif()