cmake_minimum_required(VERSION 3.21)
project(voronoi)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_GENERATOR "Unix Makefiles" CACHE INTERNAL "" FORCE)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

add_executable(voronoi voronoi.cpp)

find_package(CGAL)
target_link_libraries(voronoi CGAL::CGAL)
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(voronoi ${MPI_C_LIBRARIES})
