cmake_minimum_required(VERSION 3.10)

# Define project name
project(tcplotter)

# Define project subdirectories to compile
add_subdirectory(RDAAM_He)
add_subdirectory(ketch_aft)

# Define programs to install
install(TARGETS RDAAM_He ketch_aft DESTINATION bin)
