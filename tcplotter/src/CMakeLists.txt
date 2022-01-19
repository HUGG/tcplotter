cmake_minimum_required(VERSION 3.10)

# Define project name
project(tcplotter)

# add the MathFunctions library
add_subdirectory(RDAAM_He)
add_subdirectory(ketch_aft)

# 
install(TARGETS RDAAM_He ketch_aft DESTINATION bin)

# add the executable
#add_executable(Tutorial tutorial.cxx)

#target_link_libraries(Tutorial PUBLIC MathFunctions)

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
#target_include_directories(Tutorial PUBLIC
#                          "${PROJECT_BINARY_DIR}"
#                          "${PROJECT_SOURCE_DIR}/MathFunctions"
#                          )
