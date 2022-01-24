cmake_minimum_required(VERSION 3.10)

# set the project name
project(ketch_aft)

# add math library for build
target_link_libraries(ketch_aft PRIVATE m)

# add the executable
add_executable(ketch_aft ketch.c ketch_aft.c)
