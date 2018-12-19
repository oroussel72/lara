cmake_minimum_required(VERSION 2.6)
project(lara)

add_executable(lara main.cpp)

install(TARGETS lara RUNTIME DESTINATION bin)
