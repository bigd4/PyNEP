cmake_minimum_required(VERSION 3.15)
project(GenerateNew LANGUAGES CXX)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
add_subdirectory(../pybind11 pybind11_binary)

pybind11_add_module(nep src/pynep.cpp src/nep.cpp)
