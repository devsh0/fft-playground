cmake_minimum_required(VERSION 3.16)
project(fft)

#[[set(CMAKE_CXX_FLAGS -O2)]]
set(CMAKE_CXX_STANDARD 17)

add_executable(fft main.cpp myfft.cpp fftw.cpp signalgen.cpp fftest.hpp bench.hpp mycomplex.hpp cplxtest.hpp)
target_link_libraries(fft -lfftw3 -lm)