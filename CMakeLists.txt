cmake_minimum_required(VERSION 3.16)
project(fft)

set(CMAKE_CXX_STANDARD 14)

add_executable(fft main.cpp myfft.cpp fftw.cpp signalgen.cpp)
target_link_libraries(fft -lfftw3 -lm)