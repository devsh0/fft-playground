#include "fftw.h"
#include "fftw3.h"
#include "myfft.h"
#include "signalgen.h"
#include <cmath>
#include <iostream>
#include <vector>

#define FFTOOL MyFFT

void bench ()
{
    std::cout << "----------------- Bench -----------------\n";
    int degree = 20;
    size_t size = pow(2, degree);
    std::cout << "Size: " << size << "\n";
    const std::vector<double> samples = generate(size);
    FFTOOL fft;
    auto output = fft.forward_transform (samples);
    std::cout << "-----------------------------------------\n";
}