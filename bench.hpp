#include <iostream>
#include <cmath>
#include <vector>
#include "signalgen.h"
#include "myfft.h"

void bench ()
{
    std::cout << "----------------- Bench -----------------\n";
    int degree = 18;
    size_t size = pow(2, degree);
    std::cout << "Size: " << size << "\n";
    const std::vector<double> samples = generate(size);
    MyFFT myfft;
    auto output = myfft.transform(samples);
    std::cout << "-----------------------------------------\n";
}