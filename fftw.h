#pragma once
#include <vector>
#include <complex>

class FFTW
{
public:
    std::vector<std::complex<double>> forward_transform (const std::vector<double>& samples);
};
