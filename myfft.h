#pragma once

#include <vector>
#include "mycomplex.hpp"

class MyFFT
{
private:
    std::vector<MyComplex<double>> m_twiddles;
    std::vector<MyComplex<double>> m_samples;
    std::vector<size_t> m_bit_reverse_table;

public:
    size_t get_complement (size_t number, size_t width);
    void bit_reversal (std::vector<MyComplex<double>>& samples);
    MyComplex<double> twiddle (size_t fft_size, size_t k);
    void synth (std::vector<MyComplex<double>>& samples, size_t start, size_t end);
    void fft (std::vector<MyComplex<double>>& samples, size_t start, size_t end);
    void fft_iter (std::vector<MyComplex<double>>& samples);
    void compute_twiddles (size_t fft_size);
    std::vector<MyComplex<double>> forward_transform (const std::vector<double>& samples);

    std::vector<double> inverse_transform (const std::vector<MyComplex<double>>& samples);
};