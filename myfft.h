#pragma once

#include <vector>
#include <complex>

class MyFFT
{
private:
    std::vector<std::complex<double>> m_twiddles;
    std::vector<std::complex<double>> m_samples;
    std::vector<size_t> m_bit_reverse_table;

public:
    size_t get_complement (size_t number, size_t width);
    void bit_reversal (std::vector<double>& samples);
    std::complex<double> twiddle (size_t fft_size, size_t k);
    void synth (std::vector<std::complex<double>>& samples, size_t start, size_t end);
    void fft (std::vector<std::complex<double>>& samples, size_t start, size_t end);
    void fft_iter (std::vector<std::complex<double>>& samples);
    void compute_twiddles (size_t fft_size);
    std::vector<std::complex<double>> transform (const std::vector<double>& samples);
};