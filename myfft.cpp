#include "myfft.h"
#include <bitset>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

using u8 = unsigned char;
using u32 = unsigned int;

#define PI 3.14159265

size_t MyFFT::get_complement (size_t number, size_t width) {
    // Unpack bits
    std::vector<u8> bits;
    bits.reserve (width);
    for (size_t i = 0; i < width; i++)
        bits.push_back ((number & (1u << i)) > 0);

    // Reverse bits
    for (size_t i = 0, j = width - 1; i < j; i++, j--) {
        u8 tmp = bits.at (i);
        bits[i] = bits[j];
        bits[j] = tmp;
    }

    // Compute "complement"
    size_t accum = 0;
    for (int i = 0; i < width; i++)
        accum += bits[i] * (pow (2, i));
    return accum;
}

void MyFFT::bit_reversal (std::vector<double>& samples) {
    // TODO: check if sample vector size is a power of 2.
    int size = samples.size ();
    std::cout << "Bit reversal---------------------------------\n";
    //m_bit_reverse_index.reserve(size);
    int bit_space = (int) log2 (size);
    for (int i = 0; i < size / 2; i++) {
        int complement_i = get_complement (i, bit_space);

        std::bitset<4> ibits (i);
        std::bitset<4> cbits (complement_i);

        std::cout << ibits << "(" << i << ")" <<  ", " << cbits << "(" << complement_i << ")" << "\n";
        
        auto tmp = samples[i];
        samples[i] = samples[complement_i];
        samples[complement_i] = tmp;
        //m_bit_reverse_index.emplace_back(complement_i);
    }

    std::cout << "---------------------------------\n";

    /*for (size_t i = 0; i < size; i++)
    {
        double tmp = samples[i];
        samples[i] = samples[m_bit_reverse_index[i]];
        samples[m_bit_reverse_index[i]] = tmp;
    }*/
}

std::complex<double> MyFFT::twiddle (size_t fft_size, size_t k) {
    double numerator = 1 - pow (2, log2 (fft_size / 2));
    size_t segment = -(2 * numerator);
    return m_twiddles[segment + k];
}

void MyFFT::synth (std::vector<std::complex<double>>& samples, size_t start, size_t end) {
    size_t fft_size = (end - start) + 1;
    for (size_t k = 0, m = start; k < fft_size / 2; k++, m++) {
        auto out1 = samples[m] + (twiddle (fft_size, k) * samples[m + (fft_size / 2)]);
        auto out2 = samples[m] + (twiddle (fft_size, k + (fft_size / 2)) * samples[m + (fft_size / 2)]);
        samples[m] = out1;
        samples[m + (fft_size / 2)] = out2;
    }
}

void MyFFT::fft (std::vector<std::complex<double>>& samples, size_t start, size_t end) {
    size_t mid = (start + end) / 2;
    size_t fft_size = (end - start) + 1;
    if (fft_size == 2) {
        synth (samples, start, end);
        return;
    }
    fft (samples, start, mid);
    fft (samples, mid + 1, end);
    synth (samples, start, end);
}

void MyFFT::compute_twiddles (size_t fft_size) {
    size_t size = -(2 * ((1 - pow (2, log2 (fft_size)))));
    m_twiddles.reserve (size);
    m_twiddles.emplace_back (1, 0);
    m_twiddles.emplace_back (-1, 0);
    for (size_t m = 4; m <= fft_size; m *= 2) {
        for (size_t i = 0; i < m; i++) {
            std::complex<double> exponent = {0, (-2 * PI * i) / m};
            auto tmp = exp (exponent);
            /*if (fabs (tmp.real ()) < 1e-6)
                tmp = {0, tmp.imag ()};
            if (fabs (tmp.imag ()) < 1e-6)
                tmp = {tmp.real (), 0};*/
            m_twiddles.emplace_back (tmp);
        }
    }

    std::cout << "\n----------------------------------------\n";
        for (size_t i = 0; i < m_twiddles.size(); i++)
            std::cout << "twiddle[" << i << "] = " << m_twiddles.at(i) << "\n";
        std::cout << "\n----------------------------------------\n";
}

void MyFFT::reduce_noise () {
    for (auto& sample : m_samples) {
        if (fabs (sample.real ()) < 1e-4)
            sample = {0, sample.imag ()};
        if (fabs (sample.imag ()) < 1e-4)
            sample = {sample.real (), 0};
    }
}

std::vector<std::complex<double>> MyFFT::transform (const std::vector<double>& samples) {
    compute_twiddles (samples.size ());
    // I hope this does an element-by-element copy.
    std::vector<double> samples_copy {samples};
    bit_reversal (samples_copy);
    m_samples.reserve (samples.size ());
    for (double sample : samples_copy)
        m_samples.emplace_back (sample, 0);
    fft (m_samples, 0, samples.size () - 1);
    //reduce_noise ();
    return m_samples;
}
