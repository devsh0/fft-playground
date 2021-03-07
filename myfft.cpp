#include "myfft.h"
#include <cmath>
#include <complex>
#include <vector>

using u8 = unsigned char;
using u32 = unsigned int;

#define PI 3.14159265

size_t MyFFT::get_complement (size_t number, size_t width) {
   size_t complement = 0;
   u8 rstart = width - 1;
   u8 lstart = 0;
   for (int i = 0; i < width; i++)
   {
       u8 bit = ((1u << rstart) & number) > 0u;
       complement |= bit << lstart;
       rstart--;
       lstart++;
   }
   return complement;
}

void MyFFT::bit_reversal (std::vector<double>& samples) {
    // TODO: check if sample vector size is a power of 2.
    int size = samples.size ();
    m_bit_reverse_table.reserve(size);
    int bit_space = (int) log2 (size);
    for (int i = 0; i < size; i++) {
        unsigned complement_i = get_complement (i, bit_space);
        m_bit_reverse_table.emplace_back(complement_i);
    }

    for (size_t i = 0; i < size; i++)
    {
        size_t complement = m_bit_reverse_table[i];
        if (complement > i) {
            double tmp = samples[i];
            samples[i] = samples[complement];
            samples[complement] = tmp;
        }
    }
}

std::complex<double> MyFFT::twiddle (size_t fft_size, size_t k) {
    size_t segment = -(1 - pow (2, log2 (fft_size / 2)));
    return m_twiddles[segment + k];
}

void MyFFT::synth (std::vector<std::complex<double>>& samples, size_t start, size_t end) {
    size_t fft_size = (end - start) + 1;
    for (size_t k = 0, m = start; k < fft_size / 2; k++, m++) {
        auto rhs = twiddle (fft_size, k) * samples[m + (fft_size / 2)];
        auto out1 = samples[m] + rhs;
        auto out2 = samples[m] - rhs;
        samples[m] = out1;
        samples[m + (fft_size / 2)] = out2;
    }
}

void MyFFT::fft (std::vector<std::complex<double>>& samples, size_t start, size_t end) {
    size_t mid = start + (end - start) / 2;
    size_t fft_size = (end - start) + 1;
    if (fft_size == 2) {
        synth (samples, start, end);
        return;
    }
    fft (samples, start, mid);
    fft (samples, mid + 1, end);
    synth (samples, start, end);
}

void MyFFT::fft_iter (std::vector<std::complex<double>>& samples)
{
    size_t size = samples.size();
    auto total_stages = (unsigned)log2(size);
    for (unsigned stage = 0; stage < total_stages; stage++)
    {
        size_t fft_size = 2u << stage;
        unsigned fft_count = size / fft_size;
        size_t si = 0;
        for (unsigned i = 0; i < fft_count; i++)
        {
            size_t ops = fft_size / 2;
            for (size_t m = 0; m < ops; m++)
            {
                auto twid = twiddle(fft_size, m);
                auto rhs = twid * samples[si + (fft_size / 2)];
                auto out1 = samples[si] + rhs;
                auto out2 = samples[si] - rhs;
                samples[si] = out1;
                samples[si + (fft_size / 2)] = out2;;
                si++;
            }
            si += fft_size / 2;
        }
    }
}

void MyFFT::compute_twiddles (size_t fft_size) {
    size_t size = -(1 - pow (2, log2 (fft_size)));
    m_twiddles.reserve (size);
    m_twiddles.emplace_back (1, 0);
    for (size_t m = 4; m <= fft_size; m *= 2) {
        for (size_t i = 0; i < m / 2; i++) {
            std::complex<double> exponent = {0, (-2 * PI * i) / m};
            auto tmp = exp (exponent);
            m_twiddles.emplace_back (tmp);
        }
    }
}

std::vector<std::complex<double>> MyFFT::transform (const std::vector<double>& samples) {
    compute_twiddles (samples.size ());
    std::vector<double> samples_copy {samples};
    bit_reversal (samples_copy);
    m_samples.reserve (samples.size ());
    for (double sample : samples_copy)
        m_samples.emplace_back (sample, 0);
    //fft (m_samples, 0, samples.size () - 1);
    fft_iter(m_samples);
    return m_samples;
}
