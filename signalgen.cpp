#include "signalgen.h"
#include <iostream>
#include <cmath>

void dump_samples (const std::vector<double>& samples)
{
    std::cout << "Size: " << samples.size() << "\n";
    for (double sample : samples)
        std::cout << sample << "\n";
}

std::vector<double> generate_sine_samples (size_t size)
{
    double amplitude = 1;
    double frequency = 5000;
    double sample_rate = 3 * frequency;
    double sample_interval = 1 / sample_rate;
    // frequency * (2 * pi) * sample_interval

    std::vector<double> samples;
    samples.reserve(size);
    for (size_t i = 0; i < size; i++) {
        double arg = frequency * (2 * M_PI) * (i * sample_interval);
        samples.emplace_back (amplitude * sin(arg));
    }
    return samples;
}

std::vector<double> generate (size_t size)
{
    auto tmp = generate_sine_samples(size);
    //dump_samples(tmp);
    return tmp;
}