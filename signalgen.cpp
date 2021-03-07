#include <ctime>
#include "signalgen.h"
#include <iostream>
#include <cmath>

static bool seeded = false;

void seed ()
{
    srand(std::time(nullptr));
}

void dump_samples (const std::vector<double>& samples)
{
    std::cout << "Size: " << samples.size() << "\n";
    for (double sample : samples)
        std::cout << sample << "\n";
}

int range_rand (int min, int max)
{
    return min + (rand() % static_cast<int>(max - min + 1));
}

std::vector<double> generate_samples (size_t size)
{
    double amplitude = range_rand(1, 10);
    double frequency = range_rand(1, 25000);
    double sample_rate = range_rand(2, 4) * frequency;
    double sample_interval = 1 / sample_rate;

    std::vector<double> samples;
    samples.reserve(size);
    for (size_t i = 0; i < size; i++) {
        double arg = frequency * (2 * M_PI) * (i * sample_interval);
        samples.emplace_back (amplitude * sin(arg));
    }

    amplitude = range_rand(1, 10);
    frequency = range_rand(1, 25000);
    for (size_t i = 0; i < size; i++) {
        double arg = frequency * (2 * M_PI) * (i * sample_interval);
        samples[i] += (amplitude * cos(arg));
    }

    return samples;
}

void test_samples_not_zero(const std::vector<double>& samples)
{
    double tol = 0.01;
    double accum = 0;
    for (double sample : samples)
        accum += sample;
    if (accum < tol)
        std::cout << "It's all zero!\n";
    else std::cout << "Input okay!\n";
}

std::vector<double> generate (size_t size)
{
    if (!seeded)
    {
        seed();
        seeded = true;
    }

    std::vector<double> super(size, 0);

    for (size_t i = 0; i < 10; i++) {
        auto tmp = generate_samples (size);
        for (size_t m = 0; m < size; m++)
            super[m] += tmp[m];
    }
    test_samples_not_zero(super);
    return super;
}