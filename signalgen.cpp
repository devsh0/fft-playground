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
    std::vector<double> samples;
    samples.reserve(size);
    double scale = range_rand(1, 100);
    for (size_t i = 0; i < size; i++)
    {
        double numer = range_rand(1, 1000);
        double denom = range_rand(100, 2000);
        int sign = (range_rand(1, 4) > 2) ? 1 : -1;
        double sample = scale * sign * (numer / denom);
        samples.emplace_back(sample);
    }
    return samples;
}

std::vector<double> generate (size_t size)
{
    if (!seeded)
    {
        seed();
        seeded = true;
    }
    auto samples = generate_samples(size);
    return samples;
}