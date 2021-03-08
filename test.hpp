#include "fftw.h"
#include "myfft.h"
#include "signalgen.h"
#include <complex>
#include <iostream>
#include <vector>

#define VERBOSE 0

using fftout_t = std::vector<std::complex<double>>;

// Tolerance is introduced to compensate for noise due to floating-point math.
static constexpr double g_tolerance = 0.005;

void dump_fftout (const fftout_t& out) {
    for (size_t i = 0; i < out.size (); i++)
        std::cout << i + 1 << ": " << out[i] << "\n";
}

bool is_close (const fftout_t& out1, const fftout_t& out2) {
    if (out1.size () != out2.size ()) {
        std::cout << "Varying sizes: "
                  << "out1: " << out1.size () << ", out2: " << out2.size () << "\n";
        return false;
    }

    size_t errors = 0;
    for (size_t i = 0; i < out1.size (); i++) {
        auto out1_real = fabs (out1[i].real ());
        auto out1_imag = fabs (out1[i].imag ());

        auto out2_real = fabs (out2[i].real ());
        auto out2_imag = fabs (out2[i].imag ());

        auto real_diff = fabs (out1_real - out2_real);
        auto imag_diff = fabs (out1_imag - out2_imag);

        if (real_diff > g_tolerance) {
#if VERBOSE
            std::cout << "Index: " << i << ", Tolerance: " << g_tolerance << ", real_diff: " << real_diff << "\n";
#endif
            errors++;
        }

        if (imag_diff > g_tolerance) {
#if VERBOSE
            std::cout << "Index: " << i << ", Tolerance: " << g_tolerance << ", imag_diff: " << imag_diff << "\n";
#endif
            errors++;
        }
    }

    if (errors > 0)
        std::cout << "Errors: " << errors << "\n";
    return errors == 0;
}

// Generate random floating point sequences of size in the range [2^2, 2^16].
// Compute FFT using my implementation and FFTW and verify each output sample is close enough.
// Repeat the test `tests` times.
void test0 () {
    int tests = 5;
    bool success = true;
    for (int i = 0; i < tests; i++) {
        size_t start = 1;
        size_t steps = 16;
        size_t stop = start + steps;

        for (size_t j = start; j < stop; j++) {
            size_t size = 2u << j;
            const std::vector<double> samples = generate (size);

            MyFFT myfft;
            FFTW fftw;

            auto out1 = myfft.forward_transform (samples);
            auto out2 = fftw.forward_transform (samples);

            bool close = is_close (out1, out2);
            if (close) {
#if VERBOSE
                std::cout << "Size: " << size << " -- FFTs Match :)\n\n";
#endif
            } else {
                success = false;
#if VERBOSE
                std::cout << "Size: " << size << " -- FFTs Don't Match :(\n\n";
#endif
            }
        }
    }

    if (success)
        std::cout << "Test0 succeeded :)\n";
    else
        std::cout << "Test0 failed :(\n";
}


// Input contains 128 samples of a sine wave sampled at 512Hz, phase diff 0 and frequency 16Hz.
// Frequency spacing is fs/N = 512 / 128 = 4Hz. This setting won't bleed energy at neighbouring bins.
// Frequency components must be (very close to) zero everywhere except at bin 4 and (N - 4) = 124.
void test1 () {

    auto samples = fetch_from_file ("../sample/testcase1.sam");
    auto size = samples.size ();
    MyFFT myfft;
    auto fdomain = myfft.forward_transform (samples);

    // Peak amplitude at m = 4 and m = size-4.
    double abs_expect_real = 0;
    double abs_expect_imag = 64;
    std::complex<double> expect = {-0, -64};

    bool success = true;
    for (size_t i = 0; i < size; i++) {
        auto fsample = fdomain[i];
        double abs_real = fabs (fsample.real ());
        double abs_imag = fabs (fsample.imag ());

        if (i == 4 || i == size - 4) {
            success = (abs_real - abs_expect_real) < g_tolerance && (abs_imag - abs_expect_imag) <= g_tolerance;
#if VERBOSE
            if (!success) {
                std::cout << "Peak energy expected at m = " << i << " (Tolerance: " << g_tolerance << ")\n";
                std::cout << "Index: " << i << ", Expected: " << expect << ", Found: " << fsample << "\n";
            }
#endif
        } else {
            success = abs_real < g_tolerance && abs_imag < g_tolerance;
#if VERBOSE
            if (!success) {
                std::cout << "Energy not expected at bin " << i << " (Tolerance: " << g_tolerance << ")\n";
                std::cout << "Index: " << i << ", Expected: 0, Found: " << fsample << "\n";
            }
#endif
        }
    }

    if (success)
        std::cout << "Test1 succeeded :)\n";
    else
        std::cout << "Tes1 failed :(\n";
}

// Input contains 128 samples of a signal composed of a sine and a cosine wave (16 and 32Hz) sampled at 512Hz.
// Signals have amplitudes 1 and 2 respectively and 0 phase.
// Frequency spacing is fs/N = 512 / 128 = 4Hz. This setting won't bleed energy at neighbouring bins.
// Frequency components must be (very close to) zero everywhere except at bin 4, 8, (N - 4) = 124, and (N-8) = 122.
void test2 () {
    auto samples = fetch_from_file ("../sample/testcase2.sam");
    auto size = samples.size ();
    MyFFT myfft;
    auto fdomain = myfft.forward_transform (samples);

    // Amplitude spike at m = [4, 8, 124, 122]
    double abs_expect_real_at_4 = 0;
    double abs_expect_imag_at_4 = 64;
    double abs_expect_real_at_8 = 0;
    double abs_expect_imag_at_8 = 128;
    std::complex<double> expect_at_4 = {-0, -64};
    std::complex<double> expect_at_8 = {-0, -128};
    bool success = true;
    for (size_t i = 0; i < size; i++) {
        auto fsample = fdomain[i];
        double abs_real = fabs (fsample.real ());
        double abs_imag = fabs (fsample.imag ());

        switch (i) {
            case 4:
            case 124:
                success = (abs_real - abs_expect_real_at_4) < g_tolerance && (abs_imag - abs_expect_imag_at_4) < g_tolerance;
#if VERBOSE
                if (!success) {
                    std::cout << "Energy expected at m = " << i << " (Tolerance: " << g_tolerance << ")\n";
                    std::cout << "Index: " << i << ", Expected: " << expect_at_4 << ", Found: " << fsample << "\n";
                }
#endif
                break;
            case 8:
            case 120:
                success = (abs_real - abs_expect_real_at_8) < g_tolerance && (abs_imag - abs_expect_imag_at_8) < g_tolerance;
#if VERBOSE
                if (!success) {
                    std::cout << "Energy expected at m = " << i << " (Tolerance: " << g_tolerance << ")\n";
                    std::cout << "Index: " << i << ", Expected: " << expect_at_8 << ", Found: " << fsample << "\n";
                }
#endif
                break;
            default:
                success = abs_real < g_tolerance && abs_imag < g_tolerance;
#if VERBOSE
                if (!success) {
                    std::cout << "Energy not expected at m = " << i << "(Tolerance: " << g_tolerance << ")\n";
                    std::cout << "Index: " << i << ", Expected: 0, Found: " << fsample << "\n";
                }
#endif
        }
    }

    if (success)
        std::cout << "Test2 succeeded :)\n";
    else
        std::cout << "Test2 failed :(\n";
}

// Input contains 128 samples of a sine wave sampled at 512Hz, phase diff 0 and frequency 15Hz.
// Frequency spacing is fs/N = 512 / 128 = 4Hz. This setting WILL bleed energy at neighbouring bins.
// All bins must have significant energy, i.e. magnitude greater than then acceptable tolerance.
// Frequency components must be maximum at m = 4 and m = N-4 = 124.
// To perform this test, we just make sure that our output is close to that of FFTW.
void test3 () {

    auto samples = fetch_from_file ("../sample/testcase1.sam");
    MyFFT myfft;
    FFTW fftw;
    auto out1 = myfft.forward_transform (samples);
    auto out2 = fftw.forward_transform (samples);
    bool close = is_close (out1, out2);
    if (close)
        std::cout << "Test3 succeeded :)\n";
    else
        std::cout << "Test3 failed :(\n";
}

void test_ifft ()
{
    unsigned degree = 15;
    size_t size = 2u << degree;
    MyFFT myfft;
    auto samples = generate(size);
    auto fdomain = myfft.forward_transform (samples);
    auto tdomain = myfft.inverse_transform(fdomain);
    double maxnoise = 0;
    for (size_t i = 0; i < size; i++) {
        double noise = fabs (fabs (samples[i]) - fabs (tdomain[i]));
        if (noise > maxnoise)
            maxnoise = noise;
    }
    bool success = maxnoise < 2e-12;
    if (success)
        std::cout << "IFFT test succeeded :)\n";
    else std::cout << "IFFT test failed :(\n";
}

void test () {
    test0 ();
    test1 ();
    test2 ();
    test3 ();
    test_ifft();
}