#include "fftw.h"
#include "myfft.h"
#include <iostream>

using fft_type = std::vector<std::complex<double>>;
bool is_close (const fft_type& out1, const fft_type &out2, double tolerance) {
    if (out1.size () != out2.size ()) {
        std::cout << "Varying sizes: "
                  << "out1: " << out1.size () << ", out2: " << out2.size () << "\n";
        return false;
    }

    bool okay = true;
    for (size_t i = 0; i < out1.size (); i++) {
        auto out1_real = out1[i].real ();
        auto out1_imag = out1[i].imag ();

        auto out2_real = out2[i].real ();
        auto out2_imag = out2[i].imag ();

        auto real_diff = fabs (out1_real - out2_real);
        auto imag_diff = fabs (out1_imag - out2_imag);

        if (real_diff > tolerance) {
            std::cout << "Tolerance: " << tolerance << ", real_diff: " << real_diff;
            okay = false;
        }
        if (imag_diff > tolerance) {
            std::cout << "Tolerance: " << tolerance << ", imag_diff: " << imag_diff;
            okay = false;
        }
    }

    return okay;
}

int main () {
    std::vector<double> samples = {0.3535, 0.3535, 0.6464, 1.0607, 0.3535, -1.0607, -1.3535, -0.3535};
    MyFFT myfft;
    FFTW fftw;

    auto out1 = myfft.transform (samples);
    auto out2 = fftw.transform(samples);
    double tolerance = 0.01;
    std::cout << "Is Close: " << is_close(out1, out2, tolerance);


    return 0;
}
