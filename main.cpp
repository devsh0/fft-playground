#include "fftw.h"
#include "myfft.h"
#include "signalgen.h"
#include <iostream>

using fftout_t = std::vector<std::complex<double>>;
bool is_close (const fftout_t& out1, const fftout_t& out2, double tolerance) {
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
            std::cout << "Index: " << i << ", Tolerance: " << tolerance << ", real_diff: " << real_diff << "\n";
            okay = false;
        }
        if (imag_diff > tolerance) {
            std::cout << "Index: " << i << ", Tolerance: " << tolerance << ", imag_diff: " << imag_diff << "\n";
            okay = false;
        }
    }

    return okay;
}

void dump_fftout_side_by_side (const fftout_t& out1, const fftout_t& out2) {
    if (out1.size () != out2.size ())
        std::cerr << "Whoops! Sizes don't seem to match up!\n";
    for (size_t i = 0; i < out1.size (); i++) {
        std::cout << "[i: " << i << ", out1: {" << out1[i].real () << ", " << out1[i].imag () << "}, "
                  << "out2: {" << out2[i].real () << ", " << out2[i].imag () << "}]\n";
    }
}

int main () {
    for (size_t i = 1; i < 4; i++) {
        size_t size = 2u << i;
        //std::cout << "Size: " << size << "\n";
        const std::vector<double> samples = generate (size);
        MyFFT myfft;
        FFTW fftw;

        auto out1 = myfft.transform (samples);
        auto out2 = fftw.transform (samples);
        double tolerance = 0.0001;

        bool close = is_close (out1, out2, tolerance);
        if (close)
            std::cout << "Size: " << size << " -- FFTs Match :)\n";
        else
            std::cout << "Size: " << size << " -- FFTs Don't Match :(\n";
        //dump_fftout_side_by_side (out1, out2);
    }

    return 0;
}