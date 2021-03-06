#include "fftw.h"
#include <fftw3.h>

std::vector<std::complex<double> > FFTW::transform(const std::vector<double>& samples)
{
    size_t size = samples.size();
    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * size);
    out = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * size);
    auto plan = fftw_plan_dft_1d (size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (size_t i = 0; i < size; i++)
    {
        in[i][0] = samples[i];
        in[i][1] = 0; // Imaginary part always zero
    }

    fftw_execute(plan);

    std::vector<std::complex<double>> output;
    output.reserve(size);
    for (size_t i = 0; i < size; i++) {
        std::complex<double> tmp = {out[i][0], out[i][1]};
        output.emplace_back (tmp);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    return output;
}