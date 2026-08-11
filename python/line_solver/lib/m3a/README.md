# M3A (3rd Moment Approximation) Library

This is a Python port of the M3A library originally implemented in MATLAB
and located at `matlab/lib/m3a/` in the LINE solver repository.

M3A provides methods for fitting and compressing Markov-modulated processes
using moment matching up to the third order. It includes:

- **Fitting**: Construction of MAP/PH representations from moments
- **Compression**: Reduction of MAP/PH representations to smaller state spaces
- **M3PP utilities**: Analysis and construction of Markov-Modulated Poisson Processes
- **Moment utilities**: Computation and manipulation of joint moments

## Origin

This code is a derivative work ported from MATLAB to Python for use in the
native Python implementation of the LINE solver. The original MATLAB source
is located at `matlab/lib/m3a/` in the LINE solver repository.

## License

See `LICENSE` in this directory. BSD 3-Clause License.
