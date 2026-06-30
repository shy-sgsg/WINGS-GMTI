# range resize validation

4096 -> 11820 uses FFT, fftshift, centered zero padding, ifftshift, and FFTW backward IFFT scaled by 1/src_len.

## constant

- input rms: 1
- output rms: 1
- output max_abs: 1
- nan/inf: false/false

## single_tone

- input rms: 1
- output rms: 1
- output max_abs: 1
- nan/inf: false/false

## random

- input rms: 0.818364
- output rms: 0.818364
- output max_abs: 1.97621
- nan/inf: false/false

This is an engineering interpolation for stage 1. It does not convert old 38 us LFM physics into strict new 130 us LFM echoes.
