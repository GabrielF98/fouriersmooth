"""
A python module to perform Fourier denoising of spectra.
This module is based on the algorithm presented in Liu et al. (2016).
This paper is available at: https://ui.adsabs.harvard.edu/abs/2016ApJ...827...90L/abstract).
An IDL version of this code was made available by the original developers at: https://github.com/metal-sn/SESNspectraLib/blob/master/SNspecFFTsmooth.pro
This is the first publicly available implementation of this algorithm written in Python.
"""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def power_law(x, a, b):
    """Power law fit."""
    return a * (x**b)

def create_stddev_arr(flux, smooth_flux, wave, width):
    """
    Produces a "residual" array containing the residuals between the input spectrum and a smoothed spectrum.
    Computes the standard deviation of the residual array within a symmetric window at each pixel.
    Returns this standard deviation array.

    :param flux: Spectrum flux array.
    :type flux: numpy.ndarray
    :param smooth_flux: Smoothed spectrum flux array.
    :type smooth_flux: numpy.ndarray
    :param wave: Spectrum wavelength array.
    :type wave: numpy.ndarray
    :param width: Width of the moving window. Measured in Angstroms.
    :type width: float
    :return: An array of standard deviations of the bins centered on each wavelength in the input spectrum.
    :rtype: numpy.ndarray
    """
    f_resid = flux - smooth_flux # compute residual array.
    f = pd.Series(f_resid)
    bin_size = int(width / (wave[1] - wave[0]) + 1) # compute (integer) number of samples per width.
    g = f.rolling(bin_size).std().to_numpy() # Take rolling average standard deviation.
    return g

def fourier_smoothing(wlen, flux, k_high:int=300, k_low:int=3):
    """
    Implements Liu et al 2016 smoothing in python. See their appendix B.
    Available at: https://ui.adsabs.harvard.edu/abs/2016ApJ...827...90L/abstract
    The goal is to remove the Fourier components in which the noise dominates the spectrum.

    Procedure:
    - Rebin the spectrum on a log-wavelength axis.
    - Resample the spectrum so the bins are the same width; use the smallest dispersion as the bin width.
    - Take the FFT of the flux.
    - Define the range of wavenumbers/velocities for spectral features (see note); the FT indices are determined using k_low and k_high.
    - Fit the magnitude (M) spectrum with a power law between k_low and k_high.
    - Compute MEAN(M).
    - k_noise is the point of intersection between the power law fit an MEAN(M).
    - Set M = 0 for k>k_noise.
    - Invert FFT.
    - Resample spectrum to the original linear grid.


    Note:
    k is related to the velocity of spectral features in the SN spectrum by k = c/v.
    Thus k can be chosen to exclude high and low velocity features that are likely not due to the SN.
    The default values of k are k = 300 (3000 km/s) and k = 3 (100000 km/s) as suggested in Liu et al. 2016.

    :param wlen: Wavelength array of the spectrum.
    :type wlen: numpy.ndarray
    :param flux: Array of spectral flux.
    :type flux: numpy.ndarray
    :param k_high: Upper k value for smoothing.
    :type k_high: int
    :param k_low: Lower k value for smoothing.
    :type k_low: int
    """

    # Convert to log wavelength space
    # 1 - Convert the wavelength array to an evenly spaced array in log wavelength space.
    log_wlen = np.log(wlen)
    bin_width = min(np.diff(log_wlen))  # Smallest 'dispersion' in the spectrum
    rs_log_wlen = np.arange(log_wlen[0], log_wlen[-1], bin_width)

    # 2 - Resample the flux on the new wavelength scale
    rs_flux = interp1d(log_wlen, flux)(rs_log_wlen)

    # Compute the Fourier transform
    # 1 - Perform an FFT on the resampled flux
    flux_fft = fft(rs_flux)
    flux_fft_abs = np.abs(flux_fft)
    N = len(flux_fft_abs)
    f_samp = 1 / bin_width
    k = np.fft.fftfreq(len(rs_flux), 1 / f_samp)

    max_index = 0
    if N % 2 == 0:
        print("N is even")
        max_index = N // 2
    else:
        max_index = N // 2 + 1

    fig, ax = plt.subplots(1, 2, figsize=[15, 5])
    ax[0].set_xlabel("Wave number k [1/ln($\AA$)]")
    ax[0].set_ylabel("Magnitude")
    ax[0].plot(
        k[:max_index], flux_fft_abs[:max_index], marker="", linewidth=0.3, color="k"
    )
    ax[0].set_yscale("log")
    fig.tight_layout()

    # Find out the indices where k_low<k<k_high
    k_high_index = max(np.where(k[:max_index] < k_high)[0])
    k_low_index = min(np.where(k[:max_index] > k_low)[0])

    # Compute the mean magnitude for k_low<k<k_high
    M = np.mean(flux_fft_abs[k_low_index : k_high_index + 1])
    ax[0].hlines(M, k_low, k_high, linestyles="--", colors="g", label="Avg. mag.")
    ax[0].set_xlim([0, np.max(k)])

    # Power law fit to the array where k>3
    popt, _ = curve_fit(
        power_law, k[k_low_index:max_index], flux_fft_abs[k_low_index:max_index]
    )
    ax[0].plot(
        k[k_low_index:max_index],
        power_law(k[k_low_index:max_index], *popt),
        color="r",
        label="Power law fit",
    )

    # Determine k_noise
    k_noise_idx = np.where(power_law(k[1:max_index], *popt) < M)[0][0]
    print(M, k[k_noise_idx], power_law(k[k_noise_idx], *popt))
    ax[0].axvline(
        k[k_noise_idx],
        color="orange",
        label="k$_n$ {:.0f}".format(k[k_noise_idx]),
    )
    ax[0].legend()
    print(k)

    # Set the noise FFT coefficients to 0
    flux_fft[k_noise_idx:-k_noise_idx] = 0

    # # Compute the iFFT
    new_signal = np.real(ifft(flux_fft))

    # Rebin on the original wavelength scale
    new_flux = interp1d(rs_log_wlen, new_signal, fill_value="extrapolate")(log_wlen)
    ax[1].plot(
        np.exp(log_wlen), new_flux, alpha=0.5, color="C0", label="Original spectrum"
    )
    ax[1].step(wlen, flux, alpha=0.5, color="C1", label="Smoothed spectrum")
    ax[1].set_ylim([max(-4, np.min(flux) - 0.1), min(np.max(flux), 4)])
    ax[1].set_xlabel("Rest wavelength [$\AA$]")
    ax[1].set_ylabel("Relative Flux")
    plt.legend()

    plt.tick_params(axis="both", direction="in", which="both")

    # plt.show()
    # plt.close()
    plt.subplots_adjust(wspace=0.1)
    plt.savefig("Smoothingexample.pdf", bbox_inches="tight")

    return wlen, new_flux