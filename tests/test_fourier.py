import numpy as np
import pytest
from pyspecdenoise.fourier import fourier_smoothing

TEST_FILE = 'SN2004gq_2004-12-12_07-12-00_FLWO-1.5m_FAST_CfA-Stripped_clipped.flm'
CHECK_FILE = 'smoothed.txt'

@pytest.fixture
def load_input_spec():
    with open(TEST_FILE) as newfile:
        lines = newfile.readlines()
        i = 0
        for line in lines:
            if line[2] != " ":
                break
            i += 1

    wave, flux = np.loadtxt(TEST_FILE, usecols=range(2), skiprows=i, unpack=True)

    return wave, flux

@pytest.fixture
def load_comparison_spec():
    wave, flux, flux_smooth = np.loadtxt(CHECK_FILE)
    return wave, flux, flux_smooth

def test_fourier_smoothing():
    input_wavelength, input_flux = load_input_spec

    wlen, new_flux = fourier_smoothing(
        input_wavelength / (1 + 0.0065), input_flux / np.mean(input_flux)
    ) # Divide by the mean so that all spectra are "normalised".

    x = load_comparison_spec

    assert x == np.stack((wlen, input_flux, new_flux))