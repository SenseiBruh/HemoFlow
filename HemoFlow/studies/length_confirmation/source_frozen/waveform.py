"""Mean-preserving inlet waveforms shared by solver and viewer.

The asymmetric shape is illustrative: it is not a measured patient waveform.
Both pulse shapes have zero cycle mean and range [-1, 1].
"""
import math
import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def pulse_value(time, amplitude, frequency, shape_code):
    if shape_code == 0:
        return 1 + amplitude * math.sin(2 * math.pi * frequency * time)
    phase = (time * frequency) % 1.0
    rise = 0.14
    if phase < rise:
        wave = -math.cos(math.pi * phase / rise)
    else:
        wave = math.cos(math.pi * (phase - rise) / (1 - rise))
    return 1 + amplitude * wave


def pulse_curve(config, times):
    times = np.asarray(times)
    phase = np.mod(times * config.heart_rate / 60, 1.0)
    if config.pulse_shape == "sine":
        wave = np.sin(2 * np.pi * phase)
    else:
        wave = np.where(phase < .14, -np.cos(np.pi * phase / .14),
                        np.cos(np.pi * (phase - .14) / .86))
    return 1 + config.pulsatility_percent / 100 * wave


def phase_label(config, time):
    if config.pulsatility_percent == 0:
        return "STEADY FLOW"
    return "SYSTOLE" if (time * config.heart_rate / 60) % 1 < .35 else "DIASTOLE"
