"""Physical and numerical constants shared by the whole study."""

#: Speed of sound (m/s), about 20 degC. The simulation and the correction use the same
#: value so that the study isolates the geometry, not a temperature mismatch.
C = 343.0

#: Density of air (kg/m^3) at 20 degC.
RHO = 1.204

#: Reference pressure for dB SPL (Pa).
P_REF = 20e-6

#: Sample rate of every simulated impulse response (Hz).
FS = 48_000

#: Length of every simulated impulse response (samples): 341 ms at 48 kHz, far longer
#: than the longest propagation delay plus the decay of the drivers.
N_FFT = 16_384

#: Distance at which the driver responses are specified (m).
R_REF = 1.0
