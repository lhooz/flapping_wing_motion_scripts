"""script for tabulated 6DOF flapping wing motion"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from scipy.spatial.transform import Rotation as R

from mostafa_function import kinematic_function as kf

time_series_length = 10000
start_time = 0
end_time = 80

t = np.linspace(start_time, end_time, time_series_length)

# kinematic parameters for default kinematic functions
# ----------------------------------------------
flapping_frequency = 1

flapping_amplitude = 80
pitching_amplitude = 45

flapping_acceleration_time_coefficient = 0.5  # between 0 and 1
pitching_time_coefficient = 3.3  # between 0 and inf

flapping_delay_time_fraction = 0
pitching_delay_time_fraction = -0.1

# kinematic control parameters for mostafa_function
# ----------------------------------------------
# flapping_wing_frequency = 0.1

# flapping_angular_velocity_amplitude = 0.5
# pitching_amplitude = 45 * np.pi / 180

# flapping_acceleration_time_fraction = 0.24
# pitching_time_fraction = 0.24

# flapping_delay_time_fraction = 0
# pitching_delay_time_fraction = -0.1
# ----------------------------------------------

# default kinematic functions for phi and alf:
# -----------------------------------------------
def phi(x):
"""flapping motion function"""
if flapping_acceleration_time_coefficient == 0:
phit = flapping_amplitude * np.cos(
2 * np.pi * flapping_frequency *
(x - flapping_delay_time_fraction / flapping_frequency))
else:
phit = flapping_amplitude / np.arcsin(
flapping_acceleration_time_coefficient) * np.arcsinin(
flapping_acceleration_time_coefficient * np.cos(
2 * np.pi * flapping_frequency *
(x - flapping_delay_time_fraction / flapping_frequency)))
return phit

def alf(x):
"""pitching motion function"""
if pitching_time_coefficient == 0:
alft = pitching_amplitude * np.sin(2 * np.pi * flapping_frequency *
(x - pitching_delay_time_fraction / flapping_frequency))
else:
alft = pitching_amplitude / np.tanh(pitching_time_coefficient) * np.tanh(
pitching_time_coefficient * np.sin(2 * np.pi * flapping_frequency *
(x - pitching_delay_time_fraction / flapping_frequency))
return alft
# -------------------------------------------------

# mostafa kinematic functions
# -------------------------------------------------
# def dphi(x):
# """flapping motion angular velocity function"""
# return kf(flapping_wing_frequency, flapping_angular_velocity_amplitude,
# flapping_acceleration_time_fraction,
# flapping_delay_time_fraction, x)

# flapping_amplitude = integrate.quad(lambda x: np.abs(dphi(x)), 0,
# 1 / flapping_wing_frequency)[0]

# # print(flapping_amplitude * 57.3)

# dphi_int = []
# for ti in t:
# dphi_int.append(dphi(ti))

# def phi(x):
# """flapping motion function"""
# t_int = [tx for tx in t if tx <= x]
# time_array_length = len(t_int)
# # print(time_array_length, dphi_int[0:time_array_length])
# return -flapping_amplitude / 4 + integrate.simps(
# dphi_int[0:time_array_length], t_int)

# def alf(x):
# """pitching motion function"""
# return kf(flapping_wing_frequency, pitching_amplitude,
# pitching_time_fraction, pitching_delay_time_fraction, x)
# -------------------------------------------------

phiplot = []
alfplot = []
for ti in t:
    phiplot.append(phi(ti))
    alfplot.append(alf(ti))

plt.plot(t, phiplot, t, np.pi / 2 - np.abs(alfplot))
plt.show()

# r = R.from_euler('zyx', [-30, -20, 0], degrees=True)
# r = r.as_rotvec()
# test_vector = [0, 1, 0]
r_angle = []
for ti in t:
    euler_anglesi = [phi(ti), alf(ti), 0]

    roti = R.from_euler('YXZ', euler_anglesi, degrees=True)
    # result_vectori = roti.apply(test_vector)
    # print(result_vectori[2])

    r_anglei = roti.as_euler('XYZ', degrees=True)
    r_anglei = [str(r_anglei[0]), str(r_anglei[1]), str(r_anglei[2])]

    r_angle.append(r_anglei)

t = [str(ti) for ti in t]

motion = [str(time_series_length), '(']
for ti, angle_i in zip(t, r_angle):
    motioni = '(' + ti + ' ((0 0 0)' + '(' + ' '.join(angle_i) + ')))'

    motion.append(motioni)

motion.append(')')

# print(motion)
with open('6DoF.dat', 'w') as f:
    for item in motion:
        f.write("%s\n" % item)
