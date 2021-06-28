"""script for tabulated 6DOF flapping wing motion"""

import numpy as np

from kinematic_functions import sinu_continuous_kinematic_function as sic_f
from kinematic_functions import sinusoidal_kinematic_function as si_f
from kinematic_functions import smooth_kinematic_function as sm_f
from kinematics_write import (kf_plotter, write_2d, write_3d, write_exp,
                              write_iaoa, write_max_dphi)

# sinumation time definition and choose functions to use
time_series_length_per_cycle = 5001
start_time = 0
number_of_cycles = 2
# use_function = 'smooth'
use_function = 'sinu_continuous'
# use_function = 'sinusoidal'

# ----------------------------------------------
# common kinematic parameters
flapping_wing_frequency = 0.2

flapping_delay_time_fraction = 0
pitching_delay_time_fraction = 0

section_location = 3
# ----------------------------------------------
# additional kinematic parameters for smooth functions
half_flapping_amplitude = -80
half_pitching_amplitude = -45

flapping_acceleration_time_coefficient = 0  # between 0 and 1
pitching_time_coefficient = 100  # between 0 and inf or use ptf_function 'f'
ptf_coefficient = 1.6  # used when pitching_time_coefficient = 'f'
# ----------------------------------------------
# additional kinematic control parameters for sinu_continuous functions
flapping_acceleration_time_fraction = 0.125
pitching_time_fraction = 0.125

flapping_angular_velocity_amplitude = 391.06 * flapping_wing_frequency  # --degree/s--
pitching_angular_velocity_amplitude = 360 * flapping_wing_frequency / (
    2 * pitching_time_fraction)  # --degree/s--
# ---------------------------------------------
# additional kinematic control parameters for sinusoidal functions
flapping_angular_velocity_amplitude_sin = 28.27
pitching_amplitude_sin = 45

flapping_acceleration_time_fraction_sin = 0.5
pitching_time_fraction_sin = 0.5
# ----------------------------------------------
t1 = np.linspace(start_time, 1 / flapping_wing_frequency,
                 time_series_length_per_cycle)
t = t1
for i in range(1, number_of_cycles):
    ti = np.delete(t1, 0) + i / flapping_wing_frequency
    t = np.append(t, ti)

kinematic_parameters_smooth = [
    flapping_wing_frequency, half_flapping_amplitude, half_pitching_amplitude,
    flapping_acceleration_time_coefficient, pitching_time_coefficient,
    flapping_delay_time_fraction, pitching_delay_time_fraction, ptf_coefficient
]
kinematic_parameters_sinu_continuous = [
    flapping_wing_frequency, flapping_angular_velocity_amplitude,
    pitching_angular_velocity_amplitude, flapping_acceleration_time_fraction,
    pitching_time_fraction, flapping_delay_time_fraction,
    pitching_delay_time_fraction
]
kinematic_parameters_sinusoidal = [
    flapping_wing_frequency, flapping_angular_velocity_amplitude_sin,
    pitching_amplitude_sin, flapping_acceleration_time_fraction_sin,
    pitching_time_fraction_sin, flapping_delay_time_fraction,
    pitching_delay_time_fraction
]

if use_function == 'smooth':
    kinematic_angles = sm_f(t, kinematic_parameters_smooth)
elif use_function == 'sinu_continuous':
    kinematic_angles = sic_f(t, kinematic_parameters_sinu_continuous)
elif use_function == 'sinusoidal':
    kinematic_angles = si_f(t, kinematic_parameters_sinusoidal)
kinematic_angles = np.array(kinematic_angles)

#----------------save file name-----------
fAmp = np.abs(
    np.amax(kinematic_angles[:, 0]) - np.amin(kinematic_angles[:, 0]))
if use_function == 'sinu_continuous':
    save_file_name = use_function + '_phi' + '{0:.0f}'.format(
        fAmp) + '_at' + '{0:.3g}'.format(
            flapping_acceleration_time_fraction) + '_pt' + '{0:.3g}'.format(
                pitching_time_fraction)
elif use_function == 'smooth':
    save_file_name = use_function + '_phi' + '{0:.0f}'.format(
        fAmp) + '_atf' + '{0:.3g}'.format(
            flapping_acceleration_time_coefficient
        ) + '_ptc' + '{0:.3g}'.format(pitching_time_coefficient)
else:
    save_file_name == 'current'

#-----------------------------------------

# plotting kinematic angles
angles_to_plot = ['phi', 'alf']
# angles_to_plot = ['dphi', 'ddphi']
#----------------------------------------------------------------------
kf_plotter(t, kinematic_angles, angles_to_plot, time_series_length_per_cycle,
           'against_t', save_file_name + '.png')
# ----------------------------------------------
write_2d(t, section_location, kinematic_angles, time_series_length_per_cycle,
         save_file_name + '_2d.dat')
write_3d(t, kinematic_angles, time_series_length_per_cycle,
         save_file_name + '_3d.dat')
write_exp(t, kinematic_angles, time_series_length_per_cycle,
          save_file_name + '_exp.csv')
write_iaoa(kinematic_angles)
write_max_dphi(kinematic_angles)
