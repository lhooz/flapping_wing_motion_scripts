"""script for tabulated 6DOF flapping wing motion"""

import numpy as np
from kinematic_functions import smooth_kinematic_function as sm_f
from kinematic_functions import sinusoidal_kinematic_function as si_f
from kinematics_write import kf_plotter, write_2d, write_3d

# sinumation time definition and choose functions to use
time_series_length_per_cycle = 1000
start_time = 0
number_of_cycles = 3
use_function = 'smooth'
# use_function = 'sinusoidal'

# ----------------------------------------------
# common kinematic parameters
flapping_wing_frequency = 1

pitching_amplitude = 45

flapping_delay_time_fraction = 0
pitching_delay_time_fraction = 0

section_location = 1
# ----------------------------------------------
# additional kinematic parameters for smooth functions
flapping_amplitude = 80

flapping_acceleration_time_coefficient = 0.97  # between 0 and 1
pitching_time_coefficient = 3  # between 0 and inf

# additional kinematic control parameters for sinusoidal functions
flapping_angular_velocity_amplitude = 0.5 * 180 / np.pi

flapping_acceleration_time_fraction = 0.24
pitching_time_fraction = 0.24
# ----------------------------------------------
t1 = np.linspace(start_time, 1/flapping_wing_frequency, time_series_length_per_cycle)
t = t1
for i in range(1, number_of_cycles):
    ti = np.delete(t1, 0) + i
    t = np.append(t, ti)

kinematic_parameters_smooth = [
    flapping_wing_frequency, flapping_amplitude, pitching_amplitude,
    flapping_acceleration_time_coefficient, pitching_time_coefficient,
    flapping_delay_time_fraction, pitching_delay_time_fraction
]

kinematic_parameters_sinusoidal = [
    flapping_wing_frequency, flapping_angular_velocity_amplitude,
    pitching_amplitude, flapping_acceleration_time_fraction,
    pitching_time_fraction, flapping_delay_time_fraction,
    pitching_delay_time_fraction
]

if use_function == 'smooth':
    kinematic_angles = sm_f(t, kinematic_parameters_smooth)
elif use_function == 'sinusoidal':
    kinematic_angles = si_f(t, kinematic_parameters_sinusoidal)

# plotting kinematic angles
# angles_to_plot = ['phi', 'dphi', 'alf', 'dalf']
angles_to_plot = ['dphi', 'dalf']

kf_plotter(t, kinematic_angles, time_series_length_per_cycle, angles_to_plot)
# ----------------------------------------------
write_2d(t, section_location, time_series_length_per_cycle, kinematic_angles)
write_3d(t, time_series_length_per_cycle, kinematic_angles)
