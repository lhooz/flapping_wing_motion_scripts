"""script for tabulated 6DOF flapping wing motion"""

import numpy as np
from kinematic_functions import smooth_kinematic_function as sm_f
from kinematic_functions import sinusoidal_kinematic_function as si_f
from kinematics_write import kf_plotter, write_2d, write_3d

# sinumation time definition and choose functions to use
time_series_length = 10000
start_time = 0
end_time = 18
use_function = 'smooth'
# use_function = 'sinusoidal'
# write_mode = 'original'
write_mode = 'zero_start'

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

flapping_acceleration_time_coefficient = 0  # between 0 and 1
pitching_time_coefficient = 3.3  # between 0 and inf

# additional kinematic control parameters for sinusoidal functions
flapping_angular_velocity_amplitude = 0.5 * 180 / np.pi

flapping_acceleration_time_fraction = 0.24
pitching_time_fraction = 0.24
# ----------------------------------------------
t = np.linspace(start_time, end_time, time_series_length)

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
t_1st_cycle = [t1 for t1 in t if t1 <= 1 / flapping_wing_frequency]
# angles_to_plot = ['phi', 'dphi', 'alf', 'dalf']
angles_to_plot = ['phi', 'alf']

kf_plotter(t_1st_cycle, kinematic_angles, angles_to_plot)
# ----------------------------------------------
write_2d(t, section_location, flapping_wing_frequency, kinematic_angles,
         write_mode)
write_3d(t, flapping_wing_frequency, kinematic_angles, write_mode)
