"""script for tabulated 6DOF flapping wing motion"""

import numpy as np
import matplotlib.pyplot as plt
from kinematic_functions import smooth_kinematic_function as sm_f
from kinematic_functions import sinusoidal_kinematic_function as si_f
from kinematics_write import write_2d, write_3d

time_series_length = 10000
start_time = 0
end_time = 80

t = np.linspace(start_time, end_time, time_series_length)

# kinematic parameters for default kinematic functions
# ----------------------------------------------
flapping_wing_frequency = 0.1
flapping_amplitude = 80
pitching_amplitude = 45

flapping_acceleration_time_coefficient = 0.99  # between 0 and 1
pitching_time_coefficient = 3.3  # between 0 and inf

flapping_delay_time_fraction = 0
pitching_delay_time_fraction = 0

section_location = 1

# additional kinematic control parameters for mostafa_function
# ----------------------------------------------
flapping_angular_velocity_amplitude = 0.5 * 180 / np.pi

flapping_acceleration_time_fraction = 0.24
pitching_time_fraction = 0.24
# ----------------------------------------------

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

kinematic_angles = sm_f(t, kinematic_parameters_smooth)
# kinematic_angles = si_f(t, kinematic_parameters_sinusoidal)

# plotting kinematic_angles
t_1st_cycle = [t1 for t1 in t if t1 <= 1 / flapping_wing_frequency]
plt.plot(t_1st_cycle, kinematic_angles)
plt.show()

# ----------------------------------------------
write_2d(t, section_location, flapping_wing_frequency, kinematic_angles)
write_3d(t, flapping_wing_frequency, kinematic_angles)
