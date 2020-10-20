"""script for tabulated 6DOF revolving wing motion"""

import os
import numpy as np

from kinematic_functions import sinu_ramp_rev, smooth_linear_ramp
from kinematics_write import kf_plotter, write_2d, write_3d

# sinumation time definition and choose ramp functions to use
time_step_increment = 2e-3
start_time = 0
end_time = 14.5
steady_rotation_time = 1.62
#------------------------------
section_location = 1  #used only for 2d cases
#------------------------------
# ramp_function = 'sinu_ramp'
ramp_function = 'smooth_linear_ramp'
#------------------------------
#--parameter for sinusoidal ramp_function---
steady_rotation_frequency = 1
initial_ramp_time = 1
#-------------------------------------------
#-------------------------------------------
#--parameters for smooth_linear_ramp function---
ramp_stage_acceleration = 245.2
#-conner smoothing parameter, higher indicates shorter smooth range--
smooth_factor = 150
#--ramp time and initial zero velocity time--
ramp_time = 0.355
ramp_constant_time = 0.026
#----------------------------------------
#-decelleration process of wing motion--
ramp_mode = 'with_end_acc'
# ramp_mode = 'no_end_acc'
#----------------------------------------
#-wing pitching mode at decelleration phase--
pitch_mode = 'with_end_pitch'  #--used when ramp mode with_end_acc
# pitch_mode = 'no_end_pitch'
pitch_acceleration = 15035 / 4
pitch_time = 0.71
pitch_acc_time_fraction = 0.1  #--relative to pitch time: 0 ~ 1
pitch_delay_time_fraction = 0
#----------------------------------------
if ramp_function == 'smooth_linear_ramp':
    initial_ramp_time = ramp_time + ramp_constant_time

time_series_length = int(np.ceil(
    (end_time - start_time) / time_step_increment))
t = np.linspace(start_time, end_time, time_series_length)

i_ramp_end_time = start_time + initial_ramp_time
steady_end_time = i_ramp_end_time + steady_rotation_time
if ramp_function == 'smooth_linear_ramp':
    end_ramp_time = initial_ramp_time

    ramp_start_time = start_time + ramp_constant_time
    end_ramp_end_time = steady_end_time + end_ramp_time - ramp_constant_time

if ramp_function == 'sinu_ramp':
    kinematic_parameters = [steady_rotation_frequency, i_ramp_end_time]
    kinematic_angles = sinu_ramp_rev(t, kinematic_parameters)
elif ramp_function == 'smooth_linear_ramp':
    kinematic_parameters = [
        ramp_stage_acceleration, ramp_start_time, i_ramp_end_time,
        steady_end_time, end_ramp_end_time, smooth_factor, ramp_mode,
        ramp_constant_time, pitch_mode, pitch_time, pitch_delay_time_fraction,
        pitch_acceleration, pitch_acc_time_fraction, section_location
    ]
    kinematic_angles = smooth_linear_ramp(t, kinematic_parameters)

#--------------------------------------------------
angles_to_plot = ['dphi', 'dalf']

kf_plotter(t, kinematic_angles, angles_to_plot, 'basic', 'against_t',
           'current')

write_2d(t, section_location, kinematic_angles, 'basic', 'current')
write_3d(t, kinematic_angles, 'basic')
