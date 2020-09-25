"""script for tabulated 6DOF revolving wing motion"""

import numpy as np

from kinematic_functions import sinu_ramp_rev, smooth_linear_ramp
from kinematics_write import kf_plotter, write_2d, write_3d

# sinumation time definition and choose ramp functions to use
ramp_time_series_length = 120
steady_rotation_time_series_length = 1000
start_time = 0
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
pitch_mode = 'with_end_pitch' #--used when ramp mode with_end_acc
# pitch_mode = 'no_end_pitch'
pitch_acceleration = 7935
pitch_time = 0.355
pitch_acc_time_fraction = 0.2  #--relative to pitch time: 0 ~ 1
pitch_delay_time_fraction = 0
#----------------------------------------
#-zero velocity time after the wing stoped--
end_constant_time = 4.764
end_c_time_series_length = 2000
#------------------------------
#------------------------------
if ramp_function == 'smooth_linear_ramp':
    initial_ramp_time = ramp_time + ramp_constant_time

i_ramp_end_time = start_time + initial_ramp_time
steady_end_time = i_ramp_end_time + steady_rotation_time
t_ramp = np.linspace(start_time, i_ramp_end_time, ramp_time_series_length)
t_rev = np.linspace(i_ramp_end_time, steady_end_time,
                    steady_rotation_time_series_length)
t_rev = np.delete(t_rev, 0)
t = np.append(t_ramp, t_rev)
if ramp_function == 'smooth_linear_ramp':
    end_ramp_time = initial_ramp_time

    ramp_start_time = start_time + ramp_constant_time
    end_ramp_end_time = steady_end_time + end_ramp_time - ramp_constant_time

    if ramp_mode == 'with_end_acc':
        t_end_ramp = np.linspace(steady_end_time,
                                 steady_end_time + end_ramp_time,
                                 ramp_time_series_length)
        t_end_ramp = np.delete(t_end_ramp, 0)
        t = np.append(t, t_end_ramp)

        final_end_time = steady_end_time + end_ramp_time + end_constant_time
        t_end_final = np.linspace(steady_end_time + end_ramp_time,
                                  final_end_time, end_c_time_series_length)

        t_end_final = np.delete(t_end_final, 0)
        t = np.append(t, t_end_final)

time_series_length = len(t)
# print(t)

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

kf_plotter(t, kinematic_angles, angles_to_plot, 'basic', 'against_t')

write_2d(t, section_location, kinematic_angles, 'basic')
write_3d(t, kinematic_angles, 'basic')
