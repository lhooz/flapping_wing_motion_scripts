"""script for tabulated 6DOF revolving wing motion"""

import numpy as np

from kinematic_functions import sinu_ramp_rev, smooth_linear_ramp
from kinematics_write import kf_plotter, write_2d, write_3d

# sinumation time definition and choose ramp functions to use
ramp_time_series_length = 100
steady_rotation_time_series_length = 1000
start_time = 0
initial_ramp_time = 1
steady_rotation_time = 3
#------------------------------
steady_rotation_frequency = 1  #--not used for smooth_linear_ramp func--
section_location = 1  #used only for 2d cases
#------------------------------
# ramp_function = 'sinu_ramp'
ramp_function = 'smooth_linear_ramp'
#------------------------------
#--additional parameters for smooth_linear_ramp function---
ramp_stage_acceleration = 400
#-conner smoothing parameter, higher indicates shorter smooth range--
smooth_factor = 50
#--ramp initial zero velocity time relative to total ramp time--
ramp_constant_v_length = 0.1
#-decelleration process of wing motion--
# ramp_mode = 'with_end_acc'
ramp_mode = 'no_end_acc'
#-zero velocity time after the wing stoped--
end_constant_time = 2.1
end_c_time_series_length = 2100
#------------------------------
i_ramp_end_time = start_time + initial_ramp_time
steady_end_time = i_ramp_end_time + steady_rotation_time
t_ramp = np.linspace(start_time, i_ramp_end_time, ramp_time_series_length)
t_rev = np.linspace(i_ramp_end_time, steady_end_time,
                    steady_rotation_time_series_length)
t_rev = np.delete(t_rev, 0)
t = np.append(t_ramp, t_rev)
if ramp_function == 'smooth_linear_ramp':
    end_ramp_time = initial_ramp_time
    ramp_constant_time = ramp_constant_v_length * initial_ramp_time

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
        ramp_constant_time
    ]
    kinematic_angles = smooth_linear_ramp(t, kinematic_parameters)
#--------------------------------------------------
angles_to_plot = ['dphi', 'ddphi']

kf_plotter(t, kinematic_angles, angles_to_plot, 'revolving', 'against_t')

write_2d(t, section_location, kinematic_angles, 'revolving')
write_3d(t, kinematic_angles, 'revolving')
