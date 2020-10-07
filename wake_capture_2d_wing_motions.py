"""script for tabulated 6DOF revolving wing motion"""

import os
import shutil

import numpy as np

from kinematic_functions import (read_planning_parameters_csv,
                                 smooth_linear_ramp)
from kinematics_write import kf_plotter, write_2d

time_step_increment = 2e-3
parameters_file_name = '2d_case_parameters'
output_dir = '2d_kinematic_cases'
# sinumation time definition and choose ramp functions to use
ramp_function = 'smooth_linear_ramp'
ramp_mode = 'with_end_acc'
pitch_mode = 'with_end_pitch'  #--used when ramp mode with_end_acc
section_location = 1  #used only for 2d cases
start_time = 0
end_time = 14.5
#--------------------------------------------
ramp_constant_time = 0.1
pitch_acc_time_fraction = 0.1  #--relative to pitch time: 0 ~ 1
pitch_delay_time_fraction = 0
#-conner smoothing parameter, higher indicates shorter smooth range--
smooth_factor = 50
#------------------------------------------
cwd = os.getcwd()
output_dir_path = os.path.join(cwd, output_dir)
if os.path.exists(output_dir_path):
    shutil.rmtree(output_dir_path)
os.mkdir(output_dir_path)
parameters_file = os.path.join(cwd, parameters_file_name + '.csv')
parameters_arr = read_planning_parameters_csv(parameters_file)
for case in parameters_arr:
    file_name = 'Re' + str(case[0]) + '_stroke' + str(case[1]) + '_acf' + str(
        case[2]) + '_pf' + str(case[3])
    save_file_data = os.path.join(output_dir_path, file_name + '.dat')
    save_file_image = os.path.join(output_dir_path, file_name + '.png')
    #--------------------------------------------
    #--ramp time and initial zero velocity time--
    ramp_time = case[4]
    steady_rotation_time = case[5]
    pitch_time = case[6]
    ramp_stage_acceleration = case[8] / section_location * 180 / np.pi
    pitch_acceleration = case[9]
    #-------------------------------------------
    if ramp_function == 'smooth_linear_ramp':
        initial_ramp_time = ramp_time + ramp_constant_time
        end_constant_time = end_time - 2 * initial_ramp_time - steady_rotation_time

    ramp_time_series_length = int(
        np.ceil(initial_ramp_time / time_step_increment))
    steady_rotation_time_series_length = int(
        np.ceil(steady_rotation_time / time_step_increment))
    end_c_time_series_length = int(
        np.ceil(end_constant_time / time_step_increment))

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

    if ramp_function == 'smooth_linear_ramp':
        kinematic_parameters = [
            ramp_stage_acceleration, ramp_start_time, i_ramp_end_time,
            steady_end_time, end_ramp_end_time, smooth_factor, ramp_mode,
            ramp_constant_time, pitch_mode, pitch_time,
            pitch_delay_time_fraction, pitch_acceleration,
            pitch_acc_time_fraction, section_location
        ]
        kinematic_angles = smooth_linear_ramp(t, kinematic_parameters)

    #--------------------------------------------------
    angles_to_plot = ['dphi', 'ddphi']
    kf_plotter(t, kinematic_angles, angles_to_plot, 'basic', 'against_t',
               save_file_image)
    write_2d(t, section_location, kinematic_angles, 'basic', save_file_data)
