"""script for 2d wake capture wing motion"""

import os
import shutil

import numpy as np

from kinematic_functions import (read_planning_parameters_csv,
                                 smooth_linear_ramp)
from kinematics_write import kf_plotter, write_2d

time_step_increment = 1e-3
parameters_file_name = '2d_case_parameters'
output_dir = '2d_kinematic_cases'
# sinumation time definition and choose ramp functions to use
ramp_function = 'smooth_linear_ramp'
ramp_mode = 'with_end_acc'
pitch_mode = 'with_end_pitch'  #--used when ramp mode with_end_acc
section_location = 1  #used only for 2d cases
start_time = 0
end_time = 2.6
#--------------------------------------------
ramp_constant_time = 0.02
pitch_acc_time_fraction = 0.5  #--relative to pitch time: 0 ~ 1
pitch_delay_time_fraction = 0
#-conner smoothing parameter, higher indicates shorter smooth range--
smooth_factor = 100
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
        case[2]) + '_pf' + str(case[3]) + '_pa' + str(case[4])
    save_file_data = os.path.join(output_dir_path, file_name + '.dat')
    save_file_image = os.path.join(output_dir_path, file_name + '.png')
    save_file_cf = os.path.join(output_dir_path, file_name + '.cf')
    save_file_nu = os.path.join(output_dir_path, file_name + '.nu')

    #--------------------------------------------
    #--ramp time and initial zero velocity time--
    ramp_time = case[5]
    steady_rotation_time = case[6]
    pitch_time = case[7]
    nu = case[8]
    ramp_stage_acceleration = case[9] / section_location * 180 / np.pi
    pitch_acceleration = case[10]
    ref_area = 1
    ref_vel = case[11]
    #-------------------------------------------
    if ramp_function == 'smooth_linear_ramp':
        initial_ramp_time = ramp_time + ramp_constant_time
        # end_constant_time = end_time - 2 * initial_ramp_time - steady_rotation_time

    time_series_length = int(
        np.ceil((end_time - start_time) / time_step_increment))
    t = np.linspace(start_time, end_time, time_series_length)

    i_ramp_end_time = start_time + initial_ramp_time
    steady_end_time = i_ramp_end_time + steady_rotation_time
    if ramp_function == 'smooth_linear_ramp':
        end_ramp_time = initial_ramp_time

        ramp_start_time = start_time + ramp_constant_time
        end_ramp_end_time = steady_end_time + end_ramp_time - ramp_constant_time

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
    angles_to_plot = ['dphi', 'dalf']
    kf_plotter(t, kinematic_angles, angles_to_plot, 'basic', 'against_t',
               save_file_image)
    write_2d(t, section_location, kinematic_angles, 'basic', save_file_data)

    with open(save_file_cf, 'w') as f:
        f.write(
            '%s\n' %
            '/*--------------------------------*- C++ -*----------------------------------*\\'
        )
        f.write(
            '%s\n' %
            r'\*---------------------------------------------------------------------------*/'
        )
        f.write('%s\n%s\n' % (r'forceCoeffs_object', r'{'))
        f.write('%s\n%s\n' %
                (r'    type forceCoeffs;', r'    libs ("libforces.so");'))
        f.write('%s\n\n%s\n\n' %
                (r'    patches (wing);', r'    writeControl writeTime;'))
        f.write('%s\n%s\n\n' % (r'    pName p;', r'    Uname U;'))
        f.write('%s\n%s\n\n' % (r'    rho rhoInf;', r'    rhoInf 1;'))
        f.write('%s\n\n%s\n' % (r'    log true;', r'    CofR (0 0 0);'))
        f.write('%s\n%s\n' %
                (r'    liftDir (0 1 0);', r'    dragDir (1 0 0);'))
        f.write('%s\n' % (r'    pitchAxis (0 0 1);'))
        f.write('    magUInf %s%s\n' % (str(ref_vel), r';'))
        f.write('    lRef 1;\n')
        f.write('    Aref %s%s\n' % (str(ref_area), r';'))
        f.write(r'}')

    with open(save_file_nu, 'w') as f:
        f.write(
            '%s\n' %
            '/*--------------------------------*- C++ -*----------------------------------*\\'
        )
        f.write(
            '%s\n' %
            r'\*---------------------------------------------------------------------------*/'
        )
        f.write('%s\n%s\n' % (r'FoamFile', r'{'))
        f.write('    version     2.0;\n')
        f.write('    format      ascii;\n')
        f.write('    class       dictionary;\n')
        f.write('    object      transportProperties;\n')
        f.write('}\n')
        f.write(
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n'
        )
        f.write('transportModel  Newtonian;\n\n')
        f.write('nu              nu [ 0 2 -1 0 0 0 0 ] %s%s\n\n' %
                (str(nu), r';'))
        f.write(
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n'
        )
