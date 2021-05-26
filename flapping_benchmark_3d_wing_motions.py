"""scripts for generating 3d benchmark flapping wing motions"""

import os
import shutil

import numpy as np

from kinematic_functions import read_planning_parameters_csv
from kinematic_functions import sinu_continuous_kinematic_function as sic_f
from kinematic_functions import smooth_linear_ramp
from kinematics_write import kf_plotter, write_3d

parameters_file_name = 'wplanform_data_re100'
output_dir = '3dbm_kinematic_cases'
# sinumation time definition and choose ramp functions to use
pitching_wave_form = [0.25]
time_series_length_per_cycle = 1501
start_time = 0
number_of_cycles = 6
use_function = 'sinu_continuous'
flapping_delay_time_fraction = 0
pitching_delay_time_fraction = 0
flapping_acceleration_time_fraction = 0.25
#----------------------------------------
cwd = os.getcwd()
output_dir_path = os.path.join(cwd, output_dir)
if os.path.exists(output_dir_path):
    shutil.rmtree(output_dir_path)
os.mkdir(output_dir_path)
parameters_file = os.path.join(cwd, parameters_file_name + '.csv')
parameters_arr = read_planning_parameters_csv(parameters_file)

for pitch_tf in pitching_wave_form:
    pitching_time_fraction = pitch_tf
    for case in parameters_arr:
        Ro = case[5] * (case[1] + case[2])
        file_name = 'ar' + str(case[1]) + '_ofs' + str(case[2]) + '_r1h' + str(
            case[3]) + '__Re' + str(
                case[0]) + '_pt' + str(pitch_tf) + '_Ro' + '{0:.2f}'.format(Ro)
        save_file_data = os.path.join(output_dir_path, file_name + '.dat')
        save_file_image = os.path.join(output_dir_path, file_name + '.png')
        save_file_cf = os.path.join(output_dir_path, file_name + '.cf')
        save_file_nu = os.path.join(output_dir_path, file_name + '.nu')

        #------------------------------------
        flapping_wing_frequency = case[9]
        nu = case[10]
        ref_area = case[11]
        ref_vel = case[12]

        flapping_angular_velocity_amplitude = 391.06 * flapping_wing_frequency  # --degree/s--
        pitching_angular_velocity_amplitude = 360 * flapping_wing_frequency / (
            2 * pitching_time_fraction)  # --degree/s--
        t1 = np.linspace(start_time, 1 / flapping_wing_frequency,
                         time_series_length_per_cycle)
        t = t1
        for i in range(1, number_of_cycles):
            ti = np.delete(t1, 0) + i / flapping_wing_frequency
            t = np.append(t, ti)

        kinematic_parameters_sinu_continuous = [
            flapping_wing_frequency, flapping_angular_velocity_amplitude,
            pitching_angular_velocity_amplitude,
            flapping_acceleration_time_fraction, pitching_time_fraction,
            flapping_delay_time_fraction, pitching_delay_time_fraction
        ]

        kinematic_angles = sic_f(t, kinematic_parameters_sinu_continuous)
        # plotting kinematic angles
        angles_to_plot = ['phi', 'dphi', 'alf', 'dalf']
        # angles_to_plot = ['dphi', 'dalf', 'ddphi']

        kf_plotter(t, kinematic_angles, angles_to_plot,
                   time_series_length_per_cycle, 'against_t', save_file_image)
        # ----------------------------------------------
        write_3d(t, kinematic_angles, time_series_length_per_cycle,
                 save_file_data)

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
                    (r'    liftDir (1 0 0);', r'    dragDir (0 0 1);'))
            f.write('%s\n' % (r'    pitchAxis (0 1 0);'))
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
