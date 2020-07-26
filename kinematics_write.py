"""script for writing kinematic functions to openfoam format"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R


def kf_plotter(x_array, y_arrays, legends):
    """
    A helper function to make a graph

    Parameters
    ----------
    x_array : array
       The x data

    y_arrays : nd_array
       The y data arrays

    legends : list
       list of data name for each y array

    Returns
    -------
    out : list
        list of artists added
    """
    x_array = np.array(x_array)
    y_arrays = np.array(y_arrays)

    fig, ax = plt.subplots(1, 1)

    for legend in legends:
        if legend == 'phi':
            ax.plot(x_array, y_arrays[:, 0], label='phi')
        elif legend == 'alf':
            ax.plot(x_array, y_arrays[:, 1], label='alf')
        elif legend == 'dphi':
            ax.plot(x_array, y_arrays[:, 2], label='dphi')
        elif legend == 'dalf':
            ax.plot(x_array, y_arrays[:, 3], label='dalf')

    ax.set_xlabel('t (seconds)')
    ax.set_ylabel('angle (degrees)')
    ax.set_title('kinematics plot')
    ax.legend()

    plt.show()

    return fig


def write_2d(t, section_location, flapping_wing_frequency, kinematic_angles,
             write_mode):
    """write kinematics data for 2d wing motion"""
    t_1st_cycle = [t1 for t1 in t if t1 <= 1 / flapping_wing_frequency]
    no_of_points_per_cycle = len(t_1st_cycle)
    time_series_length = len(t)
    kinematic_angles = np.array(kinematic_angles)

    t_disp = []
    r_angle = []
    for i in range(time_series_length):
        i_moded = np.mod(i, no_of_points_per_cycle)

        t_dispi = section_location * kinematic_angles[i_moded][0] * np.pi / 180
        t_dispi = [str(t_dispi), '0', '0']
        t_disp.append(t_dispi)

        # ----------------------------------------------
        if write_mode == 'zero_start' and i <= np.where(
                kinematic_angles[:, 1] == np.amax(kinematic_angles[:, 1]))[0][0]:
            pitch_anglei = np.amax(kinematic_angles[:, 1])
        else:
            pitch_anglei = kinematic_angles[i_moded][1]

        kinematic_anglesi = [0, 0, pitch_anglei]

        roti = R.from_euler('YXZ', kinematic_anglesi, degrees=True)

        r_anglei = roti.as_euler('XYZ', degrees=True)
        r_anglei = [str(r_anglei[0]), str(r_anglei[1]), str(r_anglei[2])]

        r_angle.append(r_anglei)

    t = [str(ti) for ti in t]

    motion = [str(time_series_length), '(']
    for ti, disp_i, angle_i in zip(t, t_disp, r_angle):
        motioni = '(' + ti + ' ((' + ' '.join(disp_i) + ')' + '(' + ' '.join(
            angle_i) + ')))'

        motion.append(motioni)

    motion.append(')')

    with open('6DoF_2d.dat', 'w') as f:
        for item in motion:
            f.write("%s\n" % item)


def write_3d(t, flapping_wing_frequency, kinematic_angles, write_mode):
    """write kinematics data for 3d wing motion"""
    t_1st_cycle = [t1 for t1 in t if t1 <= 1 / flapping_wing_frequency]
    no_of_points_per_cycle = len(t_1st_cycle)
    time_series_length = len(t)
    kinematic_angles = np.array(kinematic_angles)

    t_disp = []
    r_angle = []
    for i in range(time_series_length):
        i_moded = np.mod(i, no_of_points_per_cycle)

        t_dispi = ['0', '0', '0']
        t_disp.append(t_dispi)

        #----------------------------------------------------
        if write_mode == 'zero_start' and i <= np.where(
                kinematic_angles[:, 1] == np.amax(kinematic_angles[:, 1]))[0][0]:
            pitch_anglei = np.amax(kinematic_angles[:, 1])
            # print(kinematic_angles[:, 1])
        else:
            pitch_anglei = kinematic_angles[i_moded][1]

        kinematic_anglesi = [kinematic_angles[i_moded][0], pitch_anglei, 0]

        roti = R.from_euler('YXZ', kinematic_anglesi, degrees=True)

        r_anglei = roti.as_euler('XYZ', degrees=True)
        r_anglei = [str(r_anglei[0]), str(r_anglei[1]), str(r_anglei[2])]

        r_angle.append(r_anglei)

    t = [str(ti) for ti in t]

    motion = [str(time_series_length), '(']
    for ti, disp_i, angle_i in zip(t, t_disp, r_angle):
        motioni = '(' + ti + ' ((' + ' '.join(disp_i) + ')' + '(' + ' '.join(
            angle_i) + ')))'

        motion.append(motioni)

    motion.append(')')

    with open('6DoF_3d.dat', 'w') as f:
        for item in motion:
            f.write("%s\n" % item)
