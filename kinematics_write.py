"""script for writing kinematic functions to openfoam format"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R


def kf_plotter(t, kinematic_angles, time_series_length_per_cycle, legends):
    """
    A helper function to make a graph

    Parameters
    ----------
    t : array
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
    t = np.array(t)
    kinematic_angles = np.array(kinematic_angles)

    time_series_length = len(t)
    no_of_points_per_cycle = time_series_length_per_cycle

    fig, ax = plt.subplots(1, 1)

    y_arrays = np.zeros((time_series_length, 4))
    for i in range(time_series_length):
        i_moded = np.mod(i, no_of_points_per_cycle -1)

        y_arrays[i][0] = kinematic_angles[i_moded][0]
        y_arrays[i][1] = kinematic_angles[i_moded][1]
        y_arrays[i][2] = kinematic_angles[i_moded][2]
        y_arrays[i][3] = kinematic_angles[i_moded][3]

    for legend in legends:
        if legend == 'phi':
            ax.plot(t, y_arrays[:, 0], label='phi')
        elif legend == 'alf':
            ax.plot(t, y_arrays[:, 1], label='alf')
        elif legend == 'dphi':
            ax.plot(t, y_arrays[:, 2], label='dphi')
        elif legend == 'dalf':
            ax.plot(t, y_arrays[:, 3], label='dalf')

    ax.set_xlabel('t (seconds)')
    ax.set_ylabel('angle (degrees)')
    ax.set_title('kinematics plot')
    ax.legend()

    plt.show()

    return fig


def write_2d(t, section_location, time_series_length_per_cycle, kinematic_angles):
    """write kinematics data for 2d wing motion"""
    no_of_points_per_cycle = time_series_length_per_cycle
    time_series_length = len(t)
    kinematic_angles = np.array(kinematic_angles)

    initial_phi = kinematic_angles[0][0]
    for i in range(no_of_points_per_cycle):
        kinematic_angles[i][
            0] = kinematic_angles[i][0] - initial_phi

    t_disp = []
    r_angle = []
    for i in range(time_series_length):
        i_moded = np.mod(i, no_of_points_per_cycle -1)

        t_dispi = section_location * kinematic_angles[i_moded][0] * np.pi / 180
        t_dispi = [str(t_dispi), '0', '0']
        t_disp.append(t_dispi)

        # ----------------------------------------------
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


def write_3d(t, time_series_length_per_cycle, kinematic_angles):
    """write kinematics data for 3d wing motion"""
    no_of_points_per_cycle = time_series_length_per_cycle
    time_series_length = len(t)
    kinematic_angles = np.array(kinematic_angles)

    initial_phi = kinematic_angles[0][0]
    for i in range(no_of_points_per_cycle):
        kinematic_angles[i][
            0] = kinematic_angles[i][0] - initial_phi

    t_disp = []
    r_angle = []
    for i in range(time_series_length):
        i_moded = np.mod(i, no_of_points_per_cycle -1)

        t_dispi = ['0', '0', '0']
        t_disp.append(t_dispi)

        #----------------------------------------------------
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
