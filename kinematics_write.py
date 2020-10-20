"""script for writing kinematic functions to openfoam format"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R


def kf_plotter(t, kinematic_angles, legends, time_series_length_per_cycle,
               h_axis, save_file):
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

    time_series_length_per_cycle: number or 'basic'
        flapping wing cases use number for length of time series per flapping cycle
        basic kinematic cases use 'basic' keyword to indicate

    h_axis: 'against_t' or 'against_phi'
        'against_t' for use time as horizontal axis
        'against_phi' for use stroke angle as horizontal axis

    """
    cwd = os.getcwd()
    t = np.array(t)
    kinematic_angles = np.array(kinematic_angles)
    l_width = 1.5

    time_series_length = len(t)
    if time_series_length_per_cycle == 'basic':
        kinematic_angles = kinematic_angles * -1
        no_of_points_per_cycle = time_series_length + 1
    else:
        no_of_points_per_cycle = time_series_length_per_cycle

    fig, ax = plt.subplots(1, 1)

    y_arrays = np.zeros((time_series_length, 6))
    for i in range(time_series_length):
        i_moded = np.mod(i, no_of_points_per_cycle - 1)

        y_arrays[i][0] = kinematic_angles[i_moded][0]
        y_arrays[i][1] = kinematic_angles[i_moded][1]
        y_arrays[i][2] = kinematic_angles[i_moded][2]
        y_arrays[i][3] = kinematic_angles[i_moded][3]
        y_arrays[i][4] = kinematic_angles[i_moded][4]
        y_arrays[i][5] = kinematic_angles[i_moded][5]

    y_label = []
    if h_axis == 'against_t':
        for legend in legends:
            if legend == 'phi':
                ax.plot(t,
                        y_arrays[:, 0],
                        linestyle='solid',
                        linewidth=l_width,
                        label=r'$\phi$')
                y_label.append(r'\phi')
            elif legend == 'alf':
                ax.plot(t,
                        y_arrays[:, 1],
                        linestyle='solid',
                        linewidth=l_width,
                        label=r'$\alpha$')
                y_label.append(r'\alpha')
            elif legend == 'dphi':
                ax.plot(t,
                        y_arrays[:, 2],
                        linestyle='dashed',
                        linewidth=l_width,
                        label=r'$\dot\phi$')
                y_label.append(r'\dot\phi')
            elif legend == 'dalf':
                ax.plot(t,
                        y_arrays[:, 3],
                        linestyle='dashed',
                        linewidth=l_width,
                        label=r'$\dot\alpha$')
                y_label.append(r'\dot\alpha')
            elif legend == 'ddphi':
                ax.plot(t,
                        y_arrays[:, 4],
                        linestyle='dashdot',
                        linewidth=l_width,
                        label=r'$\ddot\phi$')
                y_label.append(r'\ddot\phi')
            elif legend == 'ddalf':
                ax.plot(t,
                        y_arrays[:, 5],
                        linestyle='dashdot',
                        linewidth=l_width,
                        label=r'$\ddot\alpha$')
                y_label.append(r'\ddot\alpha')

        ax.set_xlabel('t (s)')

    elif h_axis == 'against_phi':
        x = y_arrays[:, 0]
        for legend in legends:
            if legend == 'phi':
                ax.plot(x,
                        y_arrays[:, 0],
                        linestyle='solid',
                        linewidth=l_width,
                        label=r'$\phi$')
                y_label.append(r'\phi')
            elif legend == 'alf':
                ax.plot(x,
                        y_arrays[:, 1],
                        linestyle='solid',
                        linewidth=l_width,
                        label=r'$\alpha$')
                y_label.append(r'\alpha')
            elif legend == 'dphi':
                ax.plot(x,
                        y_arrays[:, 2],
                        linestyle='dashed',
                        linewidth=l_width,
                        label=r'$\dot\phi$')
                y_label.append(r'\dot\phi')
            elif legend == 'dalf':
                ax.plot(x,
                        y_arrays[:, 3],
                        linestyle='dashed',
                        linewidth=l_width,
                        label=r'$\dot\alpha$')
                y_label.append(r'\dot\alpha')
            elif legend == 'ddphi':
                ax.plot(x,
                        y_arrays[:, 4],
                        linestyle='dashdot',
                        linewidth=l_width,
                        label=r'$\ddot\phi$')
                y_label.append(r'\ddot\phi')
            elif legend == 'ddalf':
                ax.plot(x,
                        y_arrays[:, 5],
                        linestyle='dashdot',
                        linewidth=l_width,
                        label=r'$\ddot\alpha$')
                y_label.append(r'\ddot\alpha')

        ax.set_xlabel(r'$\phi\/(\deg)$')

    y_label_units = []
    if 'phi' in legends or 'alf' in legends:
        y_label_units.append(r'\deg')
    if 'dphi' in legends or 'dalf' in legends:
        y_label_units.append(r'\deg/s')
    if 'ddphi' in legends or 'ddalf' in legends:
        y_label_units.append(r'\deg/s^2')

    y_label_str = r'$' + ',\/'.join(y_label) + r'\/(' + ',\/'.join(
        y_label_units) + r')' + r'$'
    ax.set_ylabel(y_label_str)
    ax.set_title('kinematics plot')
    ax.legend()

    if save_file == 'current':
        out_figure_file = os.path.join(cwd, 'kinematics_plot.png')
        fig.savefig(out_figure_file)
        plt.show()
    else:
        fig.savefig(save_file)
        plt.close()

    return fig


def write_2d(t, section_location, kinematic_angles,
             time_series_length_per_cycle, save_file):
    """write kinematics data for 2d wing motion"""
    kinematic_angles = np.array(kinematic_angles)
    time_series_length = len(t)

    if time_series_length_per_cycle == 'basic':
        no_of_points_per_cycle = time_series_length + 1
    else:
        no_of_points_per_cycle = time_series_length_per_cycle

        initial_phi = kinematic_angles[0][0]
        initial_alf = kinematic_angles[0][1]
        for i in range(no_of_points_per_cycle):
            kinematic_angles[i][0] = kinematic_angles[i][0] - initial_phi
            kinematic_angles[i][1] = kinematic_angles[i][1] - initial_alf

    t_disp = []
    r_angle = []
    for i in range(time_series_length):
        i_moded = np.mod(i, no_of_points_per_cycle - 1)

        t_dispi = section_location * kinematic_angles[i_moded][0] * np.pi / 180
        t_dispi = ['{:0.08g}'.format(t_dispi), '0', '0']
        t_disp.append(t_dispi)

        # ----------------------------------------------
        pitch_anglei = kinematic_angles[i_moded][1]

        kinematic_anglesi = [0, 0, pitch_anglei]  #--pitch about z axis--

        # roti = R.from_euler('XYZ', kinematic_anglesi, degrees=True)
        # r_anglei = roti.as_euler('XYZ', degrees=True)
        r_anglei = kinematic_anglesi
        r_anglei = [
            '{:0.08g}'.format(r_anglei[0]), '{:0.08g}'.format(r_anglei[1]),
            '{:0.08g}'.format(r_anglei[2])
        ]
        r_angle.append(r_anglei)

    t = ['{:0.08g}'.format(ti) for ti in t]
    motion = [str(time_series_length), '(']
    for ti, disp_i, angle_i in zip(t, t_disp, r_angle):
        motioni = '(' + ti + ' ((' + ' '.join(disp_i) + ')' + '(' + ' '.join(
            angle_i) + ')))'

        motion.append(motioni)

    motion.append(')')

    if save_file == 'current':
        with open('6DoF_2d.dat', 'w') as f:
            for item in motion:
                f.write("%s\n" % item)
    else:
        with open(save_file, 'w') as f:
            for item in motion:
                f.write("%s\n" % item)


def write_3d(t, kinematic_angles, time_series_length_per_cycle, save_file):
    """write kinematics data for 3d wing motion"""
    kinematic_angles = np.array(kinematic_angles)
    time_series_length = len(t)

    if time_series_length_per_cycle == 'basic':
        no_of_points_per_cycle = time_series_length + 1
    else:
        no_of_points_per_cycle = time_series_length_per_cycle

        initial_phi = kinematic_angles[0][0]
        initial_alf = kinematic_angles[0][1]
        for i in range(no_of_points_per_cycle):
            kinematic_angles[i][0] = kinematic_angles[i][0] - initial_phi
            kinematic_angles[i][1] = kinematic_angles[i][1] - initial_alf

    t_disp = []
    r_angle = []
    for i in range(time_series_length):
        i_moded = np.mod(i, no_of_points_per_cycle - 1)

        t_dispi = ['0', '0', '0']
        t_disp.append(t_dispi)

        #----------------------------------------------------
        pitch_anglei = kinematic_angles[i_moded][1]

        kinematic_anglesi = [kinematic_angles[i_moded][0], pitch_anglei, 0]

        # roti = R.from_euler('XYZ', kinematic_anglesi, degrees=True)
        # r_anglei = roti.as_euler('XYZ', degrees=True)
        r_anglei = kinematic_anglesi
        r_anglei = [
            '{:0.08g}'.format(r_anglei[0]), '{:0.08g}'.format(r_anglei[1]),
            '{:0.08g}'.format(r_anglei[2])
        ]

        r_angle.append(r_anglei)

    t = ['{:0.08g}'.format(ti) for ti in t]
    motion = [str(time_series_length), '(']
    for ti, disp_i, angle_i in zip(t, t_disp, r_angle):
        motioni = '(' + ti + ' ((' + ' '.join(disp_i) + ')' + '(' + ' '.join(
            angle_i) + ')))'

        motion.append(motioni)

    motion.append(')')

    if save_file == 'current':
        with open('6DoF_3d.dat', 'w') as f:
            for item in motion:
                f.write("%s\n" % item)
    else:
        with open(save_file, 'w') as f:
            for item in motion:
                f.write("%s\n" % item)


def write_iaoa(kinematic_angles):
    """write pitch angle at 1st step"""

    with open('iaoa.dat', 'w') as f:
        f.write("%s\n" % kinematic_angles[0][1])


def write_max_dphi(kinematic_angles):
    """write maximum flapping angular velocity"""

    kinematic_angles = np.array(kinematic_angles)
    mdphi = np.amax(kinematic_angles[:, 2])

    with open('max_dphi.dat', 'w') as f:
        f.write("%s\n" % mdphi)
