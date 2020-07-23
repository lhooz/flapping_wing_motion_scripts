"""script for writing kinematic functions to openfoam format"""

import numpy as np
from scipy.spatial.transform import Rotation as R

def write_2d(t, section_location, flapping_wing_frequency, kinematic_angles):
    """write kinematics data for 2d wing motion"""
    t_1st_cycle = [t1 for t1 in t if t1 <= 1/flapping_wing_frequency]
    no_of_points_per_cycle = len(t_1st_cycle)
    time_series_length = len(t)


    t_disp = []
    r_angle = []
    for i in range(time_series_length):
        i_moded = np.mod(i, no_of_points_per_cycle)

        t_dispi = section_location * kinematic_angles[i_moded][0] * np.pi / 180
        t_dispi = [str(t_dispi), '0', '0']

        t_disp.append(t_dispi)

        # ----------------------------------------------
        kinematic_anglesi = [0, 0, kinematic_angles[i_moded][1]]

        roti = R.from_euler('YXZ', kinematic_anglesi, degrees=True)

        r_anglei = roti.as_euler('XYZ', degrees=True)
        r_anglei = [str(r_anglei[0]), str(r_anglei[1]), str(r_anglei[2])]

        r_angle.append(r_anglei)

    t = [str(ti) for ti in t]

    motion = [str(time_series_length), '(']
    for ti, disp_i, angle_i in zip(t, t_disp, r_angle):
        motioni = '(' + ti + ' (('+ ' '.join(disp_i)+')' + '(' + ' '.join(angle_i) + ')))'

        motion.append(motioni)

    motion.append(')')

    with open('6DoF_2d.dat', 'w') as f:
        for item in motion:
            f.write("%s\n" % item)


def write_3d(t, flapping_wing_frequency, kinematic_angles):
    """write kinematics data for 3d wing motion"""
    t_1st_cycle = [t1 for t1 in t if t1 <= 1/flapping_wing_frequency]
    no_of_points_per_cycle = len(t_1st_cycle)
    time_series_length = len(t)


    r_angle = []
    for i in range(time_series_length):
        i_moded = np.mod(i, no_of_points_per_cycle)
        kinematic_anglesi = [kinematic_angles[i_moded][0], kinematic_angles[i_moded][1], 0]

        roti = R.from_euler('YXZ', kinematic_anglesi, degrees=True)

        r_anglei = roti.as_euler('XYZ', degrees=True)
        r_anglei = [str(r_anglei[0]), str(r_anglei[1]), str(r_anglei[2])]

        r_angle.append(r_anglei)

    t = [str(ti) for ti in t]

    motion = [str(time_series_length), '(']
    for ti, angle_i in zip(t, r_angle):
        motioni = '(' + ti + ' ((0 0 0)' + '(' + ' '.join(angle_i) + ')))'

        motion.append(motioni)

    motion.append(')')

    with open('6DoF_3d.dat', 'w') as f:
        for item in motion:
            f.write("%s\n" % item)
