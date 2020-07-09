"""script for tabulated 6DOF flapping wing motion"""

import numpy as np
from scipy.spatial.transform import Rotation as R

flapping_amplitude = 120
pitching_amplitude = 90
flapping_frequency = 0.1

time_series_length = 1000
start_time = 0
end_time = 100

t = np.linspace(start_time, end_time, time_series_length)


def phi(x):
    """flapping motion function"""
    return (flapping_amplitude /
            2) * np.sin(2 * np.pi * flapping_frequency * x - np.pi / 2)


def alf(x):
    """pitching motion function"""
    return (pitching_amplitude /
            2) * np.sin(2 * np.pi * flapping_frequency * x + np.pi)

# r = R.from_euler('zyx', [-30, -20, 0], degrees=True)
# r = r.as_rotvec()
# test_vector = [0, 1, 0]
r_angle = []
for ti in t:
    euler_anglesi = [phi(ti), alf(ti), 0]

    roti = R.from_euler('ZYX', euler_anglesi, degrees=True)
    # result_vectori = roti.apply(test_vector)
    # print(result_vectori[2])

    r_anglei = roti.as_euler('XYZ', degrees=True)
    r_anglei = [str(r_anglei[0]), str(r_anglei[1]), str(r_anglei[2])]

    r_angle.append(r_anglei)

t = [str(ti) for ti in t]

motion = ['1000', '(']
for ti, angle_i in zip(t, r_angle):
    motioni = '(' + ti + ' ((0 0 0)' + '(' + ' '.join(angle_i) + ')))'

    motion.append(motioni)

motion.append(')')

# print(motion)
with open('6DoF.dat', 'w') as f:
    for item in motion:
        f.write("%s\n" % item)
