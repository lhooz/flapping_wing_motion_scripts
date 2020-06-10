"""script for tabulated 6DOF flapping wing motion"""

import numpy as np
from scipy.spatial.transform import Rotation as R

flapping_amplitude = 90
pitching_amplitude = 90
flapping_frequency = 0.1

time_series_length = 1000
start_time = 0
end_time = 100

t = np.linspace(start_time, end_time, time_series_length)

phi = (flapping_amplitude / 2) * np.sin(2 * np.pi * flapping_frequency * t -
                                        np.pi / 2)
alf = (pitching_amplitude / 2) * np.sin(2 * np.pi * flapping_frequency * t +
                                        np.pi)
# r = R.from_euler('zyx', [-30, -20, 0], degrees=True)
# r = r.as_rotvec()

r_vector = []
for phii, alfi in zip(phi, alf):
    euler_anglesi = [phii, alfi, 0]

    roti = R.from_euler('zyx', euler_anglesi, degrees=True)
    r_vectori = roti.as_rotvec() * 180 / np.pi
    r_vectori = [str(r_vectori[0]), str(r_vectori[1]), str(r_vectori[2])]

    r_vector.append(r_vectori)

t = [str(ti) for ti in t]

motion = ['1000', '(']
for ti, angle_i in zip(t, r_vector):
    motioni = '(' + ti + ' ((0 0 0)' + '(' + ' '.join(angle_i) + ')))'

    motion.append(motioni)

motion.append(')')

# print(motion)
with open('6DoF.dat', 'w') as f:
    for item in motion:
        f.write("%s\n" % item)
