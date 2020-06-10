"""script for tabulated 6DOF flapping wing motion"""

import numpy as np

flapping_amplitude = 90
pitching_amplitude = 90
flapping_frequency = 0.1

time_series_length = 1000
start_time = 0
end_time = 100

t = np.linspace(start_time, end_time, time_series_length)

phi = (flapping_amplitude / 2) * np.sin(2 * np.pi * flapping_frequency * t -
                                        np.pi / 2)

alf = (pitching_amplitude / 2) * np.sin(2 * np.pi * flapping_frequency * t)

t = [str(ti) for ti in t]
angle_z = [str(phii) for phii in phi]
angle_y = [
    str(alfi * np.cos(phii * np.pi / 180)) for alfi, phii in zip(alf, phi)
]
angle_x = [
    str(alfi * np.sin(phii * np.pi / 180)) for alfi, phii in zip(alf, phi)
]
# print(phi)

motion = ['1000', '(']
for ti, angle_xi, angle_yi, angle_zi in zip(t, angle_x, angle_y, angle_z):
    motioni = '(' + ti + ' ((0 0 0)' + '(' + ' '.join(
        [angle_xi, angle_yi, angle_zi]) + '))'

    motion.append(motioni)

motion.append(')')

# print(motion)
with open('6DOF.dat', 'w') as f:
    for item in motion:
        f.write("%s\n" % item)
