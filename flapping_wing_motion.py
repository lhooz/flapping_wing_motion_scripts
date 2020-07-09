"""script for tabulated 6DOF flapping wing motion"""

import autograd.numpy as np
from scipy.spatial.transform import Rotation as R
from autograd import elementwise_grad as grad
import matplotlib.pyplot as plt

flapping_amplitude = 120
pitching_amplitude = 90
flapping_frequency = 0.1

time_series_length = 1000
start_time = 0
end_time = 100

t = np.linspace(start_time, end_time, time_series_length)


def phi(x):
    """flapping angle function"""
    return (flapping_amplitude / 2) * np.sin(
        2 * np.pi * flapping_frequency * x)


def alf(x):
    """pitching angle function"""
    return (pitching_amplitude /
            2) * np.sin(2 * np.pi * flapping_frequency * x + np.pi)


dphi = grad(phi)
dalf = grad(alf)

# plt.plot(t, dphi(t), t, dalf(t))
# plt.show()

omega = []
# r_vector = []
for ti in t:
    euler_anglesi = [phi(ti), alf(ti), 0]

    roti = R.from_euler('ZYX', euler_anglesi, degrees=True)

    dz = np.array([0, 0, dphi(ti)])
    # dy = roti.apply([0, dalf(ti), 0])
    dy = roti.apply(np.array([0, dalf(ti), 0]))
    omegai = dz + dy

    omega.append(omegai * np.pi / 180)

    # r_vectori = roti.as_rotvec() * 180 / np.pi
    # r_vectori = [str(r_vectori[0]), str(r_vectori[1]), str(r_vectori[2])]

    # r_vector.append(r_vectori)

t = [str(ti) for ti in t]
omegax = []
omegay = []
omegaz = []
for ti, omegai in zip(t, omega):
    omegaxi = '(' + ti + ' ' + str(omegai[0] * 2.05) + ')'
    omegayi = '(' + ti + ' ' + str(omegai[1]) + ')'
    omegazi = '(' + ti + ' ' + str(omegai[2]) + ')'

    omegax.append(omegaxi)
    omegay.append(omegayi)
    omegaz.append(omegazi)
# motion = ['1000', '(']
# for ti, angle_i in zip(t, r_vector):
# motioni = '(' + ti + ' ((0 0 0)' + '(' + ' '.join(angle_i) + ')))'

# motion.append(motioni)

# motion.append(')')

with open('omegax.dat', 'w') as f:
    for item in omegax:
        f.write("%s\n" % item)

with open('omegay.dat', 'w') as f:
    for item in omegay:
        f.write("%s\n" % item)

with open('omegaz.dat', 'w') as f:
    for item in omegaz:
        f.write("%s\n" % item)
