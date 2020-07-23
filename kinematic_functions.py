"""kinematic functions definitions: default (smooth) and mostafa (sinusiodal)"""

import autograd.numpy as np
import scipy.integrate as integrate
from autograd import grad
from scipy.misc import derivative


# default kinematic functions for phi and alf:
# -----------------------------------------------
def smooth_kinematic_function(t, kinematic_parameters):
    """definition for default (smooth) kinematic function"""
    flapping_wing_frequency = kinematic_parameters[0]
    flapping_amplitude = kinematic_parameters[1]
    pitching_amplitude = kinematic_parameters[2]
    flapping_acceleration_time_coefficient = kinematic_parameters[3]
    pitching_time_coefficient = kinematic_parameters[4]
    flapping_delay_time_fraction = kinematic_parameters[5]
    pitching_delay_time_fraction = kinematic_parameters[6]

    def phi(x):
        """flapping motion function"""
        if flapping_acceleration_time_coefficient == 0:
            phit = flapping_amplitude * np.cos(
                2 * np.pi * flapping_wing_frequency *
                (x - flapping_delay_time_fraction / flapping_wing_frequency))
        else:
            phit = flapping_amplitude / np.arcsin(
                flapping_acceleration_time_coefficient
            ) * np.arcsin(flapping_acceleration_time_coefficient * np.cos(
                2 * np.pi * flapping_wing_frequency *
                (x - flapping_delay_time_fraction / flapping_wing_frequency)))
        return phit

    def dphi(x):
        """flapping angular velocity function"""
        return grad(phi)(x)

    def alf(x):
        """pitching motion function"""
        if pitching_time_coefficient == 0:
            alft = pitching_amplitude * np.sin(
                2 * np.pi * flapping_wing_frequency *
                (x - pitching_delay_time_fraction / flapping_wing_frequency))
        else:
            alft = pitching_amplitude / np.tanh(
                pitching_time_coefficient
            ) * np.tanh(pitching_time_coefficient * np.sin(
                2 * np.pi * flapping_wing_frequency *
                (x - pitching_delay_time_fraction / flapping_wing_frequency)))
        return alft

    def dalf(x):
        """flapping angular velocity function"""
        return grad(alf)(x)

    kinematic_angles = []
    t_1st_cycle = [t1 for t1 in t if t1 <= 1 / flapping_wing_frequency]
    for ti in t_1st_cycle:
        kinematic_anglesi = [phi(ti), alf(ti), dphi(ti), dalf(ti)]
        kinematic_angles.append(kinematic_anglesi)

    return kinematic_angles


# -------------------------------------------------


# mostafa kinematic functions
# -------------------------------------------------
def sinusoidal_kinematic_function(t, kinematic_parameters):
    """definition for mostafa (sinusiodal) kinematic functiontion"""
    flapping_wing_frequency = kinematic_parameters[0]
    flapping_angular_velocity_amplitude = kinematic_parameters[1]
    pitching_amplitude = kinematic_parameters[2]
    flapping_acceleration_time_fraction = kinematic_parameters[3]
    pitching_time_fraction = kinematic_parameters[4]
    flapping_delay_time_fraction = kinematic_parameters[5]
    pitching_delay_time_fraction = kinematic_parameters[6]

    def dphi(x):
        """flapping motion angular velocity function"""
        return -kf(
            flapping_wing_frequency, flapping_angular_velocity_amplitude,
            flapping_acceleration_time_fraction, flapping_delay_time_fraction,
            x)

    flapping_amplitude = integrate.quad(lambda x: np.abs(dphi(x)), 0,
                                        1 / flapping_wing_frequency)[0]

    dphi_int = []
    for ti in t:
        dphi_int.append(dphi(ti))

    def phi(x):
        """flapping motion function"""
        t_int = [tx for tx in t if tx <= x]
        time_array_length = len(t_int)
        # print(time_array_length, dphi_int[0:time_array_length])
        return flapping_amplitude / 4 + integrate.simps(
            dphi_int[0:time_array_length], t_int)

    def alf(x):
        """pitching motion function"""
        return kf(flapping_wing_frequency, pitching_amplitude,
                  pitching_time_fraction, pitching_delay_time_fraction, x)

    def dalf(x):
        """flapping angular velocity function"""
        return derivative(alf, x, dx=1e-6)

    kinematic_angles = []
    t_1st_cycle = [t1 for t1 in t if t1 <= 1 / flapping_wing_frequency]
    for ti in t_1st_cycle:
        kinematic_anglesi = [phi(ti), alf(ti), dphi(ti), dalf(ti)]
        kinematic_angles.append(kinematic_anglesi)

    return kinematic_angles


    # original mostafa (sinusiodal) function
def kf(f, amp, tr_fac, del_fac, t):
    """kinematic_function for flapping velocity and angle of attack"""
    T = 1 / f
    delay = del_fac * T
    tr = tr_fac * T
    beta = 1 - (2 * tr / T)

    t_T1 = 0
    t_T2 = (T * (1 - beta) / 4) + delay
    t_T3 = ((T * (1 + beta) / 4)) + delay
    t_T4 = (T * (3 - beta) / 4) + delay
    t_T5 = (T * (3 + beta) / 4) + delay
    t_T6 = T

    t = np.mod(t, T)
    if t_T1 <= t < t_T2:
        f_value = amp * np.sin((2 * np.pi * (t - delay)) / (T * (1 - beta)))
    elif t_T2 <= t < t_T3:
        f_value = amp
    elif t_T3 <= t < t_T4:
        f_value = amp * np.sin(
            (2 * np.pi * (t - delay - (beta * T / 2))) / (T * (1 - beta)))
    elif t_T4 <= t < t_T5:
        f_value = -amp
    elif t_T5 <= t < t_T6:
        f_value = -amp * np.sin((2 * np.pi * (t - delay - (beta * T / 2) -
                                              (T / 2))) / (T * (1 - beta)))
    return f_value
    # ------------------------------------------
