"""kinematic functions defknitions: default (smooth) and sinusoidal (sinusiodal)"""

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
    ptf_coefficient = kinematic_parameters[7]

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

    def ddphi(x):
        """flapping angular acceleration function"""
        return derivative(dphi, x, dx=1e-6)

    def alf(x):
        """pitching motion function"""
        if pitching_time_coefficient == 0:
            alft = pitching_amplitude * np.sin(
                2 * np.pi * flapping_wing_frequency *
                (x - pitching_delay_time_fraction / flapping_wing_frequency))
        elif pitching_time_coefficient == 'f':

            ptc_x = ptf_function(x)

            alft = pitching_amplitude / np.tanh(ptc_x) * np.tanh(
                ptc_x * np.sin(2 * np.pi * flapping_wing_frequency *
                               (x - pitching_delay_time_fraction /
                                flapping_wing_frequency)))
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

    def ddalf(x):
        """pitching angular acceleration function"""
        return derivative(dalf, x, dx=1e-6)

    def ptf_function(x):
        """time varying pitching time coefficient function"""
        return ptf_coefficient + ptf_coefficient * (np.sin(
            2 * np.pi * flapping_wing_frequency *
            (x - pitching_delay_time_fraction / flapping_wing_frequency)))**2

    kinematic_angles = []
    t_1st_cycle = [t1 for t1 in t if t1 <= 1 / flapping_wing_frequency]
    for ti in t_1st_cycle:
        kinematic_anglesi = [
            phi(ti),
            alf(ti),
            dphi(ti),
            dalf(ti),
            ddphi(ti),
            ddalf(ti)
        ]
        kinematic_angles.append(kinematic_anglesi)

    return kinematic_angles


# -------------------------------------------------


# sinu_continuous kinematic functions
# -------------------------------------------------
def sinu_continuous_kinematic_function(t, kinematic_parameters):
    """definition for sinusoidal (sinusiodal) kinematic functiontion"""
    flapping_wing_frequency = kinematic_parameters[0]
    flapping_angular_velocity_amplitude = kinematic_parameters[1]
    pitching_angular_velocity_amplitude = kinematic_parameters[2]
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
    print('flapping amplitude = %s' % (flapping_amplitude / 2))

    def ddphi(x):
        """flapping angular acceleration function"""
        return derivative(dphi, x, dx=1e-6)

    initial_phi = integrate.quad(
        lambda x: dphi(x), 0,
        np.abs(flapping_delay_time_fraction) / flapping_wing_frequency)[0]
    initial_phi = -np.sign(flapping_delay_time_fraction) * initial_phi

    dphi_int = []
    for ti in t:
        dphi_int.append(dphi(ti))

    def phi(x):
        """flapping motion function"""
        t_int = [tx for tx in t if tx <= x]
        time_array_length = len(t_int)
        # print(time_array_length, dphi_int[0:time_array_length])
        return flapping_amplitude / 4 + initial_phi + integrate.simps(
            dphi_int[0:time_array_length], t_int)

    def dalf(x):
        """flapping angular velocity function"""
        return kf_continuous(flapping_wing_frequency,
                             pitching_angular_velocity_amplitude,
                             pitching_time_fraction,
                             pitching_delay_time_fraction, x)

    pitching_amplitude = integrate.quad(lambda x: np.abs(dalf(x)), 0,
                                        1 / flapping_wing_frequency)[0]
    print('pitching amplitude = %s' % (pitching_amplitude / 2))

    def ddalf(x):
        """pitching angular acceleration function"""
        return derivative(dalf, x, dx=1e-6)

    initial_alf = integrate.quad(
        lambda x: dalf(x), 0,
        np.abs(pitching_delay_time_fraction) / flapping_wing_frequency)[0]
    initial_alf = -np.sign(pitching_delay_time_fraction) * initial_alf

    dalf_int = []
    for ti in t:
        dalf_int.append(dalf(ti))

    def alf(x):
        """pitching motion function"""
        t_int = [tx for tx in t if tx <= x]
        time_array_length = len(t_int)
        # print(time_array_length, dphi_int[0:time_array_length])
        return initial_alf + integrate.simps(dalf_int[0:time_array_length],
                                             t_int)

    kinematic_angles = []
    t_1st_cycle = [t1 for t1 in t if t1 <= 1 / flapping_wing_frequency]
    for ti in t_1st_cycle:
        kinematic_anglesi = [
            phi(ti),
            alf(ti),
            dphi(ti),
            dalf(ti),
            ddphi(ti),
            ddalf(ti)
        ]
        kinematic_angles.append(kinematic_anglesi)

    return kinematic_angles


# sinusoidal kinematic functions
# -------------------------------------------------
def sinusoidal_kinematic_function(t, kinematic_parameters):
    """definition for sinusoidal (sinusiodal) kinematic functiontion"""
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
    print('flapping amplitude = %s' % (flapping_amplitude / 2))

    def ddphi(x):
        """flapping angular acceleration function"""
        return derivative(dphi, x, dx=1e-6)

    initial_phi = integrate.quad(
        lambda x: dphi(x), 0,
        np.abs(flapping_delay_time_fraction) / flapping_wing_frequency)[0]
    initial_phi = -np.sign(flapping_delay_time_fraction) * initial_phi

    dphi_int = []
    for ti in t:
        dphi_int.append(dphi(ti))

    def phi(x):
        """flapping motion function"""
        t_int = [tx for tx in t if tx <= x]
        time_array_length = len(t_int)
        # print(time_array_length, dphi_int[0:time_array_length])
        return flapping_amplitude / 4 + initial_phi + integrate.simps(
            dphi_int[0:time_array_length], t_int)

    def alf(x):
        """pitching motion function"""
        return kf(flapping_wing_frequency, pitching_amplitude,
                  pitching_time_fraction, pitching_delay_time_fraction, x)

    def dalf(x):
        """pitching angular velocity function"""
        return derivative(alf, x, dx=1e-6)

    def ddalf(x):
        """pitching angular acceleration function"""
        return derivative(dalf, x, dx=1e-6)

    kinematic_angles = []
    t_1st_cycle = [t1 for t1 in t if t1 <= 1 / flapping_wing_frequency]
    for ti in t_1st_cycle:
        kinematic_anglesi = [
            phi(ti),
            alf(ti),
            dphi(ti),
            dalf(ti),
            ddphi(ti),
            ddalf(ti)
        ]
        kinematic_angles.append(kinematic_anglesi)

    return kinematic_angles


    # original sinusoidal (sinusiodal) function
def kf(f, amp, tr_fac, del_fac, t):
    """kinematic_function for flapping velocity and angle of attack"""
    T = 1 / f
    delay = del_fac * T
    tr = tr_fac * T
    beta = 1 - (2 * tr / T)

    t_T1 = 0
    t_T2 = (T * (1 - beta) / 4)
    t_T3 = (T * (1 + beta) / 4)
    t_T4 = (T * (3 - beta) / 4)
    t_T5 = (T * (3 + beta) / 4)
    t_T6 = T

    t = np.mod(t - delay, T)
    if t_T1 <= t < t_T2:
        f_value = amp * np.sin((2 * np.pi * t) / (T * (1 - beta)))
    elif t_T2 <= t < t_T3:
        f_value = amp
    elif t_T3 <= t < t_T4:
        f_value = amp * np.sin(
            (2 * np.pi * (t - (beta * T / 2))) / (T * (1 - beta)))
    elif t_T4 <= t < t_T5:
        f_value = -amp
    elif t_T5 <= t < t_T6:
        f_value = amp * np.sin((2 * np.pi * (t - beta * T)) / (T * (1 - beta)))
    return f_value
    # ------------------------------------------


    # original sinusoidal (sinusiodal) function
def kf_continuous(f, amp, tr_fac, del_fac, t):
    """kinematic_function for flapping velocity and angle of attack"""
    T = 1 / f
    T_sinu = 0.5 * T
    delay = del_fac * T
    tr = tr_fac * T
    beta = 1 - (2 * tr / T)

    t_T1 = 0
    t_T2 = (T * (1 - beta) / 4)
    t_T3 = (T * (1 + beta) / 4)
    t_T4 = (T * (3 - beta) / 4)
    t_T5 = (T * (3 + beta) / 4)
    t_T6 = T

    t = np.mod(t - delay, T)
    if t_T1 <= t < t_T2:
        f_value = amp / 2 + amp / 2 * np.cos(
            (2 * np.pi * t) / (T_sinu * (1 - beta)))
    elif t_T2 <= t < t_T3:
        f_value = 0
    elif t_T3 <= t < t_T4:
        f_value = -amp / 2 - amp / 2 * np.cos(
            (2 * np.pi * (t - (beta * T / 2))) / (T_sinu * (1 - beta)))
    elif t_T4 <= t < t_T5:
        f_value = 0
    elif t_T5 <= t < t_T6:
        f_value = amp / 2 + amp / 2 * np.cos(
            (2 * np.pi * (t - beta * T)) / (T_sinu * (1 - beta)))
    return f_value
    # ------------------------------------------


    # revolving wing kinematics with sinusiodal ramp function
def sinu_ramp_rev(t, kinematic_parameters):
    """revolving wing kinematics function with sinusoidal ramp at start"""
    steady_rotation_frequency = kinematic_parameters[0]
    initial_ramp_time = kinematic_parameters[1]

    steady_rotation_omega = 360 * steady_rotation_frequency
    omega_print = steady_rotation_omega * np.pi / 180
    print('steady revolving omega = %s' % omega_print)

    def omega(x):
        """rotation speed function"""
        if x <= initial_ramp_time:
            omega = steady_rotation_omega / 2 + steady_rotation_omega / 2 * np.sin(
                2 * np.pi * x / (2 * initial_ramp_time) - np.pi / 2)
        else:
            omega = steady_rotation_omega

        return omega

    def ddphi(x):
        """flapping angular acceleration function"""
        return derivative(omega, x, dx=1e-6)

    ramp_angle = integrate.quad(lambda x: np.abs(omega(x)), 0,
                                initial_ramp_time)[0]
    print('initial sinu ramp angle = %s' % ramp_angle)

    omega_int = []
    for ti in t:
        omega_int.append(omega(ti))

    def phi(x):
        """rotation angle function"""
        if x <= initial_ramp_time:
            t_int = [tx for tx in t if tx <= x]
            time_array_length = len(t_int)
            # print(time_array_length, dphi_int[0:time_array_length])
            return integrate.simps(omega_int[0:time_array_length], t_int)
        else:
            return ramp_angle + steady_rotation_omega * (x - initial_ramp_time)

    kinematic_angles = []
    for ti in t:
        kinematic_anglesi = [-phi(ti), 0, -omega(ti), 0, -ddphi(ti), 0]
        kinematic_angles.append(kinematic_anglesi)

    return kinematic_angles


#---linear ramp function smoothed at conner for revolving or 2d translating wing--
def smooth_linear_ramp(t, kinematic_parameters):
    """smoothed linear ramp function"""
    ramp_stage_acceleration = kinematic_parameters[0]
    ramp_start_time = kinematic_parameters[1]
    i_ramp_end_time = kinematic_parameters[2]
    steady_end_time = kinematic_parameters[3]
    end_ramp_end_time = kinematic_parameters[4]
    smooth_factor = kinematic_parameters[5]
    ramp_mode = kinematic_parameters[6]
    ramp_constant_time = kinematic_parameters[7]
    pitch_mode = kinematic_parameters[8]
    pitch_time = kinematic_parameters[9]
    pitch_delay_time_fraction = kinematic_parameters[10]
    pitch_acceleration = kinematic_parameters[11]
    pitch_acc_time_fraction = kinematic_parameters[12]
    section_location = kinematic_parameters[13]

    def omega(x):
        """linear ramp rotation speed function"""
        if x <= end_ramp_end_time + 2 * ramp_constant_time:
            f_t0 = np.cosh(smooth_factor * (x - ramp_start_time))
            f_t1 = np.cosh(smooth_factor * (x - i_ramp_end_time))
            if ramp_mode == 'with_end_acc':
                f_t2 = np.cosh(smooth_factor * (x - steady_end_time))
                f_t3 = np.cosh(smooth_factor * (x - end_ramp_end_time))
            elif ramp_mode == 'no_end_acc':
                f_t2 = np.cosh(smooth_factor * ramp_start_time)
                f_t3 = np.cosh(smooth_factor * i_ramp_end_time)

            omegax = (ramp_stage_acceleration / 2) / smooth_factor * np.log(
                f_t0 / f_t1 * f_t3 / f_t2)
        else:
            omegax = 0

        return omegax

    steady_rotation_omega = omega((i_ramp_end_time + steady_end_time) / 2)
    omega_print = steady_rotation_omega * np.pi / 180
    print('steady revolving omega = %s' % omega_print)

    def ddphi(x):
        """flapping angular acceleration function"""
        return derivative(omega, x, dx=1e-6)

    ramp_angle = integrate.quad(lambda x: np.abs(omega(x)),
                                ramp_start_time - ramp_constant_time,
                                i_ramp_end_time + ramp_constant_time)[0]
    print('initial linear ramp angle = %s' % ramp_angle)

    if ramp_mode == 'with_end_acc':
        end_ramp_angle = integrate.quad(
            lambda x: np.abs(omega(x)), steady_end_time - ramp_constant_time,
            end_ramp_end_time + ramp_constant_time)[0]
        print('end linear ramp angle = %s' % end_ramp_angle)

    omega_int = []
    for ti in t:
        omega_int.append(omega(ti))

    def phi(x):
        """rotation angle function"""
        t_int = [tx for tx in t if tx <= x]
        time_array_length = len(t_int)
        return integrate.simps(omega_int[0:time_array_length], t_int)

    st_dist = np.abs(phi(t[-1])) * np.pi / 180 * section_location
    print('2d wing travel distance = %s\n' % st_dist)

    #--pitching motion functions--
    if pitch_mode == 'with_end_pitch':
        pitch_delay_time = pitch_time * pitch_delay_time_fraction
        pitch_acc_time = pitch_time * pitch_acc_time_fraction / 2

        pitch_start_time = end_ramp_end_time - pitch_time + pitch_delay_time
        p_acc_end_time = pitch_start_time + pitch_acc_time
        pitch_end_time = pitch_start_time + pitch_time
        p_decc_start_time = pitch_end_time - pitch_acc_time

        def dalf(x):
            """linear ramp pitch speed function"""
            if x <= pitch_end_time + 2 * ramp_constant_time and x >= pitch_start_time - 2 * ramp_constant_time:
                f_t0 = np.cosh(smooth_factor * (x - pitch_start_time))
                f_t1 = np.cosh(smooth_factor * (x - p_acc_end_time))
                f_t2 = np.cosh(smooth_factor * (x - p_decc_start_time))
                f_t3 = np.cosh(smooth_factor * (x - pitch_end_time))

                dalfx = (pitch_acceleration / 2) / smooth_factor * np.log(
                    f_t0 / f_t1 * f_t3 / f_t2)
            else:
                dalfx = 0

            return dalfx

        steady_pitching_omega = dalf((pitch_start_time + pitch_end_time) / 2)
        omega_print = steady_pitching_omega * np.pi / 180
        print('steady wing pitch omega = %s' % omega_print)

        def ddalf(x):
            """flapping angular acceleration function"""
            return derivative(dalf, x, dx=1e-6)

        dalf_int = []
        for ti in t:
            dalf_int.append(dalf(ti))

        def alf(x):
            """rotation angle function"""
            t_int = [tx for tx in t if tx <= x]
            time_array_length = len(t_int)
            return integrate.simps(dalf_int[0:time_array_length], t_int)

        print('wing pitch angle = %s' % alf(t[-1]))

    kinematic_angles = []
    for ti in t:
        if pitch_mode == 'no_end_pitch':
            kinematic_anglesi = [-phi(ti), 0, -omega(ti), 0, -ddphi(ti), 0]
        elif pitch_mode == 'with_end_pitch':
            kinematic_anglesi = [
                -phi(ti), -alf(ti), -omega(ti), -dalf(ti), -ddphi(ti),
                -ddalf(ti)
            ]
        kinematic_angles.append(kinematic_anglesi)

    return kinematic_angles
