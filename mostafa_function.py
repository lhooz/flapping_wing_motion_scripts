"""sinusiodal motion kinematic_function definition"""
import numpy as np


def kinematic_function(f, amp, tr_fac, del_fac, t):
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
