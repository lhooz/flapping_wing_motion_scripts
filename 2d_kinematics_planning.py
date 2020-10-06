"""script for 2d kinematics planning for fundamental cases"""

import numpy as np

Re = [100]
stroke_dist = [3, 5, 7, 9]
acc_dist_fraction = [0.125, 0.25, 0.5]  #-include both acc and decc dist-
pitch_dist_fraction = [0.125, 0.25, 0.5]
#---------------------------------------
wing_chord = 1
c_U = 1  #--chord length time scale always = 1--
steady_velocity = wing_chord / c_U
pitch_angle = 90
pitch_acc_fraction = 0.1
#--------------------------------------
with open('2d_case_parameters.csv', 'w') as f:
    f.write('wing_chord = %s (m)\n' % str(wing_chord))
    f.write('steady_velocity = %s (m/s)\n' % str(steady_velocity))
    f.write('pitch_angle = %s (deg)\n' % str(pitch_angle))
    f.write('pitch_acc_time_fraction = %s\n' % str(pitch_acc_fraction))
    f.write(
        '%s\n' %
        r'Re, stroke_dist, acc_dist_fraction, pitch_dist_fraction, acc_t (s), steady_t (s), pitch_t (s), nu (m^2/s), acc (m/s^2), pitch_acc (deg/s^2)'
    )
#--------------------------------------
for Rei in Re:
    for strokei in stroke_dist:
        for acc_df in acc_dist_fraction:
            acc_dfi = acc_df / 2
            for pitch_dfi in pitch_dist_fraction:
                acc_t = acc_dfi * strokei * wing_chord / (0.5 *
                                                          steady_velocity)
                steady_t = (
                    1 - 2 * acc_dfi) * strokei * wing_chord / steady_velocity
                acc = steady_velocity / acc_t
                if pitch_dfi <= acc_dfi:
                    pitch_t = (pitch_dfi * strokei * wing_chord /
                               (0.5 * acc))**0.5
                else:
                    pitch_t = acc_t
                    pitch_t += (pitch_dfi - acc_dfi
                                ) * strokei * wing_chord / steady_velocity

                pitch_avel = pitch_angle / (
                    (1 - 0.5 * pitch_acc_fraction) * pitch_t)
                pitch_acc = pitch_avel / (0.5 * pitch_acc_fraction * pitch_t)
                nu = steady_velocity * wing_chord / Rei
                parameters = [
                    Rei, strokei, acc_df, pitch_dfi, acc_t, steady_t, pitch_t,
                    nu, acc, pitch_acc
                ]
                parameters = [str(x) for x in parameters]

                with open('2d_case_parameters.csv', 'a') as f:
                    f.write('%s\n' % ', '.join(parameters))
