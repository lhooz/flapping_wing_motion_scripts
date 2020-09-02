"""script for tabulated 6DOF revolving wing motion"""

import numpy as np
from kinematic_functions import sinu_ramp_rev
from kinematics_write import kf_plotter, write_2d, write_3d

# sinumation time definition and choose ramp functions to use
ramp_time_series_length = 100
steady_rotation_time_series_length = 5000
start_time = 0
initial_ramp_time = 0.1
steady_rotation_time = 4.5
ramp_function = 'sinu_ramp'
#------------------------------
steady_rotation_frequency = 1
section_location = 1  #used only for 2d cases
#------------------------------
end_time = start_time + initial_ramp_time + steady_rotation_time
t_ramp = np.linspace(start_time, initial_ramp_time, ramp_time_series_length)
t_rev = np.linspace(start_time + initial_ramp_time, end_time,
                    steady_rotation_time_series_length)
t_rev = np.delete(t_rev, 0)

t = np.append(t_ramp, t_rev)
time_series_length = len(t)

kinematic_parameters = [steady_rotation_frequency, initial_ramp_time]

if ramp_function == 'sinu_ramp':
    kinematic_angles = sinu_ramp_rev(t, kinematic_parameters)
#--------------------------------------------------
angles_to_plot = ['phi', 'dphi']

kf_plotter(t, kinematic_angles, angles_to_plot, 'revolving')

write_2d(t, section_location, kinematic_angles, 'revolving')
write_3d(t, kinematic_angles, 'revolving')
