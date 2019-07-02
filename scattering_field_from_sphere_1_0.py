# ------------------------------------------------------
#################################################
# Code Dependencies: Library: Holopy 3.0
# The unit for distance in this code is microns #
#################################################

## If the figures do not pop out in a different window,
## use this command in the console window:
## %matplotlib qt



import holopy as hp
from holopy.scattering import calc_holo, Sphere
import numpy as np
from holopy.scattering import calc_scat_matrix, calc_intensity, calc_field
import matplotlib.pyplot as plt

medium_index = 1.33
illum_wavelen = 0.750 # 750nm
illum_polarization = (1,0) # x polarized
sphere_z_position = 2 # in microns
detector = hp.detector_grid(shape = 201, spacing = .1)
distant_sphere = Sphere(r=0.5, n=1.59, center = (10.05, 10.05, sphere_z_position))

## Calculate the field here given above paramters
scat = calc_field(detector, distant_sphere, medium_index, illum_wavelen, illum_polarization);
scattered_matrix = scat.data
scat_intensity_x = scattered_matrix[0,:,:,0]
scat_intensity_y = scattered_matrix[1,:,:,0]
scat_intensity_z = scattered_matrix[2,:,:,0]

# real value
scat_x_real = scat_intensity_x.real
# imaginary value
scat_x_imag = scat_intensity_x.imag
# absolute value
scat_x = np.multiply(scat_x_real,scat_x_real) + np.multiply(scat_x_imag, scat_x_imag)
np.sqrt(scat_x)
plt.figure("Scattered Intensity X")
plt.imshow(scat_x)

# real value
scat_y_real = scat_intensity_y.real
# imaginary value
scat_y_imag = scat_intensity_y.imag
# absolute value
scat_y = np.multiply(scat_y_real,scat_y_real) + np.multiply(scat_y_imag, scat_y_imag)
np.sqrt(scat_y)
plt.figure("Scattered Intensity Y")
plt.imshow(scat_y)

# real value
scat_z_real = scat_intensity_z.real
# imaginary value
scat_z_imag = scat_intensity_z.imag
# absolute value
scat_z = np.multiply(scat_z_real,scat_z_real) + np.multiply(scat_z_imag, scat_z_imag)
np.sqrt(scat_z)
plt.figure("Scattered Intensity Z")
plt.imshow(scat_z)