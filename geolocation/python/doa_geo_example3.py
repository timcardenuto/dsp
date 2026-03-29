# Passive Geolocation Homework 7 Problem 3 - Tim Cardenuto
# TODO set plot axes to be based on min/max values that are being plotted. The matplotlib autosizing doesn't seem to work when I need it to be a box

import numpy as np
import matplotlib.pyplot as plt
from geolocation import *

plt.ion()
fig,axs = plt.subplots(1,1)
axs.set_xlabel('west/east (nm)')
axs.set_ylabel('south/north (nm)')
axs.set_title('DOA Geolocation Estimation Using Iterated Least Squares')

# Sensor 1
loc1 = np.vstack([[0],[0]])                                 # Sensor location x,y (meters)
theta1 = 50/180*np.pi                                       # DOA angle measurement (radians)
sigma1 = 3*np.pi/180                                        # DOA angle error, i.e. standard deviation (radians)
axs.plot(loc1[0], loc1[1], 'ok', label='Sensor 1')          # Plot sensor location
axs.plot([loc1[0], 100*np.cos(theta1+sigma1)+loc1[0]],[loc1[1], 100*np.sin(theta1+sigma1)+loc1[1]], '--k') # Plot DOA line from sensor toward target +/- sigma
axs.plot([loc1[0], 100*np.cos(theta1-sigma1)+loc1[0]],[loc1[1], 100*np.sin(theta1-sigma1)+loc1[1]], '--k')

# Sensor 2
loc2 = np.vstack([[20],[0]])
theta2 = 70/180*np.pi
sigma2 = 2*np.pi/180
axs.plot(loc2[0], loc2[1], 'ok', label='Sensor 2')
axs.plot([loc2[0], 100*np.cos(theta2+sigma2)+loc2[0]],[loc2[1], 100*np.sin(theta2+sigma2)+loc2[1]], '--k')
axs.plot([loc2[0], 100*np.cos(theta2-sigma2)+loc2[0]],[loc2[1], 100*np.sin(theta2-sigma2)+loc2[1]], '--k')

# # Sensor 3
# loc3 = np.vstack([[60],[0]])
# theta3 = 110/180*np.pi
# sigma3 = 1*np.pi/180
# plt.plot(loc3[0], loc3[1], 'ok', label='Sensor 3')
# plt.plot([loc3[0], 100*np.cos(theta3)+loc3[0]],[loc3[1], 100*np.sin(theta3)+loc3[1]], 'k')

# build arrays so we can support N sensor locations w/o hardcoding the math
# loc_array = [loc1, loc2, loc3]
# theta_array = [theta1, theta2, theta3]
# sigma_array = [sigma1, sigma2, sigma3]

loc_array = [loc1, loc2]
theta_array = [theta1, theta2]
sigma_array = [sigma1, sigma2]

ellipse = geolocate(loc_array, theta_array, sigma_array, 0.95)

# test = np.rot90(ellipse['shape'])#, k=-1)
# print(test[0].tolist())

# Plot estimated target location and EEP
axs.plot(ellipse['x'], ellipse['y'], '*g', label='Final Target Location Estimate')
axs.plot(ellipse['shape'][:,0],ellipse['shape'][:,1],'-g')
axs.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
axs.set(xlim=(-5, 100), ylim=(-5, 100))
axs.set_aspect('equal')
input("Press Enter to continue...")
