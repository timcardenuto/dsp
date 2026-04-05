# Passive Geolocation Homework 7 Problem 5 - Tim Cardenuto
# TDOA h(x) = r1 - r2 = ((x-a1)'*(x-a1))^.5 - ((x-a2)'*(x-a2))^.5
# H = (1/r1)*(x-a1)' - (1/r2)*(x-a2)'

# NOTE: a range *difference* (r1 - r2) to a target from two sensors is equivalent to the time difference of arrival between the two sensors
#       let r1 be the notation for the distance from a1 to x (target)
#       this is equivalent to r1 = c*(ta1 - tx), where c is speed of light (propogation) and ta1-tx is the time an RF wave emitting from the target took to get to sensor a1
#       then, r2 = c*(ta2 - tx)
#       and, r1-r2 = c*(ta1 - tx) - c*(ta2 - tx)
#                  = c*(ta1-ta2)
#       No knowlege of the time of emission (tx) is necessary for the 'range difference'. 

import numpy as np
import matplotlib.pyplot as plt
from geolocation import *

# set up plotting
plt.ion()
plt.xlabel('y (meters)')
plt.ylabel('x (meters)')
plt.title("TDOA Hyperbolas")

sigma = 50
z = np.array([13000,-5600])

# 1 nautical mile = 1852 meters
a1 = np.array([-10*1852, 20*1852])
a2 = np.array([0, 20*1852])
a3 = np.array([10*1852, 20*1852])
plt.plot(a1[0]/1852, a1[1]/1852, 'ob', label='Aircraft 1')
plt.plot(a2[0]/1852, a2[1]/1852, 'ob', label='Aircraft 2')
plt.plot(a3[0]/1852, a3[1]/1852, 'ob', label='Aircraft 3')

#########################################
# Using the Geolocation library function
# TODO:  Something wrong with the measurements vs. sensor locations.... not sure the above example makes sense? not sure which measurement goes with which and sigma is large...

# request = [
#     {'mtype': 'tdoa_range', 'value':  z[0], 'sigma': sigma, 's1loc': a1, 's2loc': a2},
#     {'mtype': 'tdoa_range', 'value':  z[1], 'sigma': sigma, 's1loc': a1, 's2loc': a3}
#     ]
# ellipse = geolocate(request, 0.95)

# plt.plot(ellipse['x'], ellipse['y'], '*g', label='Final Target Location Estimate')
# plt.plot(ellipse['shape'][:,0],ellipse['shape'][:,1],'-g', label='Error Ellipse')
# plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
# plt.tight_layout()
# plt.draw()

# input("Press Enter to continue...")

#########################################
# Original code here

print("sigma: "+str(sigma))
print("z: "+str(z))
print("a1:"+str(a1))
print("a2:"+str(a2))
print("a3:"+str(a3))
print("\n")

R = (sigma*sigma)*np.identity(2)
xhat = np.array([0,0])
plt.plot(xhat[0]/1852, xhat[1]/1852, '*g', label='Estimated location')

count = 1
maxcount = 10
error = [100,100]
maxerror = .1
while (count < maxcount and (error[0]+error[1]/2) > maxerror):

    h = np.array([
        np.power(np.transpose(xhat-a1).dot(xhat-a1),0.5) - np.power(np.transpose(xhat-a2).dot(xhat-a2),0.5), 
        np.power(np.transpose(xhat-a3).dot(xhat-a3),0.5) - np.power(np.transpose(xhat-a2).dot(xhat-a2),0.5)
        ])

    H = np.array([
        (1/np.power(np.transpose(xhat-a1).dot(xhat-a1),0.5))*np.transpose(xhat-a1) - (1/np.power(np.transpose(xhat-a2).dot(xhat-a2),0.5))*np.transpose(xhat-a2), 
        (1/np.power(np.transpose(xhat-a3).dot(xhat-a3),0.5))*np.transpose(xhat-a3) - (1/np.power(np.transpose(xhat-a2).dot(xhat-a2),0.5))*np.transpose(xhat-a2)
        ])

    P = np.linalg.inv(np.transpose(H).dot(np.linalg.inv(R)).dot(H))

    xhatnew = xhat + P.dot(np.transpose(H)).dot(np.linalg.inv(R)).dot(z-h)
    error = abs(xhatnew - xhat)
    xhat = xhatnew
   
    print("count = ",count) 
    print("Covariance = ",P)
    print("xhat = ",xhat)
    print("error = ",error)
    print("")
    count = count+1

# Create and scale ellipse
w,v = np.linalg.eig(P)
print("Eigenvalues = ",w)
k = 5.9915; # 95% EEP
semimajor = np.sqrt(k*max(w))
semiminor = np.sqrt(k*min(w))
theta = np.linspace(0,2*np.pi)
ellipse = np.transpose(np.vstack([semimajor*np.cos(theta), semiminor*np.sin(theta)]))

# Rotate ellipse
eigenvector = v[:,1]
print("Eigenvector = ",eigenvector)
angle = np.arctan2(eigenvector[1], eigenvector[0])
if(angle < 0):
    angle = angle + 2*np.pi
R = np.vstack([[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]])
ellipse = ellipse.dot(R)

# Shift ellipse and plot 
plt.plot((ellipse[:,0]+xhat[0])/1852,(ellipse[:,1]+xhat[1])/1852,'-', label='Error Ellipse')
plt.xlabel('south/north (nm)')
plt.ylabel('west/east (nm)')
plt.title('Geolocation Estimate')
plt.show()

input("Press Enter to continue...")
