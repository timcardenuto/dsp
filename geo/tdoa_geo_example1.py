# Passive Geolocation Homework 7 Problem 5 - Tim Cardenuto
# TDOA h(x) = r1 - r2 = ((x-a1)'*(x-a1))^.5 - ((x-a2)'*(x-a2))^.5
# H = (1/r1)*(x-a1)' - (1/r2)*(x-a2)'

import numpy as np
import matplotlib.pyplot as plt

sigma = 50
R = (sigma*sigma)*np.identity(2)
xhat = np.array([0,0])
z = np.array([13000,-5600])
# 1 nautical mile = 1852 meters
a1 = np.array([-10*1852, 20*1852])
a2 = np.array([0, 20*1852])
a3 = np.array([10*1852, 20*1852])
count = 1
maxcount = 10
error = [100,100]
maxerror = .1


while (count < maxcount and (error[0]+error[1]/2) > maxerror):

    h = np.array([np.power(np.transpose(xhat-a1).dot(xhat-a1),0.5) - np.power(np.transpose(xhat-a2).dot(xhat-a2),0.5), np.power(np.transpose(xhat-a3).dot(xhat-a3),0.5) - np.power(np.transpose(xhat-a2).dot(xhat-a2),0.5)])

    H = np.array([(1/np.power(np.transpose(xhat-a1).dot(xhat-a1),0.5))*np.transpose(xhat-a1) - (1/np.power(np.transpose(xhat-a2).dot(xhat-a2),0.5))*np.transpose(xhat-a2), (1/np.power(np.transpose(xhat-a3).dot(xhat-a3),0.5))*np.transpose(xhat-a3) - (1/np.power(np.transpose(xhat-a2).dot(xhat-a2),0.5))*np.transpose(xhat-a2)])

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

    
# 95% EEP
w,v = np.linalg.eig(P)
print("Eigenvalues = ",w)
k = 5.9915;
semimajor = np.sqrt(k*max(w))
semiminor = np.sqrt(k*min(w))


plt.plot(a1[0]/1852, a1[1]/1852, 'ob', label='Aircraft 1')
plt.plot(a2[0]/1852, a2[1]/1852, 'ob', label='Aircraft 2')
plt.plot(a3[0]/1852, a3[1]/1852, 'ob', label='Aircraft 3')
plt.plot(xhat[0]/1852, xhat[1]/1852, '*g', label='Estimated location')

# Create and scale ellipse
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
