# Passive Geolocation Homework 8 Problem 2 - Tim Cardenuto

import numpy as np
import matplotlib.pyplot as plt

plt.ion()
plt.xlabel('west/east (nm)')
plt.ylabel('south/north (nm)')
plt.title('Geolocation Estimation Using Iterated Least Squares')



def calculateEllipse(xhat, P, alpha):
    """ Computes an Elliptical Error Probable (EEP), i.e. an error ellipse
        The center of the ellipse is at the measured or estimated location (xhat)
        The shape (how elliptical or circular), orientation (angular rotation), and size (distance of semimajor/semiminor axis) is based on the estimated covariance matrix (P) and the desired containment (alpha). 
        The covariance matrix (P) incorporates the errors of the original measurements used to estimate the center point. If there was zero error in the original measurements then the ellipse would be zero size (a point). 
        The desired containment (alpha) is interpreted as meaning the ellipse should be scaled big enough to provide an alpha probability of the target residing within the ellipse. The higher you want this to be, the bigger the ellipse will be. Common examples of alpha are 0.50 or 0.95. 
    """
    # Calculate containment critical value k from alpha
    k = -2 * np.log(1-alpha)    # for 95% EEP, k = 5.9915 
    
    # Calculate ellipse shape/size
    eigenvalues,eigenvectors = np.linalg.eig(P)
    semimajor = np.sqrt(k*max(eigenvalues))
    semiminor = np.sqrt(k*min(eigenvalues))

    # Create ellipse 
    theta = np.linspace(0,2*np.pi)
    ellipse = np.transpose(np.vstack([semimajor*np.cos(theta), semiminor*np.sin(theta)]))
    
    # Rotate ellipse
    eigenvector = eigenvectors[:,1]
    angle = np.arctan2(eigenvector[1], eigenvector[0])
    if(angle < 0):
        angle = angle + 2*np.pi
    R = np.vstack([
        [np.cos(angle), np.sin(angle)],
        [-np.sin(angle), np.cos(angle)]
        ])
    ellipse = ellipse @ R
    
    # Shift ellipse to center point
    ellipse = np.hstack([
        np.vstack(ellipse[:,0]+xhat[0]), 
        np.vstack(ellipse[:,1]+xhat[1])
        ])
        
    print("Eigen Values = ",eigenvalues)
    print("Eigen Vectors = ",eigenvectors)
    print("Semi-major Axis = ",semimajor)
    print("Semi-minor Axis = ",semiminor)
    print("")
    
    return ellipse
    
    

# Plot True Target location
plt.plot(0,0,'or', label='True Target Location')

# Sensor 1 - DOA measurement
a1 = np.vstack([[(np.sin(60/180*np.pi)*15*1852)], [(np.cos(60/180*np.pi)*15*1852)]]) # Sensor location x,y (meters)
theta1 = -(90+60)/180*np.pi                                                         # DOA angle measurement (radians)
sigma1 = 2/180*np.pi                                                                # DOA angle error, i.e. standard deviation (radians)
u1 = np.vstack([[np.cos(theta1)], [np.sin(theta1)]])                                 # unit vector along DOA
M1 = np.identity(2) - u1@(u1.conj().transpose())                                    # Moore-Penrose matrix
plt.plot(a1[0]/1852, a1[1]/1852, 'ok', label="Sensor 1")
plt.plot([a1[0][0]/1852, (100*np.cos(theta1)/1852)],[a1[1][0]/1852, (100*np.sin(theta1)/1852)], 'k')

# Sensor 2 - range measurement
a2 = np.vstack([[-(np.sin(60/180*np.pi)*10*1852)], [(np.cos(60/180*np.pi)*10*1852)]])
range2 = 10*1852
theta2 = -(90-60)/180*np.pi
sigma2 = 100
u2 = np.vstack([[np.cos(theta2)], [np.sin(theta2)]])
M2 = np.identity(2) - u2@(u2.conj().transpose())
plt.plot(a2[0]/1852, a2[1]/1852, 'ob', label="Sensor 2")
th = np.arange(0, 2*np.pi, np.pi/50)
xunit = range2/1852 * np.cos(th) + a2[0]/1852
yunit = range2/1852 * np.sin(th) + a2[1]/1852
plt.plot(xunit, yunit, 'b')

# Initial estimate of target position (initial guess) given by Moore-Penrose algorithm
xhat = np.linalg.inv(M1+M2)@(M1@a1+M2@a2)
plt.plot(xhat[0]/1852, xhat[1]/1852, '*r', label='Initial Target Location Estimate w/ Moore-Penrose')


# Final estimation by Iterated Least Squares 
z = np.vstack([[theta1], [range2]])                     # measurements
R = np.vstack([[sigma1*sigma1, 0], [0, sigma2*sigma2]]) # covariance matrix (measurement errors)

count = 1
maxcount = 10               # limit number of runs w/ a maximum count
error = np.hstack([100,100])
maxerror = 0.1              # when to stop iterating (acceptable error)
while (count < maxcount and (error[0]+error[1]/2) > maxerror):
    print("Iteration: ",count)
    # DOA h and H
    # h(x) = theta = arctan((x2-a2),(x1-a1))
    # H = (1/r)*[-sin(h) cos(h)]
    # Range h and H
    # h = r = ((x-a)'*(x-a))^.5
    # H = (1/r)*(x-a)'
    
    h = np.vstack([
        np.arctan2((xhat[1]-a1[1]),(xhat[0]-a1[0])),
        np.power((xhat-a2).conj().transpose()@(xhat-a2),0.5)
        ])

    H = np.vstack([
        (1/np.power((xhat-a1).conj().transpose()@(xhat-a1),0.5))*np.hstack([-np.sin(h[0]), np.cos(h[0])]), 
        (1/np.power((xhat-a2).conj().transpose()@(xhat-a2),0.5))@(xhat-a2).conj().transpose()
        ])
        
    P = np.linalg.inv(H.conj().transpose()@np.linalg.inv(R)@H)
    xhatnew = xhat + P@H.conj().transpose()@np.linalg.inv(R)@(z-h)
    error = abs(xhatnew - xhat)
    xhat = xhatnew
    count = count+1
    
    print("Error: ",error)
    plt.plot(xhat[0], xhat[1], '*k', label="Estimate Update, Error={}".format(error))
    plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
    plt.tight_layout()
    plt.draw()
    print("")



# Calculate the Elliptical Error Probable (EEP), i.e. error ellipse
ellipse = calculateEllipse(xhat, P, 0.95) 

# Plot estimated target location and EEP
plt.plot(xhat[0]/1852, xhat[1]/1852, '*g', label='Final Location Estimate')
plt.plot(ellipse[:,0]/1852, ellipse[:,1]/1852, '-g', label='95% Ellipse')
plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
input("Press Enter to continue...")


