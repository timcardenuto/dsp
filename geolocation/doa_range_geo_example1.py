# Passive Geolocation Homework 8 Problem 2 - Tim Cardenuto

import numpy as np
import matplotlib.pyplot as plt


def calculateEllipse(xhat,P):
    # Create and scale ellipse
    eigenvalues,eigenvectors = np.linalg.eig(P)
    k = 5.9915 # 95% EEP
    semimajor = np.sqrt(k*max(eigenvalues))
    semiminor = np.sqrt(k*min(eigenvalues))
    theta = np.linspace(0,2*np.pi)
    ellipse = np.transpose(np.vstack([semimajor*np.cos(theta), semiminor*np.sin(theta)]))
    print("Eigen Values = ",eigenvalues)
    print("Eigen Vectors = ",eigenvectors)
    print("Semi-major Axis = ",semimajor)
    print("Semi-minor Axis = ",semiminor)
    print("")
    # Rotate ellipse
    eigenvector = eigenvectors[:,1]
    angle = np.arctan2(eigenvector[1], eigenvector[0])
    if(angle < 0):
        angle = angle + 2*np.pi
    R = np.vstack([[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]])
    ellipse = ellipse @ R
    
    # Shift ellipse 
    #ellipse[:,0] = ellipse[:,0]+xhat[0]
    #ellipse[:,1] = ellipse[:,1]+xhat[1]
    return ellipse
    
    


# Plot True Target location
plt.plot(0,0,'or', label='True Target Location')


# Buoy 1 - DOA
a1 = np.array([[(np.sin(60/180*np.pi)*15*1852)], [(np.cos(60/180*np.pi)*15*1852)]])
range1 = 15*1852
theta1 = -(90+60)/180*np.pi
sigma1 = 2/180*np.pi
u1 = np.array([[np.cos(theta1)], [np.sin(theta1)]])
M1 = np.identity(2) - u1@(u1.conj().transpose())
plt.plot(a1[0]/1852, a1[1]/1852, 'ok', label="Buoys")
plt.plot([a1[0][0]/1852, (100*np.cos(theta1)/1852)],[a1[1][0]/1852, (100*np.sin(theta1)/1852)], 'k')

# Buoy 2 - range
a2 = np.array([[-(np.sin(60/180*np.pi)*10*1852)], [(np.cos(60/180*np.pi)*10*1852)]])
range2 = 10*1852
theta2 = -(90-60)/180*np.pi
sigma2 = 100
u2 = np.array([[np.cos(theta2)], [np.sin(theta2)]])
M2 = np.identity(2) - u2@(u2.conj().transpose())
plt.plot(a2[0]/1852, a2[1]/1852, 'ok')
th = np.arange(0, 2*np.pi, np.pi/50)
xunit = range2/1852 * np.cos(th) + a2[0]/1852
yunit = range2/1852 * np.sin(th) + a2[1]/1852
plt.plot(xunit, yunit, 'k')

# Initial position estimate given by Moore-Penrose solution
xhat = np.linalg.inv(M1+M2)@(M1@a1+M2@a2)
plt.plot(xhat[0]/1852, xhat[1]/1852, '*k', label='Initial Estimate')

# Final Estimation by Iterated Least Squares 
# DOA h and H
# h(x) = theta = arctan((x2-a2),(x1-a1))
# H = (1/r)*[-sin(h) cos(h)]
# Range h and H
# h = r = ((x-a)'*(x-a))^.5
# H = (1/r)*(x-a)'
count = 1
maxcount = 10
error = [100,100]
maxerror = 0.1
xhatold = xhat
# Using two measurements - 1 DOA and 1 range
z = np.array([[theta1], [range2]])
R = np.array([[sigma1*sigma1, 0], [0, sigma2*sigma2]])
while (count < maxcount and (error[0]+error[1]/2) > maxerror):
    print("Iteration: ",count)
    #h = [atan2((xhatold(2)-a1(2)),(xhatold(1)-a1(1))); ...
    #    ((xhatold-a2)'*(xhatold-a2))^.5]
    h = np.array([[np.arctan2((xhatold[1]-a1[1]),(xhatold[0]-a1[0])).item()],  [np.power((xhatold-a2).conj().transpose()@(xhatold-a2),0.5).item()]])
    #H = [(1/(((xhatold-a1)'*(xhatold-a1))^.5))*[-np.sin(h(1)) np.cos(h(1))]; ...
    #    (1/(((xhatold-a2)'*(xhatold-a2))^.5))*(xhatold-a2)']
    H = np.array([
        [(1/np.power((xhatold-a1).conj().transpose()@(xhatold-a1),0.5))*np.array([-np.sin(h[0].item()), np.cos(h[0].item())])][0][0], 
        [(1/np.power((xhatold-a2).conj().transpose()@(xhatold-a2),0.5))@(xhatold-a2).conj().transpose()][0][0]
        ])
    P = np.linalg.inv(H.conj().transpose()@np.linalg.inv(R)@H)
    xhatnew = xhatold + P@H.conj().transpose()@np.linalg.inv(R)@(z-h)
    error = abs(xhatnew - xhatold)
    xhatold = xhatnew
    count = count+1
    print("Error: ",error)
    print("")




# Plot estimated target location and EEP
plt.plot(xhatnew[0]/1852, xhatnew[1]/1852, '*b', label='Final Location Estimate')
ellipse = calculateEllipse(xhatnew, P) 
plt.plot((ellipse[:,0]+xhatnew[0])/1852,(ellipse[:,1]+xhatnew[1])/1852,'-b', label='95% Ellipse')
plt.xlabel('west/east (nm)')
plt.ylabel('south/north (nm)')
plt.title('Geolocation Estimate')
plt.legend(loc='upper right', shadow=True)
plt.show()

