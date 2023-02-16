# Passive Geolocation Homework 7 Problem 3 - Tim Cardenuto

import time
import numpy as np
import matplotlib.pyplot as plt

plt.ion()
plt.xlabel('west/east (nm)')
plt.ylabel('south/north (nm)')
plt.title('DOA Geolocation Estimation Using Iterated Least Squares')


    
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
    
    
    
a1 = np.array([[0],[0]])
theta1 = 50/180*np.pi
sigma1 = 3*np.pi/180
u1 = np.array([[np.cos(theta1)], [np.sin(theta1)]])
M1 = np.identity(2) - u1@(u1.conj().transpose())
plt.plot(a1[0], a1[1], 'ok', label='Sensor 1')
plt.plot([a1[0], 100*np.cos(theta1)+a1[0]],[a1[1], 100*np.sin(theta1)+a1[1]], 'k')

a2 = np.array([[20],[0]])
theta2 = 70/180*np.pi
sigma2 = 2*np.pi/180
u2 = np.array([[np.cos(theta2)], [np.sin(theta2)]])
M2 = np.identity(2) - u2@(u2.conj().transpose())
plt.plot(a2[0], a2[1], 'ok', label='Sensor 2')
plt.plot([a2[0], 100*np.cos(theta2)+a2[0]],[a2[1], 100*np.sin(theta2)+a2[1]], 'k')

a3 = np.array([[60],[0]])
theta3 = 110/180*np.pi
sigma3 = 1*np.pi/180
u3 = np.array([[np.cos(theta3)], [np.sin(theta3)]])
M3 = np.identity(2) - u3@(u3.conj().transpose())
plt.plot(a3[0], a3[1], 'ok', label='Sensor 3')
plt.plot([a3[0], 100*np.cos(theta3)+a3[0]],[a3[1], 100*np.sin(theta3)+a3[1]], 'k')

# Initial position estimate given by Moore-Penrose solution
xhat = np.linalg.inv(M1+M2+M3)@(M1@a1+M2@a2+M3@a3)
plt.plot(xhat[0], xhat[1], '*r', label='Initial Target Location Estimate w/ Moore-Penrose')


# Final Estimation by Iterated Least Squares 
# h(x) = theta = arctan((x2-a2)/(x1-a1))
# H = (1/r)*[-sin(h) cos(h)]
count = 1
maxcount = 10
error = [100,100]
maxerror = 0.1
xhatold = xhat
# Using 3 measurements, 3 DOA's
theta =  np.array([[theta1], [theta2], [theta3]])
R =  np.array([[sigma1*sigma1, 0, 0], [0, sigma2*sigma2, 0], [0, 0, sigma3*sigma3]])
while (count < maxcount and (error[0]+error[1]/2) > maxerror):
    print("Iteration: ",count)
    #thetahat = atan2((xhatold[1]-a3[1]),(xhatold[0]-a3[0]))
    thetahat = np.array([
        [np.arctan2((xhatold[1]-a1[1]),(xhatold[0]-a1[0])).item()], 
        [np.arctan2((xhatold[1]-a2[1]),(xhatold[0]-a2[0])).item()],
        [np.arctan2((xhatold[1]-a3[1]),(xhatold[0]-a3[0])).item()]
        ])
        
    #H = (1/(((xhatold-a3)'*(xhatold-a3))^.5))*[-sin(thetahat) cos(thetahat)]
    H = np.array([
        [(1/np.power((xhatold-a1).conj().transpose()@(xhatold-a1),0.5))*np.array([-np.sin(thetahat[0].item()), np.cos(thetahat[0].item())])][0][0],
        [(1/np.power((xhatold-a2).conj().transpose()@(xhatold-a2),0.5))*np.array([-np.sin(thetahat[1].item()), np.cos(thetahat[1].item())])][0][0],
        [(1/np.power((xhatold-a3).conj().transpose()@(xhatold-a3),0.5))*np.array([-np.sin(thetahat[2].item()), np.cos(thetahat[2].item())])][0][0]
        ])

    P = np.linalg.inv(H.conj().transpose()@np.linalg.inv(R)@H)
    xhatnew = xhatold + P@H.conj().transpose()@np.linalg.inv(R)@(theta-thetahat)
    error = abs(xhatnew - xhatold)
    xhatold = xhatnew
    count = count+1
    print("Error: ",error)
    plt.plot(xhatnew[0], xhatnew[1], '*k', label="Estimate Update, Error={}".format(error))
    plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
    plt.tight_layout()
    plt.draw()
    input("Press Enter to continue...")
    print("")



# Plot estimated target location and EEP
plt.plot(xhatnew[0], xhatnew[1], '*g', label='Final Target Location Estimate')
ellipse = calculateEllipse(xhatnew, P) 
plt.plot((ellipse[:,0]+xhatnew[0]),(ellipse[:,1]+xhatnew[1]),'-g')
plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
input("Press Enter to continue...")


