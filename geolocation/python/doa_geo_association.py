# Passive Geolocation Homework9 Problem 1 - Tim Cardenuto

import numpy as np
import matplotlib.pyplot as plt


# NOTE: sqrt in matlab will accept negative and complex numbers, the equivalent to numpy sqrt in matlab is realsqrt. 
# TODO: can the eigenvalues ever be negative or does that not make sense for these cases? Does that imply that there is an error farther up, like these measurements cannot be part of the same target?
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


# Import data and assumptions
#m = matfile('data_association_hw')
import scipy.io
m = scipy.io.loadmat('data_association_hw.mat')
#print(m)
xhats = np.array([m['xhat1'], m['xhat2'], m['xhat3']])  # apriori target locations
Ps = np.array([m['P1'], m['P2'], m['P3']])              # apriori target covariances
ac = m['ac']
z = (m['z']).item()/180*np.pi           # new measurement
sigma = (m['sigma']).item()/180*np.pi   # new measurement error
R = sigma*sigma
B=.05;
a=.95;
k=3.8415;
largestLL=0;


# Plot sensor location and measurements
plt.plot(ac[0][0], ac[1][0], 'ok', label='Sensor Location')
plt.plot([ac[0][0], 10*np.cos(z)], [ac[1][0], 10*np.sin(z)], 'k', label='DOA Measurement')
plt.plot([ac[0][0], 10*np.cos(z+sigma)], [ac[1][0], 10*np.sin(z+sigma)], '--k', label='DOA Error')
plt.plot([ac[0][0], 10*np.cos(z-sigma)], [ac[1][0], 10*np.sin(z-sigma)], '--k', label='DOA Error')

# Plot the apriori target estimated locations
plt.plot(xhats[:,0], xhats[:,1], '*r', label="A-Priori Target Locations")

# For each apriori target, determine which is most likely associated with the new measurement
for i,val in enumerate(xhats):
    print("Target ",i)
    xhat = xhats[i]     # Estimated location for this target
    P = Ps[i]           # Covariance for this estimated target location
    
    # Plot the apriori target EEP
    ellipse = calculateEllipse(xhat, P)
    plt.plot((ellipse[:,0]+xhat[0]),(ellipse[:,1]+xhat[1]),'-r')

    # Calculate ?
    zhat = np.arctan2(xhat[1]-ac[1], xhat[0]-ac[0]).item()  # I think this is always suppose to be scalar, so need the .item(), otherwise H ends up as vertical array and the H*P*H' matrix manipulation will fail
    H = (1/np.power((xhat-ac).conj().transpose()@(xhat-ac),0.5)) * (np.array([-np.sin(zhat), np.cos(zhat)])) # I think the first half of this before the multiply always ends up as a scalar, therefore have to use regular multiply not dot multiply (matlab must just do this conversion implicitly)
    md = np.power((z-zhat),2) / (H@P@(H.conj().transpose())+R)
    print(" md = ",md)
    if md < k:
        LL = -0.5*(np.log(2*np.pi)+np.log(H@P@(H.conj().transpose())+R)+md)
        print(" LL = ",LL)
        if LL > largestLL:
            largestLL = LL
            target = i
            targetxhat = xhat
            targetP = P
            targetzhat = zhat
            targetH = H
            
    print("")
    #i=i+2
    #j=j+4

print("Target ",target," has the biggest LL, is most likely source")
print(" xhat = ",targetxhat)
print(" P = ",targetP)
print(" zhat = ",targetzhat)
print(" H = ",targetH)
print("")

# Plot ellipse and most likely target 
ellipse = calculateEllipse(targetxhat, targetP)
plt.plot(targetxhat[0], targetxhat[1], '.b')
plt.plot((ellipse[:,0]+targetxhat[0]),(ellipse[:,1]+targetxhat[1]),'-b', label='Associated Target')

# Add new measurement to most likely target
# TODO ended up with a negative min(eigenvalue) which you can't take a sqrt of to get the semiminor axis... stuck an abs() on that but I feel like something is wrong, like can we not have negative eigenvalues? This also doesn't match the picture in my original homework...
xhat = targetxhat
h = targetzhat
H = targetH
P = targetP
K = P @ (H.conj().transpose()) @ np.linalg.inv(H@P@(H.conj().transpose())+R)
xhat = xhat + K*(z-h)
P = (1 - K@H) @ P
print("xhat = ",xhat)
print("P = ",P)
ellipse = calculateEllipse(xhat, P)
plt.plot(xhat[0], xhat[1], '.g', label='New Target Location')
plt.plot((ellipse[:,0]+xhat[0]),(ellipse[:,1]+xhat[1]),'-g', label='New Target EEP')

plt.xlabel('west/east (nm)')
plt.ylabel('south/north (nm)')
plt.title('DOA Target Association')
plt.legend(loc='upper left', shadow=True)
plt.show()

