# Passive Geolocation Homework 7 Problem 3 - Tim Cardenuto
# TODO set plot axes to be based on min/max values that are being plotted. The matplotlib autosizing doesn't seem to work when I need it to be a box

import numpy as np
import matplotlib.pyplot as plt

plt.ion()
fig,axs = plt.subplots(1,1)
axs.set_xlabel('west/east (nm)')
axs.set_ylabel('south/north (nm)')
axs.set_title('DOA Geolocation Estimation Using Iterated Least Squares')


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



def estimateGeolocation(loc_array, theta_array, sigma_array):

    # Initial estimate of target position (initial guess) given by Moore-Penrose algorithm
    # u1 = np.vstack([[np.cos(theta1)], [np.sin(theta1)]])        # unit vector along DOA
    # M1 = np.identity(2) - u1@(u1.conj().transpose())            # Moore-Penrose matrix
    # u2 = np.vstack([[np.cos(theta2)], [np.sin(theta2)]])
    # M2 = np.identity(2) - u2@(u2.conj().transpose())
    # u3 = np.vstack([[np.cos(theta3)], [np.sin(theta3)]])
    # M3 = np.identity(2) - u3@(u3.conj().transpose())
    # xhat = np.linalg.inv(M1+M2+M3)@(M1@loc1+M2@loc2+M3@loc3)

    Msum = []
    MTasum = []
    for a,t in zip(loc_array,theta_array):
        u = np.vstack([[np.cos(t)], [np.sin(t)]])           # unit vector along DOA
        M = np.identity(2) - u@(u.conj().transpose())       # Moore-Penrose matrix
        Msum = M if len(Msum) == 0 else Msum + M
        MTasum = M@a if len(MTasum) == 0 else MTasum + M@a
    xhat = np.linalg.inv(Msum)@MTasum
    print(xhat)
    print("")

    plt.plot(xhat[0], xhat[1], '*r', label='Initial Target Location Estimate w/ Moore-Penrose')

    # Final estimation by Iterated Least Squares algorithm
    # theta =  np.vstack([        # measurements
    #     [theta1],
    #     [theta2],
    #     [theta3]
    #     ])
    # R =  np.vstack([            # covariance matrix (measurement errors)
    #     [sigma1*sigma1, 0, 0],
    #     [0, sigma2*sigma2, 0],
    #     [0, 0, sigma3*sigma3]
    #     ])

    # theta = np.vstack([[theta_array[0]]])
    # for t in theta_array[1:]:
        # theta = np.vstack([theta, [t]])

    # measurements
    theta = np.vstack(theta_array)

    # covariance matrix (measurement errors)
    R = np.identity(len(sigma_array))*np.vstack(np.hstack(sigma_array)*np.hstack(sigma_array))


    count = 1
    maxcount = 10               # limit number of runs w/ a maximum count
    error = np.hstack([100,100])
    maxerror = 0.1              # when to stop iterating (acceptable error)
    while (count < maxcount and max(error) > maxerror):
        print("Iteration: ",count)
        # General solution:
        #   h(x) = theta = arctan((x2-loc2)/(x1-loc1))
        #   H = (1/r)*[-sin(h) cos(h)]
        # MATLAB specific solution:
        #   thetahat = atan2((xhat(2)-a(2)),(xhat(1)-a(1)))
        #   H = (1/(((xhat-a)'*(xhat-a))^.5))*[-sin(thetahat) cos(thetahat)]
        #   P = inv(H'*inv(R)*H)
        #   xhat = xhat + P*H'*inv(R)*(theta-thetahat)

        # thetahat = np.vstack([
        #     [np.arctan2((xhat[1]-loc1[1]), (xhat[0]-loc1[0]))],
        #     [np.arctan2((xhat[1]-loc2[1]), (xhat[0]-loc2[0]))],
        #     [np.arctan2((xhat[1]-loc3[1]), (xhat[0]-loc3[0]))]
        #     ])

        xtemp = xhat - loc_array[0]
        thetahat = np.arctan2(xtemp[1], xtemp[0])
        for a in loc_array[1:]:
            xtemp = xhat - a
            thetahat = np.vstack([thetahat, [np.arctan2(xtemp[1], xtemp[0])]])

        # H = np.vstack([
        #     (1/np.power((xhat-loc1).conj().transpose()@(xhat-loc1),0.5))*np.hstack([-np.sin(thetahat[0]), np.cos(thetahat[0])]),
        #     (1/np.power((xhat-loc2).conj().transpose()@(xhat-loc2),0.5))*np.hstack([-np.sin(thetahat[1]), np.cos(thetahat[1])]),
        #     (1/np.power((xhat-loc3).conj().transpose()@(xhat-loc3),0.5))*np.hstack([-np.sin(thetahat[2]), np.cos(thetahat[2])])
        #     ])

        H = (1/np.power((xhat-loc_array[0]).conj().transpose()@(xhat-loc_array[0]),0.5))*np.hstack([-np.sin(thetahat[0]), np.cos(thetahat[0])])
        for a,th in zip(loc_array[1:],thetahat[1:]):
            H = np.vstack([H, (1/np.power((xhat-a).conj().transpose()@(xhat-a),0.5))*np.hstack([-np.sin(th), np.cos(th)])])

        P = np.linalg.inv(H.conj().transpose()@np.linalg.inv(R)@H)

        xhatnew = xhat + P@H.conj().transpose()@np.linalg.inv(R)@(theta-thetahat)
        error = abs(xhatnew - xhat)
        xhat = xhatnew
        count = count+1

        print("Error: ",error)
        axs.plot(xhat[0], xhat[1], '*k', label="Estimate Update, Error={}".format(error))
        axs.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
        # plt.tight_layout()
        axs.set(xlim=(-5, 100), ylim=(-5, 100))
        axs.set_aspect('equal')
        plt.draw()
        #input("Press Enter to continue...")
        print("")

    return xhat,P



# Sensor 1
loc1 = np.vstack([[0],[0]])                                 # Sensor location x,y (meters)
theta1 = 50/180*np.pi                                       # DOA angle measurement (radians)
sigma1 = 3*np.pi/180                                        # DOA angle error, i.e. standard deviation (radians)
axs.plot(loc1[0], loc1[1], 'ok', label='Sensor 1')
axs.plot([loc1[0], 100*np.cos(theta1)+loc1[0]],[loc1[1], 100*np.sin(theta1)+loc1[1]], 'k')

# Sensor 2
loc2 = np.vstack([[20],[0]])
theta2 = 70/180*np.pi
sigma2 = 2*np.pi/180
axs.plot(loc2[0], loc2[1], 'ok', label='Sensor 2')
axs.plot([loc2[0], 100*np.cos(theta2)+loc2[0]],[loc2[1], 100*np.sin(theta2)+loc2[1]], 'k')

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

# Use Moore-Penrose and Iterated Least Squares (ILS) to estimate the geolocation and covariance matrix
xhat,P = estimateGeolocation(loc_array, theta_array, sigma_array)

# Calculate the Elliptical Error Probable (EEP), i.e. error ellipse
ellipse = calculateEllipse(xhat, P, 0.95)

# Plot estimated target location and EEP
axs.plot(xhat[0], xhat[1], '*g', label='Final Target Location Estimate')
axs.plot(ellipse[:,0],ellipse[:,1],'-g')
axs.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
axs.set(xlim=(-5, 100), ylim=(-5, 100))
axs.set_aspect('equal')
input("Press Enter to continue...")
