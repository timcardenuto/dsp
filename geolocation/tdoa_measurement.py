import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import sys


plt.ion()
fig,axs = plt.subplots(1,1)


# Known sensor locations
A = np.array([1,1])
B = np.array([77,31])
dAB = np.sqrt(np.power(B[0]-A[0],2)+np.power(B[1]-A[1],2))
plt.plot(A[0], A[1], '*k')
plt.plot(B[0], B[1], '*k')

# make up a tB and tA for the scenario - in real life these would be the input measurements to this algorithm
# set tA = 0 as the reference point
# set tB such that tA < tB < dAB/c+tA  (physical constraint)
# SHIFT = variable you use in this script to shift how close the emitter will be to point A or B, affects the hyperbola (also called isochrone). A value of 0 would result in a perfect line, perpendicular to the line between A and B. A value of -1 or 1 will put the emitter on top of either point A or B respectively (no line no hyperbola). So this value should be -1 < SHIFT < 1. 
SHIFT = .5
c = 299792458  # speed of RF wave through vacuum (close enough)
tA = 0
tB = (dAB/c+tA)*SHIFT
TDOA = tB - tA
print("TDOA: "+str(TDOA))

# Initial guess at emitter location and time
# we don't know what point the emitter is at (x,y), and we don't have range measurements to the emitter (dAE or dBE)
# we do know that 1 of the possible positions of E lies on the line between A and B, specifically at the point where the wave could reach point A prior to point B with a time difference given by the TDOA measurement
# therefore *1* solution of dAE is: 
dAE = (dAB - c*(TDOA))/2
print("Initial Guess dAE: "+str(dAE))

# from this distance dAE, and the given points A, B and distance dAB, we can calculate *1* possible emitter point E
E = A + dAE*(B-A)/dAB 
Ex = E[0]
Ey = E[1]
print("Initial Guess E(x,y): "+str(Ex)+","+str(Ey))
plt.plot(Ex, Ey, '.r')

# now we can find the time that the emitter transmitted, if it was at this point E, such that the TDOA is correct
# this is only *1* of the (technically infinite) solutions
# using a definition of the distance dAE = c(tA - tE):
tE_guess = tA - dAE/c
print("Initial Guess tE: "+str(tE_guess))

# plt.draw()
# sys.exit(0)

for i in range(10):
    print("")
    print("Iteration: "+str(i))
    tE_guess = tE_guess + tE_guess
    dAE = c*(tA-tE_guess)
    dBE = c*(tB-tE_guess)
    print("tE_guess: "+str(tE_guess))
    print("dAE: "+str(dAE))
    print("dBE: "+str(dBE))

    # This is for testing, in real world E is unknown and dAE/dBE are calculated by measuring tA and tB, and then estimating tE iteratively from the basic case (on the line A->B)
    # E = np.array([32,41])
    # plt.plot(E[0], E[1], '^k')
    # dAE = np.sqrt(np.power(E[0]-A[0],2)+np.power(E[1]-A[1],2))  
    # dBE = np.sqrt(np.power(E[0]-B[0],2)+np.power(E[1]-B[1],2))

    # plot real triangle between sensors and emitter
    # pts = np.array([A, B, E])
    # p = Polygon(pts, closed=False)
    # ax = plt.gca()
    # ax.add_patch(p)

    # # My attempt to solve for real
    # # TODO rotate AB line to have y=0, can't just zero it without rotation first
    # xE = ((dAE*dAE) - (dBE*dBE) + (B[0]*B[0])) / (2*(B[0]-A[0]))
    # yE = np.sqrt((dBE*dBE) - np.power(xE-B[0],2)) + B[1]
    # # TODO rotate back to get 'real' xE,yE
    # print("xE: "+str(xE))
    # print("yE: "+str(yE))


    # cheating based on equations found online
    
    # Find distance from point A to an intermediate point (C) on line between points A and B that lies perpendicular to point E
    dAC = ((dAE*dAE) - (dBE*dBE) + (dAB*dAB))/(2*dAB)   
    
    # Find distance from intermediate point C to point E (creates right triangle, this variable typically called 'h' for height)
    dCE = np.sqrt((dAE*dAE) - (dAC*dAC))
    
    # Now find the coordinates for point C
    C = A + dAC*(B-A)/dAB  
    # print(C)
    # plt.plot(C[0], C[1], 'or')

    # Finally find the coordinates of the emitter for this tE - this is only *1* possible location! 
    # Oh and actually for each possible tE there's actually *2* symmetrical possible locations E
    # TODO tried to use matrices but .... needed the  values flipped and didn't work?
    # E = C + dCE*(np.flip(B)-np.flip(A))/dAB
    Ex1 = C[0] - dCE*(B[1]-A[1])/dAB
    Ey1 = C[1] + dCE*(B[0]-A[0])/dAB
    Ex2 = C[0] + dCE*(B[1]-A[1])/dAB
    Ey2 = C[1] - dCE*(B[0]-A[0])/dAB

    print("E1: "+str(Ex1)+","+str(Ey1))
    print("E2: "+str(Ex2)+","+str(Ey2))
    plt.plot(Ex1, Ey1, '.r')
    plt.plot(Ex2, Ey2, '.r')

    # ax.set_xlim(1,7)
    # ax.set_ylim(1,7)
    axs.set_aspect('equal')
    plt.draw()
    input("Press Enter to continue...")

sys.exit(0)

# uhhhh where is the TDOA here??

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

plt.plot(a1[0]/1852, a1[1]/1852, 'ob', label='Aircraft 1')
plt.plot(a2[0]/1852, a2[1]/1852, 'ob', label='Aircraft 2')
plt.plot(a3[0]/1852, a3[1]/1852, 'ob', label='Aircraft 3')
plt.plot(xhat[0]/1852, xhat[1]/1852, '*g', label='Estimated location')
 
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
