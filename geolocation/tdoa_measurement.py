import numpy as np
import random as rand
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import sys

# constant variables
c = 299792458   # speed of RF wave through vacuum (close enough)
nm = 1852       # 1 nautical mile = 1852 meters

# set up plotting
plt.ion()
plt.xlabel('y (meters)')
plt.ylabel('x (meters)')
plt.title("TDOA Hyperbolas")


def get_tdoa_hyperbola(A, tA, B, tB, length):

    TDOA = tB - tA
    print("TDOA: "+str(TDOA))
    print("Range difference: "+str(c*TDOA))

    # Initial guess at emitter location and time
    # we don't know what point the emitter is at (x,y), and we don't have range measurements to the emitter (dAE or dBE)
    # we do know that 1 of the possible positions of E lies on the line between A and B, specifically at the point where the wave could reach point A prior to point B with a time difference given by the TDOA measurement
    # therefore *1* solution of dAE is: 
    dAB = np.sqrt(np.power(B[0]-A[0],2)+np.power(B[1]-A[1],2))
    dAE = (dAB - c*(TDOA))/2
    print("Initial Guess dAE: "+str(dAE))

    # from this distance dAE, and the given points A, B and distance dAB, we can calculate *1* possible emitter point E
    E = A + dAE*(B-A)/dAB 
    print("Initial Guess E(x,y): "+str(E[0])+","+str(E[1]))

    # now we can find the time that the emitter transmitted, if it was at this point E, such that the TDOA is correct
    # this is only *1* of the (technically infinite) solutions
    # using a definition of the distance dAE = c(tA - tE):
    tE_guess = tA - dAE/c
    print("Initial Guess tE: "+str(tE_guess))

    # Batch processor
    # explanation for math steps given in the iterative version of this below (same math just w/o matrices)
    resolution = 100*length      # determines how many points along the line to use, 10x seems to be reasonable
    tE = np.vstack(np.linspace(tE_guess, (tE_guess*length), resolution))
    dAE = c*(tA-tE)
    dBE = c*(tB-tE)
    dAC = ((dAE*dAE) - (dBE*dBE) + (dAB*dAB))/(2*dAB)   
    dCE = np.sqrt((dAE*dAE) - (dAC*dAC))
    C = A + dAC*(B-A)/dAB

    # Ep(x,y) possible emitter location positive side
    Exp = np.vstack(C[:,0]) - dCE*(B[1]-A[1])/dAB
    Eyp = np.vstack(C[:,1]) + dCE*(B[0]-A[0])/dAB
    # En(x,y) possible emitter location negative side
    Exn = np.vstack(C[:,0]) + dCE*(B[1]-A[1])/dAB
    Eyn = np.vstack(C[:,1]) - dCE*(B[0]-A[0])/dAB

    E = np.hstack([np.vstack([np.flip(Exp), Exn]), np.vstack([np.flip(Eyp), Eyn])])
    print(E)
    print("")
    return TDOA,E



if __name__ == "__main__":

    # Known sensor locations
    s1 = np.array([-10*nm, 20*nm])
    s2 = np.array([0, 20*nm])
    s3 = np.array([10*nm, 20*nm])

    # target emitter (truth)
    e1 = np.array([3*nm, 40*nm])

    # calculate times of arrival (TOA) at the sensors and add error
    d1e = np.sqrt(np.power(s1[0]-e1[0],2)+np.power(s1[1]-e1[1],2))
    d2e = np.sqrt(np.power(s2[0]-e1[0],2)+np.power(s2[1]-e1[1],2))
    d3e = np.sqrt(np.power(s3[0]-e1[0],2)+np.power(s3[1]-e1[1],2))
    # TODO what is a realistic error prob distribution? Using the following
    # NOTE assume all sensors have same variance for timing error (1 us), not necessarily same actual error value
    sigma_time = 0.000001
    t1 = 0              # TOA at sensor 1 is our reference point
    te = t1 - d1e/c     # emission time will be negative since we're using TOA at sensor 1 for the zero reference point
    t2 = te + d2e/c + rand.normalvariate(mu=0.0, sigma=sigma_time)
    t3 = te + d3e/c + rand.normalvariate(mu=0.0, sigma=sigma_time)
    t1 = te + d1e/c + rand.normalvariate(mu=0.0, sigma=sigma_time)    # realistic t1 has to include error, but couldn't add it to the other terms above

    length = 10  # how long the line should be, in terms of the # of tE's lengths..... hard to quantify but 10 gives a good picture for 2 points, would need to calculate this based on the overall scenario map limits
    tdoa21,hyperbola21 = get_tdoa_hyperbola(s1, t1, s2, t2, length)
    tdoa31,hyperbola31 = get_tdoa_hyperbola(s1, t1, s3, t3, length)
    tdoa32,hyperbola32 = get_tdoa_hyperbola(s2, t2, s3, t3, length)

    # plot stuff
    plt.plot(e1[0], e1[1], '*r', label='True Emitter')
    plt.plot(s1[0], s1[1], '^b', label='Sensor 1')
    plt.plot(s2[0], s2[1], '^k', label='Sensor 2')
    plt.plot(s3[0], s3[1], '^g', label='Sensor 3')
    plt.plot(hyperbola21[:,0], hyperbola21[:,1], '-k', label='S1>S2 TDOA {}'.format(tdoa21))
    plt.plot(hyperbola31[:,0], hyperbola31[:,1], '-k', label='S1>S3 TDOA {}'.format(tdoa31))
    plt.plot(hyperbola32[:,0], hyperbola32[:,1], '-k', label='S1>S3 TDOA {}'.format(tdoa32))
    plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
    plt.tight_layout()
    plt.draw()


    input("Press Enter to continue...")


    #### Geolocation ####
    print("\n##############################")
    print("Geolocation")
    print("##############################\n")

    # TDOA measurements, in terms of "range differences"
    # r1 - r2 = c*(ta1-ta2)
    # r1 - r3 = c*(ta1-ta3)
    z = np.array([(c*np.abs(tdoa21)),(c*np.abs(tdoa31))])
    # z = np.array([(c*tdoa21),(c*tdoa31)])

    # TDOA measurement error (std deviation), in terms of distance (m/s * s) = m
    sigma_dist = c*sigma_time
    R = (sigma_dist*sigma_dist)*np.identity(2)

    xhat = np.array([0,0])          # initial guess at location of target, in real life would need to be smarter about it
    a1 = s1 
    a2 = s2
    a3 = s3
    print("sigma: "+str(sigma_dist))
    print("z: "+str(z))
    print("a1:"+str(a1))
    print("a2:"+str(a2))
    print("a3:"+str(a3))
    print("\n")

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
    plt.plot(xhat[0], xhat[1], '*g', label='Estimated location')
    plt.plot((ellipse[:,0]+xhat[0]),(ellipse[:,1]+xhat[1]),'-', label='Error Ellipse')
    plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
    plt.tight_layout()
    plt.draw()

    input("Press Enter to continue...")



    #### Iterative TDOA hyperbola calc ####
    # comment out to see the iterative version overlayed (only for 1 hyperbola)
    sys.exit()

    A = s1
    B = s2
    tA = t1
    tB = t2

    # # TEST ONLY - make up a tB and tA for the scenario - in real life these would be input measurements from a prior stage in the receivers
    # # this won't work for more than 2 points
    # # set tA = 0 as the reference point
    # # set tB such that tA < tB < dAB/c+tA  (physical constraint)
    # # SHIFT = variable you use in this script to shift how close the emitter will be to point A or B, affects the hyperbola (also called isochrone). A value of 0 would result in a perfect line, perpendicular to the line between A and B. A value of -1 or 1 will put the emitter on top of either point A or B respectively (no line no hyperbola). So this value should be -1 < SHIFT < 1. 
    # dAB = np.sqrt(np.power(B[0]-A[0],2)+np.power(B[1]-A[1],2))
    # SHIFT = 0.75
    # tA = 0
    # tB = (dAB/c+tA)*SHIFT

    TDOA = tB - tA
    print("TDOA: "+str(TDOA))
    print("Range difference: "+str(c*TDOA))

    dAB = np.sqrt(np.power(B[0]-A[0],2)+np.power(B[1]-A[1],2))
    dAE = (dAB - c*(TDOA))/2
    print("Initial Guess dAE: "+str(dAE))

    # from this distance dAE, and the given points A, B and distance dAB, we can calculate *1* possible emitter point E
    E = A + dAE*(B-A)/dAB 
    print("Initial Guess E(x,y): "+str(E[0])+","+str(E[1]))

    # now we can find the time that the emitter transmitted, if it was at this point E, such that the TDOA is correct
    # this is only *1* of the (technically infinite) solutions
    # using a definition of the distance dAE = c(tA - tE):
    tE_guess = tA - dAE/c
    print("Initial Guess tE: "+str(tE_guess))

    # This would be a interative way to estimate, instead of batched (matrix)
    # as you step through the iteration you can see the red estimated points sit on the same line as the batched process above. 
    for i in range(4):
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

        # # My attempt to solve using the assumption that Y axis = zero. This is a common math problem but not the generic case for all possible x,y coordinates. Doing it this way would require some kind of rotation from original x,y coordinates to a y=0 and then back.
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
        Exp = C[0] - dCE*(B[1]-A[1])/dAB
        Eyp = C[1] + dCE*(B[0]-A[0])/dAB
        Exn = C[0] + dCE*(B[1]-A[1])/dAB
        Eyn = C[1] - dCE*(B[0]-A[0])/dAB

        print("E(x,y) positive: "+str(Exp)+","+str(Eyp))
        print("E(x,y) negative: "+str(Exn)+","+str(Eyn))
        plt.plot(Exp, Eyp, '.r')
        plt.plot(Exn, Eyn, '.r')

        # ax.set_xlim(1,7)
        # ax.set_ylim(1,7)
        # axs.set_aspect('equal')

        plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
        plt.tight_layout()
        plt.draw()
        input("Press Enter to continue...")