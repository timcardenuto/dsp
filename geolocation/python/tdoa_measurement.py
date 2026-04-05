import numpy as np
import random as rand
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import sys
from geolocation import *

# constant variables
c = 299792458   # speed of RF wave through vacuum (close enough)
nm = 1852       # 1 nautical mile = 1852 meters

# set up plotting
plt.ion()
plt.xlabel('y (meters)')
plt.ylabel('x (meters)')
plt.title("TDOA Hyperbolas")


def get_tdoa_hyperbola(locA, toaA, locB, toaB, hyperbola_length):
    """
    locA = sensor A position [x y]
    toaA = time of arrival (TOA) of the signal at sensor A
    locB = sensor B position [x y]
    toaB = time of arrival (TOA) of the signal at sensor B
    hyperbola_length = length of the hyperbola to draw (for visual representation only, does not affect the measurement)
    """

    print("")
    print("######################################")
    print(">> get_tdoa_hyperbola")
    print("######################################")
    print("locA: "+str(locA))
    print("toaA: "+str(toaA))
    print("locB: "+str(locB))
    print("toaB: "+str(toaB))
    print("hyperbola_length: "+str(hyperbola_length))

    TDOA = toaB - toaA
    print("TDOA: "+str(TDOA))
    print("Range difference: "+str(c*TDOA))

    # Initial guess at emitter location and time
    # We don't know what point the emitter is at (x,y), and we don't have range measurements to the emitter (dAE or dBE)
    # We do know that *1* of the possible positions of E lies between A and B, specifically at the point where the wave could reach point A prior to point B with a time difference given by the TDOA measurement
    # Since the wave from E must reach either A or B first, there is always a valid solution that lies between them providing the time difference (in the extreme case, E is directly on the far side of A or B on a line connecting all 3 points, but what this would look like from a Time of Arrival (TOA) perspective is that the wave arrives a A (or B in reverse) first and then B a the speed of light time afterwards which would make it look like E is *at* A). 
    # Therefore *1* solution of dAE is: 
    dAB = np.sqrt(np.power(locB[0]-locA[0],2)+np.power(locB[1]-locA[1],2))
    dAE = (dAB - c*(TDOA))/2
    print("Initial Guess dAE: "+str(dAE))

    # from this distance dAE, and the given points A, B and distance dAB, we can calculate *1* possible emitter point E
    locE = locA + dAE*(locB-locA)/dAB 
    print("Initial Guess E(x,y): "+str(locE[0])+","+str(locE[1]))

    # now we can find the time that the emitter transmitted, if it was at this point E, such that the TDOA is correct
    # this is only *1* of the (technically infinite) solutions
    # using a definition of the distance dAE = c(toaA - toaE):
    toaE_guess = toaA - dAE/c
    print("Initial Guess toaE: "+str(toaE_guess))

    # Batch processor
    # explanation for math steps given in the iterative version of this below (same math just w/o matrices)
    resolution = 100*hyperbola_length      # determines how many points along the line to use, this is just for drawing purposes. 10x seems to be reasonable
    toaE = np.vstack(np.linspace(toaE_guess, (toaE_guess*hyperbola_length), resolution))
    dAE = c*(toaA-toaE)
    dBE = c*(toaB-toaE)
    dAC = ((dAE*dAE) - (dBE*dBE) + (dAB*dAB))/(2*dAB)
    dCE = np.sqrt((dAE*dAE) - (dAC*dAC))

    print("toaE: "+str(toaE[0]))
    print("dAE: "+str(dAE[0]))
    print("dBE: "+str(dBE[0]))
    print("dAC: "+str(dAC[0]))
    print("dCE: "+str(dCE[0]))

    locC = locA + dAC*(locB-locA)/dAB

    # Ep(x,y) possible emitter location positive side
    Exp = np.vstack(locC[:,0]) - dCE*(locB[1]-locA[1])/dAB
    Eyp = np.vstack(locC[:,1]) + dCE*(locB[0]-locA[0])/dAB
    # En(x,y) possible emitter location negative side
    Exn = np.vstack(locC[:,0]) + dCE*(locB[1]-locA[1])/dAB
    Eyn = np.vstack(locC[:,1]) - dCE*(locB[0]-locA[0])/dAB

    locE = np.hstack([np.vstack([np.flip(Exp), Exn]), np.vstack([np.flip(Eyp), Eyn])])
    print("locE: "+str(locE[0]))
    print("")
    return TDOA,locE



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
    t1 = te + d1e/c + np.random.normal(0, sigma_time) # rand.normalvariate(mu=0.0, sigma=sigma_time)    # realistic t1 has to include error, but couldn't add it to the other terms above
    t2 = te + d2e/c + np.random.normal(0, sigma_time) # rand.normalvariate(mu=0.0, sigma=sigma_time)
    t3 = te + d3e/c + np.random.normal(0, sigma_time) # rand.normalvariate(mu=0.0, sigma=sigma_time)


    length = 10  # how long the line should be, in terms of the # of toaE's lengths..... hard to quantify but 10 gives a good picture for 2 points, would need to calculate this based on the overall scenario map limits
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
    plt.plot(hyperbola32[:,0], hyperbola32[:,1], '-k', label='S2>S3 TDOA {}'.format(tdoa32))
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
    # z = np.array([7298.46419703, 5368.64643455])
    # z = np.array([8000, -50])     # was uncommented ...
    # z = np.array([(c*tdoa21),(c*tdoa31)])
    # z = np.array([13000,-5600])

    # TDOA measurement error (std deviation), in terms of distance (m/s * s) = m
    sigma_dist = c*sigma_time
    

    #########################################
    # Using the Geolocation library function
    request = [
        {'mtype': 'tdoa_range', 'value':  z[0], 'sigma': sigma_dist, 's1loc': s1, 's2loc': s2},
        {'mtype': 'tdoa_range', 'value':  z[1], 'sigma': sigma_dist, 's1loc': s1, 's2loc': s3}
        ]
    ellipse = geolocate(request, 0.95)
    plt.plot(ellipse['x'], ellipse['y'], '*g', label='Final Target Location Estimate')
    plt.plot(ellipse['shape'][:,0],ellipse['shape'][:,1],'-g', label='Error Ellipse')
    plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
    plt.tight_layout()
    plt.draw()

    input("Press Enter to continue...")
    
    #########################################
    # Original code here

    # Covariance matrix
    R = (sigma_dist*sigma_dist)*np.identity(2)
    
    # Initial estimate of target position (initial guess) given by Moore-Penrose algorithm
    # NOTE: This does work, but see below note. This initial guess is based on measurements which have random error, and depending on the result here there are errors resolving the geo.
    # TODO: Which or how many of the sensor locations should I use here? The measurements don't align 1-to-1 with sensor locations ... 
    #       if I'm only use two range measurements but they are related to measurements where I used all 3 sensors???
    #       This algorithm works regardless of whether you use s1,s2 or s1,s2,s3; and either way sometimes produces singular matrix errors later 
    Msum = []
    MTasum = []
    for a,t in zip([s1,s2],z):
        u = np.vstack([[np.cos(t)], [np.sin(t)]])           # unit vector along TDOA?
        M = np.identity(2) - u@(u.conj().transpose())       # Moore-Penrose matrix
        Msum = M if len(Msum) == 0 else Msum + M
        MTasum = M@a if len(MTasum) == 0 else MTasum + M@a
    xhat = np.linalg.pinv(Msum)@MTasum

    # TODO apparently Moore-Penrose doesn't work for ranges so this is supposed to but it doesn't here, need to debug...
    # uniquelocs = np.concatenate((s1, s2, s3))
    # xhat = np.mean(uniquelocs, axis=0)

    print("xhat: "+str(xhat))
    print("")

    # TODO: When you comment out the normalvariate error for the time measurements, and use Moore-Penrose above this is what you get every time 
    #       When you include the error, this changes every time (as expected) but sometimes it gives Singular Matrix errors... Sometimes it converges on one or the other cross over points
    # NOTE: Even though I *have definitely seen this provide the correct result multiple times* ... now it's not? So there's still some randomness below that determines which cross over point it converges to?
    # xhat = np.array([251.44757026, 31248.34867717])     # initial guess at location of target that seems to convergense on the wrong side...
    # xhat = np.array([-1895.33108933, 43617.27057437])   # another guess that **should** result in convergense on the right side... I saw it!

    a1 = s1 
    a2 = s2
    a3 = s3
    print("sigma: "+str(sigma_dist))
    print("R: "+str(R))
    print("z: "+str(z))
    print("a1:"+str(a1))
    print("a2:"+str(a2))
    print("a3:"+str(a3))
    print("\n")


    sref = [a1, a1]
    srel = [a2, a3]

    count = 1
    maxcount = 10
    error = [100,100]
    maxerror = .1
    while (count < maxcount and (error[0]+error[1]/2) > maxerror):

        # Range h and H
        # h = r = ((x-a)'*(x-a))^.5
        # H = (1/r)*(x-a)'

        h = np.zeros(len(sref))
        H = np.zeros((len(sref),2))
        # TODO: Probably a way to do this w/o loop using matrices but too lazy to figure it out...
        for i in range(len(sref)):
            h[i] = (np.power(np.transpose(xhat-sref[i])@(xhat-sref[i]),0.5) - np.power(np.transpose(xhat-srel[i])@(xhat-srel[i]),0.5))
            H[i] = ((1/np.power(np.transpose(xhat-sref[i])@(xhat-sref[i]),0.5))*np.transpose(xhat-sref[i]) - (1/np.power(np.transpose(xhat-srel[i])@(xhat-srel[i]),0.5))*np.transpose(xhat-srel[i]))

        # Old way, had to manually add every sensor...
        # h = np.array([
        #     np.power(np.transpose(xhat-a1)@(xhat-a1),0.5) - np.power(np.transpose(xhat-a2)@(xhat-a2),0.5), 
        #     np.power(np.transpose(xhat-a1)@(xhat-a1),0.5) - np.power(np.transpose(xhat-a3)@(xhat-a3),0.5)
        #     ])
        # H = np.array([
        #     (1/np.power(np.transpose(xhat-a1)@(xhat-a1),0.5))*np.transpose(xhat-a1) - (1/np.power(np.transpose(xhat-a2)@(xhat-a2),0.5))*np.transpose(xhat-a2), 
        #     (1/np.power(np.transpose(xhat-a1)@(xhat-a1),0.5))*np.transpose(xhat-a1) - (1/np.power(np.transpose(xhat-a3)@(xhat-a3),0.5))*np.transpose(xhat-a3)
        #     ])

        # print("h: "+str(h))
        # print("H: "+str(H))
        # print("")

        P = np.linalg.pinv(np.transpose(H)@(np.linalg.pinv(R))@(H))

        xhatnew = xhat + P@(np.transpose(H))@(np.linalg.pinv(R))@(z-h)
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
    ellipse = ellipse@(R)

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

    locA = s1
    locB = s2
    toaA = t1
    toaB = t2

    # # TEST ONLY - make up a toaB and toaA for the scenario - in real life these would be input measurements from a prior stage in the receivers
    # # this won't work for more than 2 points
    # # set toaA = 0 as the reference point
    # # set toaB such that toaA < toaB < dAB/c+toaA  (physical constraint)
    # # SHIFT = variable you use in this script to shift how close the emitter will be to point A or B, affects the hyperbola (also called isochrone). A value of 0 would result in a perfect line, perpendicular to the line between A and B. A value of -1 or 1 will put the emitter on top of either point A or B respectively (no line no hyperbola). So this value should be -1 < SHIFT < 1. 
    # dAB = np.sqrt(np.power(locB[0]-locA[0],2)+np.power(locB[1]-locA[1],2))
    # SHIFT = 0.75
    # toaA = 0
    # toaB = (dAB/c+toaA)*SHIFT

    TDOA = toaB - toaA
    print("TDOA: "+str(TDOA))
    print("Range difference: "+str(c*TDOA))

    dAB = np.sqrt(np.power(locB[0]-locA[0],2)+np.power(locB[1]-locA[1],2))
    dAE = (dAB - c*(TDOA))/2
    print("Initial Guess dAE: "+str(dAE))

    # from this distance dAE, and the given points A, B and distance dAB, we can calculate *1* possible emitter point E
    locE = locA + dAE*(locB-locA)/dAB 
    print("Initial Guess E(x,y): "+str(locE[0])+","+str(locE[1]))

    # now we can find the time that the emitter transmitted, if it was at this point E, such that the TDOA is correct
    # this is only *1* of the (technically infinite) solutions
    # using a definition of the distance dAE = c(toaA - toaE):
    toaE_guess = toaA - dAE/c
    print("Initial Guess toaE: "+str(toaE_guess))

    # This would be a interative way to estimate, instead of batched (matrix)
    # as you step through the iteration you can see the red estimated points sit on the same line as the batched process above. 
    for i in range(4):
        print("")
        print("Iteration: "+str(i))
        toaE_guess = toaE_guess + toaE_guess
        dAE = c*(toaA-toaE_guess)
        dBE = c*(toaB-toaE_guess)
        print("toaE_guess: "+str(toaE_guess))
        print("dAE: "+str(dAE))
        print("dBE: "+str(dBE))

        # This is for testing, in real world E is unknown and dAE/dBE are calculated by measuring toaA and toaB, and then estimating toaE iteratively from the basic case (on the line A->B)
        # locE = np.array([32,41])
        # plt.plot(locE[0], locE[1], '^k')
        # dAE = np.sqrt(np.power(locE[0]-locA[0],2)+np.power(locE[1]-locA[1],2))  
        # dBE = np.sqrt(np.power(locE[0]-locB[0],2)+np.power(locE[1]-locB[1],2))

        # plot real triangle between sensors and emitter
        # pts = np.array([locA, locB, locE])
        # p = Polygon(pts, closed=False)
        # ax = plt.gca()
        # ax.add_patch(p)

        # # My attempt to solve using the assumption that Y axis = zero. This is a common math problem but not the generic case for all possible x,y coordinates. Doing it this way would require some kind of rotation from original x,y coordinates to a y=0 and then back.
        # # TODO rotate AB line to have y=0, can't just zero it without rotation first
        # xE = ((dAE*dAE) - (dBE*dBE) + (locB[0]*locB[0])) / (2*(locB[0]-locA[0]))
        # yE = np.sqrt((dBE*dBE) - np.power(xE-locB[0],2)) + locB[1]
        # # TODO rotate back to get 'real' xE,yE
        # print("xE: "+str(xE))
        # print("yE: "+str(yE))


        # cheating based on equations found online
        
        # Find distance from point A to an intermediate point (C) on line between points A and locB that lies perpendicular to point E
        dAC = ((dAE*dAE) - (dBE*dBE) + (dAB*dAB))/(2*dAB)   
        
        # Find distance from intermediate point C to point E (creates right triangle, this variable typically called 'h' for height)
        dCE = np.sqrt((dAE*dAE) - (dAC*dAC))
        
        # Now find the coordinates for point C
        locC = locA + dAC*(locB-locA)/dAB  
        # print(locC)
        # plt.plot(locC[0], locC[1], 'or')

        # Finally find the coordinates of the emitter for this toaE - this is only *1* possible location! 
        # Oh and actually for each possible toaE there's actually *2* symmetrical possible locations E
        # TODO tried to use matrices but .... needed the  values flipped and didn't work?
        # locE = locC + dCE*(np.flip(locB)-np.flip(locA))/dAB
        Exp = locC[0] - dCE*(locB[1]-locA[1])/dAB
        Eyp = locC[1] + dCE*(locB[0]-locA[0])/dAB
        Exn = locC[0] + dCE*(locB[1]-locA[1])/dAB
        Eyn = locC[1] - dCE*(locB[0]-locA[0])/dAB

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
