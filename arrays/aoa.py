#!/usr/bin/env python

import scipy as sp
import matplotlib.pyplot as plt


"""
    num_elements
    type
    seperation
    ref_element
"""
def buildAntennaArray(num_elements,ref_element,seperation,type='linear'):
    if type == 'linear':
        array = sp.zeros(num_elements, dtype = sp.double)
        for idx in range(ref_element,num_elements,1):     # calculates distances on right side of reference element
            array[idx] = array[idx-1] + seperation
        for idx in range(ref_element-2,-1,-1):              # calculates distances on left side of reference element
            array[idx] = array[idx+1] + seperation
        return array
    elif type == 'circular':
        return 'not supported'
    else:
        return 'not supported'

""" Calculate phase delay using equation here
    (http://www.radartutorial.eu/06.antennas/Phased%20Array%20Antenna.en.html)
    This works for both calculating phase delay needed to TX at a certain angle,
    and calculating phase delay of an incoming signal angle
    array -
    angle - relative to boresight of the reference element, in radians. 0 degrees is boresight, and supports from left to right angles of incidence -90:0:90 degrees.
    wavelength - of carrier frequency in meters.
"""
def calculatePhaseDelay(array, angle, wavelength):
    phasedelay = 2*sp.pi*array*sp.sin(angle)/wavelength
    phasedelay = abs(phasedelay)%(2*sp.pi)*phasedelay/abs(phasedelay)   # this adds realism, you'd never be able to detect the overall phasedelays beyond the -pi to pi range
    return phasedelay

""" Based on wikipedia article for Phase Interferometry
    NOTE: this is on the receiver end, the tricky part is handling the ambiguities for elements beyond 1 wavelength
          I'm kind of cheating here? because I'm assuming the second element (first after the reference) is always correct
          and I'm just correcting the others based on it.... but maybe not, b/c this might just lock them all together,
          meaning the first one could be wrong, but once adjusted the others are all adjusted relative to it...
"""
def calculateAOA(array, phasedelay, wavelength):
    # so, we know the first two elements closest to the reference element will always be correct
    # (assuming array elements are spaced 1/2 wavelength of the signal being detected)
    # begging the question of whether we need the other elements....
    # but for now we'll use that to know that the others need wavelengths added

    # remove phase ambiguities
    # TODO these corrections are only valid for linear array, and with array[0] as the reference element
    aoa1 = (sp.arcsin((wavelength*phasedelay[1])/(2*sp.pi*array[1])))*180/sp.pi # calculate aoa of first and second element first, as a reference point
    print aoa1
    if aoa1 > -90 and aoa1 <= -42:
        phasedelay[3] = phasedelay[3]-(1*(2*sp.pi))    # [-89:-42]
    elif aoa1 > -42 and aoa1 <= 41:
        phasedelay[3] = phasedelay[3]-(0*(2*sp.pi))    # [-41:41]
    elif aoa1 > 41 and aoa1 < 90:
        phasedelay[3] = phasedelay[3]+(1*(2*sp.pi))    # [42:89]
    else:
        print "phase delay out-of-bounds"

    if aoa1 > -90 and aoa1 <= -30:
        phasedelay[4] = phasedelay[4]-(1*(2*sp.pi))    # [-89:-31]
    elif aoa1 > -30 and aoa1 <= 30:
        phasedelay[4] = phasedelay[4]-(0*(2*sp.pi))    # [-30:30]
    elif aoa1 > 30 and aoa1 < 90:
        phasedelay[4] = phasedelay[4]+(1*(2*sp.pi))    # [31:89]
    else:
        print "phase delay out-of-bounds"

    if aoa1 > -90 and aoa1 <= -54:
        phasedelay[5] = phasedelay[5]-(2*(2*sp.pi))    # [-89:-54]
    elif aoa1 > -54 and aoa1 <= -24:
        phasedelay[5] = phasedelay[5]-(1*(2*sp.pi))    # [-53:-24]
    elif aoa1 > -24 and aoa1 <= 23:
        phasedelay[5] = phasedelay[5]-(0*(2*sp.pi))    # [-23:23]
    elif aoa1 > 23 and aoa1 <= 53:
        phasedelay[5] = phasedelay[5]+(1*(2*sp.pi))    # [24:53]
    elif aoa1 > 53 and aoa1 < 90:
        phasedelay[5] = phasedelay[5]+(2*(2*sp.pi))    # [54:89]
    else:
        print "phase delay out-of-bounds"

    if aoa1 > -90 and aoa1 <= -42:
        phasedelay[6] = phasedelay[6]-(2*(2*sp.pi))    # [-89:-42]
    elif aoa1 > -42 and aoa1 <= -20:
        phasedelay[6] = phasedelay[6]-(1*(2*sp.pi))    # [-41:-20]
    elif aoa1 > -20 and aoa1 <= 19:
        phasedelay[6] = phasedelay[6]-(0*(2*sp.pi))    # [-19:19]
    elif aoa1 > 19 and aoa1 <= 41:
        phasedelay[6] = phasedelay[6]+(1*(2*sp.pi))    # [20:41]
    elif aoa1 > 41 and aoa1 < 90:
        phasedelay[6] = phasedelay[6]+(2*(2*sp.pi))    # [42:89]
    else:
        print "phase delay out-of-bounds"

    if aoa1 > -90 and aoa1 <= -58.9:
        phasedelay[7] = phasedelay[7]-(3*(2*sp.pi))    # [-89:-59]
    elif aoa1 > -59 and aoa1 <= -35:
        phasedelay[7] = phasedelay[7]-(2*(2*sp.pi))    # [-58:-35]
    elif aoa1 > -35 and aoa1 <= -17:
        phasedelay[7] = phasedelay[7]-(1*(2*sp.pi))    # [-34:-17]
    elif aoa1 > -17 and aoa1 <= 16:
        phasedelay[7] = phasedelay[7]-(0*(2*sp.pi))    # [-16:16]
    elif aoa1 > 16 and aoa1 <= 34:
        phasedelay[7] = phasedelay[7]+(1*(2*sp.pi))    # [17:34]
    elif aoa1 > 34 and aoa1 <= 58.9:
        phasedelay[7] = phasedelay[7]+(2*(2*sp.pi))    # [35:58]
    elif aoa1 > 58 and aoa1 < 90:
        phasedelay[7] = phasedelay[7]+(3*(2*sp.pi))    # [59:89]
    else:
        print("phase delay out-of-bounds")

    # Observations: every element has a non-ambigious range about 0 degrees. It shrinks as you get farther away from the reference element
    #               every 2 elements after reference have 2 additional ambigious ranges, one on either side of the 0 degree mark

    aoa = sp.arcsin((wavelength*phasedelay)/(2*sp.pi*array))
    aoa = aoa*180/sp.pi
    return aoa


if __name__ == "__main__":
    Fc = 1e9                          # carrier frequency
    wavelength = 3*pow(10,8) / Fc
    N=8
    array = buildAntennaArray(num_elements=N, ref_element=1, seperation=wavelength/2)

    plt.axis([-2,2,-2,2])
    plt.ion()
    plt.show()

    for angle in range(-89,89,1):
        plt.clf()
        plt.axis([-2,2,-2,2])
        plt.plot(array,sp.zeros(N, dtype = sp.double),'or')

        pd = calculatePhaseDelay(array, (angle*sp.pi/180), wavelength)
        aoa = calculateAOA(array, pd, wavelength)
        #print aoa

        for i in range(N):
            mag = 1
            ang = aoa[i]
            x = array[i]
            y = 0
            # remember, 0 degrees should be boresight
            plt.plot([x, mag*sp.cos((90-ang)*sp.pi/180)+x],[y, mag*sp.sin((90-ang)*sp.pi/180)+y], 'r')

        plt.draw()
        plt.pause(0.001)
        #raw_input("paused...")
