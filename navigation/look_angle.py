#! /usr/bin/env python

#------------------------------------------------------------------------------
# MIT License
#
# Copyright (C) 2021  Tim Cardenuto
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# For questions contact the originator at timcardenuto@gmail.com
#------------------------------------------------------------------------------

import sys
import csv
import argparse
import yaml
from math import pi, pow, sqrt, cos, sin, tan, atan2, asin


def convert_csv_to_yaml_list(filename):
    """
    Only works on copied table from https://www.n2yo.com/satellites/?c=10&srt=11&dir=1
    """
    dictdata = []
    with open(filename, "r") as stream:
        data = csv.reader(stream)
        for row in data:
            if row[4][-1:] == "W":
                longitude = float(row[4][:-1]) * -1
            elif row[4][-1:] == "E":
                longitude = float(row[4][:-1])
            else:
                longitude = float(row[4]) # Assume it's just a number
            r = {"name": str(row[0]), "norad_id": int(row[1]), "international_code": str(row[2]), "period": float(row[3]), "longitude": longitude}
            print(r)
            dictdata.append(r)
            
    yamldata = yaml.dump(dictdata, explicit_start=True, default_flow_style=False)
    f = open("geo-satellites.yaml", "w")
    f.write(yamldata) 
    f.close() 


def angle_to_gps_satellite(earth_lat, earth_lon):
    """
    Find look angle from point on earth to GPS satellite
    """

    # Earth station (radians)
    earth_lat = earth_lat * pi / 180.0
    earth_lon = earth_lon * pi / 180.0


    # Time of interest (s) for the GPS week of interest
    toi = 500000

    # Constants
    R = 6378137.00          # WGS84 mean Earth radius (m)
    u = 3.986004418e14      # WGS84 value of the earth’s gravitational constant for GPS user (m^3/s^2)
    sidereal_day = 23 + (56/60.0) + (4.1/3600)    # length of sidereal day in hours

    # Parameters given by GPS Almanac week 801, PRN-01
    a = pow(5153.614258,2)     # Almanac: Semi-major length (m), given by SQRT(A)^2
    e = 0.3765106201e-2             # Almanac: Eccentricity
    b = a * sqrt(1 - (e * e))     # Semi-minor length (m)
    rs = a * (1 - e)           # TODO: distance to satellite at perigree from earth center?
    M = 0.4480223834e0              # Almanac: Mean anomaly (rad)
    toa = 503808.0000               # Almanac: Time of applicability (s)
    w = 0.434909394                 # Almanac: Argument of perigree (rad)
    omega_week = 0.7017688714e0     # Almanac: Right ascension at week (rad)
    omega_rate = -0.7817468486e-8   # Almanac: Rate of right ascension (r/s)
    inc = 0.9617064849              # Almanac: Orbital inclination (rad)

    # Excel example
    # a = pow(5153.621582,2)
    # e = 9.066581726e-3
    # b = a * sqrt(1 - (e * e))
    # M = -0.1083742499e1
    # toa = 319488                  # Almanac: Time of applicability (s)
    # w = 2.875033498               # Almanac: Argument of perigree (rad)
    # omega_week = -0.1663788676e1  # Almanac: Right ascension at week (rad)
    # omega_rate = -0.8236384019e-8 # Almanac: Rate of right ascension (r/s)
    # inc = 0.9405422211        # Almanac: Orbital inclination (rad)
    # earth_lat = 0.680964
    # earth_lon = -1.346714
    # toi = 103338.990          # TODO Time of interest (s) in week?


    # 1) Estimate correct mean anomaly adjusted for time of interest
    T = (2 * pi / sqrt(u)) * (pow(sqrt(a),3)) # Period of satellite orbit
    Mm = M + (toi - toa) * 2 * pi / T

    # 2) Estimate eccentric anomaly E using Newtons method
    E = 0                         # initial guess
    for i in range(10):           # iterate
        M = E - e * sin(E) - Mm   # f(E)
        Md = 1 - e * cos(E)       # derivative of f(E)
        E = E - M / Md            # update step

    print("E estimate: "+str(E))

    # 3) Calculate position of satellite in orbital plane
    x0 = a*(cos(E-e)) # NOTE I originally had a.*(cos(E)-e) , not sure what's right
    y0 = b*sin(E)
    r0 = sqrt(x0*x0 + y0*y0) # radius in orbital plane, not from earth station

    # Estimate correct right ascension adjusted for time of interest
    omega = omega_week + omega_rate*toi

    # 4) Transform orbital to earth centered coord
    Px = cos(w)*cos(omega) - sin(w)*sin(omega)*cos(inc)
    Py = cos(w)*sin(omega) + sin(w)*cos(omega)*cos(inc)
    Pz = sin(w)*sin(inc)
    Qx = -sin(w)*cos(omega) - cos(w)*sin(omega)*cos(inc)
    Qy = -sin(w)*sin(omega) + cos(w)*cos(omega)*cos(inc)
    Qz = cos(w)*sin(inc)

    # Satellite position in geocentric coord
    x = Px*x0 + Qx*y0
    y = Py*x0 + Qy*y0
    z = Pz*x0 + Qz*y0
    r = sqrt(x*x + y*y + z*z)

    asc = atan2(y, x) # right ascension
    # Not sure the offset is calculated correctly... using mod to fix fractional day and longitude values
    asc_off = (((((toi/3600)%sidereal_day)/sidereal_day)*360)%180) * pi/180 # right ascension longitude offset (rad)
    sat_lon = asc - asc_off # convert satellite right ascension into longitude
    dec = atan2(z, sqrt(x*x + y*y)) # declination
    ns = asin(sin(dec)*sin(earth_lat) + cos(dec)*cos(earth_lat)*cos((sat_lon-earth_lon))) # NOTE excel swaps sat_lon-earth_lon
    el = atan2((sin(ns)-R/r), cos(ns))
    az = atan2(sin((sat_lon-earth_lon)), (cos(earth_lat)*tan(dec) - sin(earth_lat)*cos((sat_lon-earth_lon))))

    # convert from radians to degrees
    el = el * 180/pi
    az = az * 180/pi
    print("Elevation (deg): "+str(el))
    print("Azimuth (deg):   "+str(az))

    sat_lat = asin(z/r) * 180/pi
    sat_lon = sat_lon * 180/pi
    #TODO gives different results than the above method
    # sat_lon = atan2(y,x)*180/pi 

    print("Sat Lat (deg):   "+str(sat_lat))
    print("Sat Lon (deg):   "+str(sat_lon))



def angle_to_geo_satellite(earth_lat, earth_lon, sat_lat, sat_lon):
    """ 
    Find look angle from point on earth to geostationary satellite
    """
    # Constants
    R = 6378137.00      # WGS84 mean Earth radius (m)
    f = 1/298.257223563 # WGS84 flattening value
    e = 2*f - f*f

    # Earth Station
    elat = earth_lat*pi/180
    elon = earth_lon*pi/180

    # Geostationary Satellite
    slat = sat_lat*pi/180  # geostationary has zero latitude (above equator)
    slon = sat_lon*pi/180
    salt = 35787000         # geostationary satellites have an average altitude of 35,786 km above sea level
    N = R/sqrt(1 - pow((e*sin(slat)),2))
    xs = (N+salt)*cos(slat)*cos(slon)
    ys = (N+salt)*cos(slat)*sin(slon)
    zs = (N*(1-e*e)+salt)*sin(slat)
    rs = sqrt(xs*xs + ys*ys + zs*zs) # also = R+sat_alt

    asc = atan2(ys, xs)                     # right ascension
    inc = atan2(zs, sqrt(xs*xs + ys*ys))    # inclination (called declination when natural satellite)
    ns = asin(sin(inc)*sin(elat) + cos(inc)*cos(elat)*cos((slon-elon)))
    el = atan2((sin(ns)-R/rs), cos(ns))
    az = atan2(sin((slon-elon)), (cos(elat)*tan(inc) - sin(elat)*cos((slon-elon))))
    if az < 0:
        az = az + 2*pi # convert back to a positive reference from true north
    r = sqrt(rs*rs - R*R*cos(el)*cos(el)) - R*sin(el)

    print("Calculated")
    print("  Sat Azimuth (deg):   "+str(az*180/pi))
    print("  Sat Elevation (deg): "+str(el*180/pi))  # Ho in nautical almanac 
    print("  Sat Range (m):       "+str(r))




def print_list_as_columns(alist):
    # sort them alphabetically
    alist.sort()
    
    # print them in 4 columns (should fit on normal half screen width terminal) 
    count = 0
    for a,b,c,d in zip(alist[::4],alist[1::4],alist[2::4],alist[3::4]):
        count = count + 4
        print('{:<25}{:<25}{:<25}{:<}'.format(a,b,c,d))

    # TODO the above fancy thing won't print a line if it's less than 4 items
    leftover = len(alist) - count
    if leftover == 3:
        print('{:<25}{:<25}{:<}'.format(alist[-3:-2][0], alist[-2:-1][0], alist[-1:][0]))
    elif leftover == 2:
        print('{:<25}{:<}'.format(alist[-2:-1][0], alist[-1:][0]))
    elif leftover == 1:
        print('{:<}'.format(alist[-1:][0]))
    else:
        print("How did this happen??")
        sys.exit(1)
        
        
        
def list_known_satellites(afilter=""):
    print("Known satellites:")
    print("----------------------------------------------------------------------------------------------")
    with open("geo-satellites.yaml", "r") as stream:
        try:
            data = yaml.safe_load(stream)
        except yaml.YAMLError as ex:
            print(ex)
        
    # get list of satellite names
    names = []    
    for sat in data:
        if afilter in sat['name']:
            names.append(sat['name'])
            
    print_list_as_columns(names)
    print("----------------------------------------------------------------------------------------------")
    sys.exit(0)
    
    

def get_satellite_longitude(sat_name):
    with open("geo-satellites.yaml", "r") as stream:
        try:
            data = yaml.safe_load(stream)
        except yaml.YAMLError as ex:
            print(ex)
            
    found = []
    for sat in data:
        if sat_name in sat['name']:
            sat_lon = float(sat['longitude'])
            found.append(sat['name'])

    if len(found) == 0:
        print("[ERROR] Could not find entry for satellite `"+sat_name+"`")
        sys.exit(1)
    elif len(found) > 1:
        print("[ERROR] Found duplicate entries for satellite `"+sat_name+"`  Perhaps you meant one of these?")
        print_list_as_columns(found)
        sys.exit(1)
    else:
        return sat_lon
    


def PLANGA(IRLE, ILNK, XE, YE, ZE, IPS, PAS, XS, YS, ZS, XA, ya, za, ipse, pase, xse, yse, zse, xae, yae, zae, parw, dpaa):
    """
    Args
        IRLE - index for the reference line, relative to which the polarization angle of the radio wave is measured at the earth point
             = 1 for a line parallel to the equatorial plane
             = 2 for a line parallel to the local horizontal plane (default)
        ILNK - index for uplink or downlink case
             = 1 for uplink
             = 2 for downlink
        XE, YE, ZE - earth center coordinates of the earth point
        IPS - index of the polarization type of the antenna of the satellite in question
            = 0 for circular
            = 1 for linear when the angle is specified relative to a line parallel to the equatorial plane at the aim point
            = 2 for linear when the angle is specified relative to a line parallel to the local horizontal plane at the aim point
        PAS - polarization angle (degrees) of the antenna of the satellite in question specified at its aim point (must be +45 or -45 when IPS = 0)
        XS, YS, ZS - earth center coordinates of the satellite in question
        XA, YA, ZA - earth center coordinate of the aim point (on surface of the earth) of the atenna of the satellite in question
        IPSE - index of the polarization type of the antenna of the satellite that servers the earth point
             = 0 for circular
             = 1 for linear when the angle is specified relative to a line parallel to the equatorial plane at the aim point
             = 2 for linear when the angle is specified relative to a line parallel to the local horizontal plane at the aim point
         PASE - polarization angle (degrees) of the antenna of the satellite that serves the earth point specified at its aim point (must be +45 or -45 when IPSE = 0
         XSE, YSE, ZSE - earth center coordinates of the satellite that serves the earth point
         XAE, YAE, 
         
         ...
    """

	
if __name__ == "__main__":   
    parser = argparse.ArgumentParser(description='Calculate look angle from Earth to Satellite')
    parser.add_argument('--elat', type=float, default=0, help='Earth station latitude')
    parser.add_argument('--elon', type=float, default=0, help='Earth station longitude')
    parser.add_argument('--slat', type=float, default=0, help='Satellite latitude')
    parser.add_argument('--slon', type=float, default=0, help='Satellite longitude')
    parser.add_argument('--sat', type=str, help='Satellite name.')
    parser.add_argument('--list', nargs="*", help='List known satellites. Can use a substring as an argument to filter for.')
    args = parser.parse_args()

    # NOTE: satellite list is based on what's listed at https://www.n2yo.com/satellites/?c=10&srt=11&dir=1
    if args.list is not None:
        try:
            filter_for = args.list[0]
        except:
            filter_for = ""
        list_known_satellites(filter_for)

    if args.sat is not None:
        args.slon = get_satellite_longitude(args.sat)
    
    
    print("Provided")
    print("  Earth Station lat/lon: "+str(args.elat)+","+str(args.elon))   
    print("  Satellite lat/lon:     "+str(args.slat)+","+str(args.slon))
    
    angle_to_geo_satellite(args.elat, args.elon, args.slat, args.slon)
        
    
