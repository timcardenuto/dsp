import numpy as np



def estimateGeolocationTDOA(mtype, z, sigma, sref, srel):
    """
        All measurements should use the same reference sensor A, I don't think you can combine other ones
        e.g. if you have 3 sensors A, B, C calculate the TDOA between A and B (tdoaAB) and between A and C (tdoaAC) 
        and include those in the measurements structure. Don't include TDOA between B and C (tdoaBC). You should 
        end up with N - 1 measurements where N is the number of sensors.

        The tags are kept common for parsing, "sensorRef" is the reference sensor location (sensorA in the above example)
        "sensorRel" is the second sensor that the measurement is related to (sensorB or sensorC in the above example)
        "tdoa" is the TDOA measurement between sensorRef and sensorRel
        
       measurements = [
            {'type': 'tdoa_range', 'value': tdoa_range_value, 'sigma': value, 'sensorRef': [x,y], 'sensorRel': [x,y]},
            {'type': 'tdoa_range', 'value': tdoa_range_value, 'sigma': value, 'sensorRef': [x,y], 'sensorRel': [x,y]},
            ...
       ] 
    """
    print("############################")
    print(">> estimateGeolocationTDOA()")
    print("############################")
    print("Measurement Types:        "+str(mtype))
    print("Measurement Values (z):   "+str(z))
    print("Measurement Stds (sigma): "+str(sigma))

    a = np.unique(sref, axis=0)
    b = np.unique(srel, axis=0)
    uniquelocs = np.concatenate((a, b))
    print("Sensor Locations [x,y]:   "+str(uniquelocs))
    # R = np.identity(len(sigma))*np.vstack(np.hstack(sigma)*np.hstack(sigma))       # was a more complicated version of the next line, not sure why I was using this in some places....
    R = (sigma*sigma)*np.identity(len(sigma))    # covariance matrix (measurement errors)
    print("Covariance Matrix (R):    "+str(R))

    # Initial estimate: centroid of all unique sensor locations.
    # The Moore-Penrose formula requires DOA angles, but z here are signed TDOA range
    # differences (meters), so that approach produces nonsense for TDOA. The centroid
    # is a neutral starting point for the ILS to converge from.
    xhat = np.mean(uniquelocs, axis=0)
    print("xhat: "+str(xhat))
    print("")

    count = 1
    maxcount = 10
    error = [100,100]
    maxerror = .1
    while (count < maxcount and (error[0]+error[1]/2) > maxerror):
        h = np.zeros(len(sref))
        H = np.zeros((len(sref),2))
        # TODO: Probably a way to do this w/o loop using matrices but too lazy to figure it out...
        for i in range(len(sref)):
            h[i] = (np.power(np.transpose(xhat-sref[i])@(xhat-sref[i]),0.5) - np.power(np.transpose(xhat-srel[i])@(xhat-srel[i]),0.5))
            H[i] = ((1/np.power(np.transpose(xhat-sref[i])@(xhat-sref[i]),0.5))*np.transpose(xhat-sref[i]) - (1/np.power(np.transpose(xhat-srel[i])@(xhat-srel[i]),0.5))*np.transpose(xhat-srel[i]))

        print("h: "+str(h))
        print("H: "+str(H))
        print("")

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

    return xhat,P


def _intersect_doa_tdoa(doa_loc, theta, sref_tdoa, srel_tdoa, z_tdoa, t_max):
    """Find the point along the DOA ray that lies on a TDOA hyperbola.

    Parametrises the DOA ray as P(t) = doa_loc + t*[cos θ, sin θ] and solves
    the TDOA equation dist(P, sref_tdoa) - dist(P, srel_tdoa) = z_tdoa for t.

    This turns two 2-D constraints into one 1-D root-finding problem, giving an
    exact (noise-free) initial estimate for the ILS when there is exactly one DOA
    and at least one TDOA measurement.

    Uses a scan-then-bisect approach (no scipy required).

    Returns the [x, y] intersection point, or None if no root found in [1, t_max].
    """
    cx, cy   = float(doa_loc[0]), float(doa_loc[1])
    cos_t    = np.cos(theta)
    sin_t    = np.sin(theta)
    sx1, sy1 = float(sref_tdoa[0]), float(sref_tdoa[1])
    sx2, sy2 = float(srel_tdoa[0]), float(srel_tdoa[1])

    def residual(t):
        x  = cx + t * cos_t
        y  = cy + t * sin_t
        d1 = np.sqrt((x - sx1)**2 + (y - sy1)**2)
        d2 = np.sqrt((x - sx2)**2 + (y - sy2)**2)
        return (d1 - d2) - z_tdoa

    # Scan the positive ray direction for a sign change, then bisect
    t_vals = np.linspace(1.0, t_max, 600)
    f_prev = residual(t_vals[0])
    for k in range(1, len(t_vals)):
        f_curr = residual(t_vals[k])
        if f_prev * f_curr < 0:
            a, b = t_vals[k - 1], t_vals[k]
            for _ in range(60):            # 60 bisection steps → ~1e-18 relative error
                mid = (a + b) / 2
                if residual(mid) * residual(a) < 0:
                    b = mid
                else:
                    a = mid
            t_root = (a + b) / 2
            return np.array([cx + t_root * cos_t, cy + t_root * sin_t])
        f_prev = f_curr
    return None


def estimateGeolocationMulti(mtype, z, sigma, sref, srel):
    """ Function to allow geolocation when multiple different measurement types are combined

    """
    print("#############################")
    print(">> estimateGeolocationMulti()")
    print("#############################")
    print("Measurement Types:        "+str(mtype))
    print("Measurement Values (z):   "+str(z))
    print("Measurement Stds (sigma): "+str(sigma))

    a = np.unique(sref, axis=0)
    # srel has [None, None] placeholder entries for DOA measurements; filter those out
    srel_valid = np.array([s for s in srel if s[0] is not None], dtype=float)
    uniquelocs = np.unique(np.concatenate((a, srel_valid)), axis=0) if len(srel_valid) else a
    print("Sensor Locations [x,y]:   "+str(uniquelocs))

    R = (sigma*sigma)*np.identity(len(sigma))    # covariance matrix (measurement errors)
    print("Covariance Matrix (R):    "+str(R))

    doa_idxs  = [i for i, mt in enumerate(mtype) if mt == 'doa_angle']
    tdoa_idxs = [i for i, mt in enumerate(mtype) if mt == 'tdoa_range']

    if len(doa_idxs) >= 2:
        # Moore-Penrose on ≥2 DOA measurements gives a good triangulated initial guess.
        # (Requires ≥2 rays; a single ray gives a degenerate projection, not a point.)
        Msum = []
        MTasum = []
        for i in doa_idxs:
            a_loc = sref[i]
            t = z[i]
            u = np.vstack([[np.cos(t)], [np.sin(t)]])
            M = np.identity(2) - u@(u.conj().transpose())
            Msum = M if len(Msum) == 0 else Msum + M
            MTasum = M@a_loc if len(MTasum) == 0 else MTasum + M@a_loc
        xhat = np.linalg.pinv(Msum)@MTasum

    elif len(doa_idxs) == 1 and len(tdoa_idxs) >= 1:
        # Exactly 1 DOA + ≥1 TDOA: substitute the ray P(t) = doa_loc + t*[cos θ, sin θ]
        # into the TDOA equation to get a 1-D root-finding problem.  This gives the
        # exact intersection of the DOA ray with the TDOA hyperbola, which is an ideal
        # starting point for the ILS.
        max_sep = max(
            (np.linalg.norm(uniquelocs[i] - uniquelocs[j])
             for i in range(len(uniquelocs)) for j in range(i + 1, len(uniquelocs))),
            default=1e4
        )
        t_max = max(max_sep * 50, 1e5)   # generous range: 50× sensor baseline or 100 km

        xhat = None
        for tdoa_idx in tdoa_idxs:
            candidate = _intersect_doa_tdoa(
                sref[doa_idxs[0]], z[doa_idxs[0]],
                sref[tdoa_idx], srel[tdoa_idx],
                z[tdoa_idx], t_max
            )
            if candidate is not None:
                xhat = candidate
                break
        if xhat is None:
            print("  WARNING: DOA/TDOA ray-hyperbola intersection not found, falling back to centroid")
            xhat = np.mean(uniquelocs, axis=0)

    else:
        # No DOA measurements (pure TDOA) or no TDOA to pair with: use sensor centroid.
        xhat = np.mean(uniquelocs, axis=0)

    print("xhat: "+str(xhat))
    print("")

    # DOA h and H
    # h(x) = theta = arctan((x2-a2),(x1-a1))
    # H = (1/r)*[-sin(h) cos(h)]
    # Range h and H
    # h = r = ((x-a)'*(x-a))^.5
    # H = (1/r)*(x-a)'
    # TODO: Check that changes to how the old DOA algorithm worked are still working, I changed syntax to simplify and match the TDOA version (no vstack/hstack)
    # TODO: Check that transpose type is consistent and correct, not sure if I have or should have complex numbers...
    # TODO: Check that pinv vs. solve works or which is better
    count = 1
    maxcount = 10               # limit number of runs w/ a maximum count
    error = np.hstack([100,100])
    maxerror = 0.1              # when to stop iterating (acceptable error)
    while (count < maxcount and max(error) > maxerror):
        h = np.zeros(len(sref))
        H = np.zeros((len(sref),2))
        for i in range(len(sref)):
            if mtype[i] == 'tdoa_range':
                h[i] = (np.power(np.transpose(xhat-sref[i])@(xhat-sref[i]),0.5) - np.power(np.transpose(xhat-srel[i])@(xhat-srel[i]),0.5))
                H[i] = ((1/np.power(np.transpose(xhat-sref[i])@(xhat-sref[i]),0.5))*np.transpose(xhat-sref[i]) - (1/np.power(np.transpose(xhat-srel[i])@(xhat-srel[i]),0.5))*np.transpose(xhat-srel[i]))
            elif mtype[i] == 'doa_angle':
                xtemp = xhat - sref[i]
                h[i] = (np.arctan2(xtemp[1], xtemp[0]))
                H[i] = ((1/np.power(np.transpose(xhat-sref[i])@(xhat-sref[i]),0.5))*np.transpose([-np.sin(z[i]), np.cos(z[i])]))
            else:
                print("ERROR unknown type? "+str(mtype[i]))

        # NOTE: np.transpose(H) performs a standard matrix transpose, which only swaps rows and columns.
        # P = np.linalg.inv(np.transpose(H).dot(np.linalg.inv(R)).dot(H))             # TDOA version
        # NOTE: H.conj().transpose() (or the equivalent H.conj().T) performs a conjugate transpose (also called a Hermitian transpose), which swaps rows and columns and replaces each entry with its complex conjugate.
        # P = np.linalg.inv(H.conj().transpose()@np.linalg.inv(R)@H)                  # DOA version
        # NOTE: If your matrices are large or nearly singular, np.linalg.inv can be numerically unstable. Consider using np.linalg.solve or np.linalg.pinv (pseudo-inverse) for better accuracy.
        P = np.linalg.pinv(H.conj().transpose()@np.linalg.pinv(R)@H)

        # xhatnew = xhat + P.dot(np.transpose(H)).dot(np.linalg.inv(R)).dot(z-h)      # TDOA version
        # xhatnew = xhat + P@H.conj().transpose()@np.linalg.inv(R)@(theta-thetahat)   # DOA version
        xhatnew = xhat + P@H.conj().transpose()@np.linalg.pinv(R)@(z-h)

        error = abs(xhatnew - xhat)
        xhat = xhatnew
        print("count = ",count) 
        print("Covariance = ",P)
        print("xhat = ",xhat)
        print("error = ",error)
        print("")
        count = count+1

    return xhat,P



def estimateGeolocationDOA(mtype, z, sigma, s1loc):
    print("###########################")
    print(">> estimateGeolocationDOA()")
    print("###########################")
    print("Measurement Types:        "+str(mtype))
    print("Measurement Values (z):   "+str(z))
    print("Measurement Stds (sigma): "+str(sigma))

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
    for a,t in zip(s1loc,z):
        u = np.vstack([[np.cos(t)], [np.sin(t)]])           # unit vector along DOA
        M = np.identity(2) - u@(u.conj().transpose())       # Moore-Penrose matrix
        Msum = M if len(Msum) == 0 else Msum + M
        MTasum = M@a if len(MTasum) == 0 else MTasum + M@a
    xhat = np.linalg.pinv(Msum)@MTasum
    print("xhat: "+str(xhat))
    print("")

    # plt.plot(xhat[0], xhat[1], '*r', label='Initial Target Location Estimate w/ Moore-Penrose')

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

    # theta = np.vstack([[z[0]]])
    # for t in z[1:]:
        # theta = np.vstack([theta, [t]])

    # measurements (z)
    theta = np.vstack(z)

    # covariance matrix (measurement errors)
    R = np.identity(len(sigma))*np.vstack(np.hstack(sigma)*np.hstack(sigma))
    print("Covariance Matrix (R): "+str(R))

    count = 1
    maxcount = 10               # limit number of runs w/ a maximum count
    error = np.hstack([100,100])
    maxerror = 0.1              # when to stop iterating (acceptable error)
    while (count < maxcount and max(error) > maxerror):
        # print("Iteration: ",count)
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

        xtemp = xhat - s1loc[0]
        thetahat = np.arctan2(xtemp[1], xtemp[0])
        for a in s1loc[1:]:
            xtemp = xhat - a
            thetahat = np.vstack([thetahat, [np.arctan2(xtemp[1], xtemp[0])]])

        # H = np.vstack([
        #     (1/np.power((xhat-loc1).conj().transpose()@(xhat-loc1),0.5))*np.hstack([-np.sin(thetahat[0]), np.cos(thetahat[0])]),
        #     (1/np.power((xhat-loc2).conj().transpose()@(xhat-loc2),0.5))*np.hstack([-np.sin(thetahat[1]), np.cos(thetahat[1])]),
        #     (1/np.power((xhat-loc3).conj().transpose()@(xhat-loc3),0.5))*np.hstack([-np.sin(thetahat[2]), np.cos(thetahat[2])])
        #     ])

        H = (1/np.power((xhat-s1loc[0]).conj().transpose()@(xhat-s1loc[0]),0.5))*np.hstack([-np.sin(thetahat[0]), np.cos(thetahat[0])])
        for a,th in zip(s1loc[1:],thetahat[1:]):
            H = np.vstack([H, (1/np.power((xhat-a).conj().transpose()@(xhat-a),0.5))*np.hstack([-np.sin(th), np.cos(th)])])

        print("h: "+str(thetahat))
        print("H: "+str(H))
        print("")

        P = np.linalg.pinv(H.conj().transpose()@np.linalg.pinv(R)@H)

        xhatnew = xhat + P@H.conj().transpose()@np.linalg.pinv(R)@(theta-thetahat)
        error = abs(xhatnew - xhat)
        xhat = xhatnew
        print("count = ",count) 
        print("Covariance = ",P)
        print("xhat = ",xhat)
        print("error = ",error)
        print("")
        count = count+1

        # axs.plot(xhat[0], xhat[1], '*k', label="Estimate Update, Error={}".format(error))
        # axs.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
        # axs.set(xlim=(-5, 100), ylim=(-5, 100))
        # axs.set_aspect('equal')
        # plt.draw()
        # input("Press Enter to continue...")

    return xhat,P



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

    # print("Eigen Values = ",eigenvalues)
    # print("Eigen Vectors = ",eigenvectors)
    # print("Semi-major Axis = ",semimajor)
    # print("Semi-minor Axis = ",semiminor)
    # print("")

    return semimajor,semiminor,angle,ellipse


# def geolocate(loc_array, theta_array, sigma_array, containment):

def geolocate(measurements, containment):
    """
        measurements is a list of dicts that look like this:

        measurements = [
            {'mtype': 'tdoa_range', 'value': tdoa_range_value, 'sigma': value, 's1loc': [x,y], 's2loc': [x,y]},
            {'mtype': 'doa_angle', 'value': doa_angle_value, 'sigma': value, 's1loc': [x,y]},
            ...
       ] 

        mtype = type of measurement, expected values are "doa_angle" or "tdoa_range"
        s1loc = primary sensor location [x,y] in meters, this is the sensor the measurement is relative to. 
        value = sensor measurement (float). They can include different measurement types: DOA angle (radians), TDOA range (meters)
        sigma = sensor measurement standard deviations (float). For DOA it's expected to be in radians. For TDOA ranges it's expected to be in meters.
        s2loc = secondary sensor location [x,y] in meters. This is for TDOA measurements only, it's the second sensor location that the time difference was calulated relative to the primary sensor location.
        containment = float between 0 and 100 exclusive, scales the ellipse to ensure this value is the probability that the target is within the ellipse 
    """
    mtype = []        # Measurement types
    z = []            # All measurements (DOA, TDOA range, etc)
    sigma = []        # Measurement errors (std deviation)
    s1loc = []        # Reference sensor1 location for TDOA or DOA measurement
    s2loc = []        # Reference sensor2 location for TDOA measurement (only relevant to TDOA measurement types)

    # unpack the incoming measurements
    for measurement in measurements:
        mtype.append(measurement['mtype'])
        z.append(measurement['value'])
        sigma.append(measurement['sigma'])
        s1loc.append(measurement['s1loc'])
        if measurement['mtype'] == 'tdoa_range':
            s2loc.append(measurement['s2loc'])
        else:
            s2loc.append([None, None])   # There is no second sensor for DOA measurement type, but we need the indexing to match for Multi geo request

    # z = np.array(z)
    sigma = np.array(sigma)
    s1loc = np.array(s1loc)
    s2loc = np.array(s2loc)

    print("  |--mtype:        "+str(mtype))
    print("  |--measurement:  "+str(z))
    print("  |--sigma:        "+str(sigma))
    print("  |--sensor1:      "+str(s1loc))
    print("  |--sensor2:      "+str(s2loc))
    print("  |--containment:  "+str(containment))

    # Use Moore-Penrose and Iterated Least Squares (ILS) to estimate the geolocation and covariance matrix
    if 'tdoa_range' not in mtype:
        xhat,P = estimateGeolocationDOA(mtype, z, sigma, s1loc)
        # TODO: Why is this function output an extra [[]] around the measurement?
        xhat = [xhat[0][0], xhat[1][0]]
    elif 'doa_angle' not in mtype:
        xhat,P = estimateGeolocationTDOA(mtype, z, sigma, s1loc, s2loc)
    else:
        xhat,P = estimateGeolocationMulti(mtype, z, sigma, s1loc, s2loc)

    print("  |--xhat:         "+str(xhat))
    print("  |--P:            "+str(P))

    # Calculate the Elliptical Error Probable (EEP), i.e. error ellipse
    semimajor,semiminor,orientation,shape = calculateEllipse(xhat, P, containment)

    print("  |--semimajor:    "+str(semimajor))
    print("  |--semiminor:    "+str(semiminor))
    print("  |--orientation:  "+str(orientation))

    # Plot estimated target location and EEP
    # axs.plot(xhat[0], xhat[1], '*g', label='Final Target Location Estimate')
    # axs.plot(ellipse[:,0],ellipse[:,1],'-g')
    # axs.legend(bbox_to_anchor=(1, 0.5), loc='center left', shadow=True, fontsize='small')
    # axs.set(xlim=(-5, 100), ylim=(-5, 100))
    # axs.set_aspect('equal')

    ellipse = {'x': xhat[0], 'y': xhat[1], 'semimajor': semimajor, 'semiminor': semiminor, 'orientation': orientation, 'containment': containment, 'shape': shape}

    return ellipse
