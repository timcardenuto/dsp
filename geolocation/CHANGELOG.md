# Changelog



### Added
❯ now that we have a geographic coord system are there any improvements that can be made to make it     
more realistic                                                                                          
                                                                                                        
● Good question. Let me read the measurement simulation code carefully to identify the gaps before      
  proposing anything.
                                                                                                        
  Read 1 file (ctrl+o to expand)                                                                        
                                                                                                        
● Now I have a clear picture of all the gaps. Here's what the current simulation does that's physically 
  incorrect when using the geo estimator:                                                               
                                                                                                        
  1. DOA bearing: atan2(Δy, Δx) — flat-earth math angle. Vincenty gives the true geodesic bearing, which
   diverges noticeably beyond ~10 km especially at oblique angles.                                      
  2. TDOA distances: Euclidean sqrt(Δx²+Δy²) — should be geodesic. Flat-earth diverges at long range.   
  3. Speed of propagation: vacuum c = 299,792,458 m/s. In the troposphere it's c/n ≈ 299,705,000 m/s    
  (n≈1.000293). Matters for the sigma-time → sigma-range conversion.                                    
  4. Convention mismatch: the existing calculate_measurements produces math-angle DOA, and              
  calculate_geolocation_geo converts it. If we fix the simulation to use Vincenty, we should keep math  
  convention in storage so that conversion path stays consistent.                                       
                                                                                                        
  The fix is a calculate_measurements_geo that uses _vincenty_inverse (already in geolocation.py,       
  imported via *) for both DOA bearing and TDOA distances, uses c_air, and wires into btn-measure when  
  in geo mode.                                                                                          
                  

### Added
❯ Can you add a new estimateGeolocationMulti() function/algorithm that uses a geographic coordinate     
system (lat/lon) instead of cartesian, using a WGS84 model and True North as the reference for DOA's?   
                                      
● Now I have a thorough understanding of the existing algorithm. The existing estimateGeolocationMulti
  uses:                                                  
  - Cartesian (x, y) in meters
  - DOA in math angle convention (0° = East, counterclockwise, from atan2)                              
  - Flat-earth geometry                                                   
                                                                                                        
  I need to add a new version that uses:                                                              
  - WGS84 (lat/lon in degrees)                                                                          
  - DOA as True North azimuth (0° = North, clockwise)                                                   
  - Vincenty geodesics for distance and bearing                                                         
  - Jacobians derived from WGS84 radii of curvature                                                     
                                                                                                        
  The ILS structure stays the same — only the measurement model h(x) and Jacobian H change. 

  What was added (geolocation.py)                                                                       
                                                                                                        
  Module-level WGS84 constants                                                                          
                                                                                                        
  _WGS84_A, _WGS84_F, _WGS84_B, _WGS84_E2 — shared by all geo helpers.                                  
   
  _wgs84_radii(lat_rad) → (M, N)                                                                        
                  
  Returns the WGS84 meridian (M) and prime-vertical (N) radii of curvature. These are the metric        
  conversion factors: M = metres per radian of latitude, N·cos(lat) = metres per radian of longitude.
                                                                                                        
  _vincenty_inverse(lat1, lon1, lat2, lon2) → (dist_m, az_rad)                                          
   
  Vincenty's inverse formula on WGS84. Returns geodesic distance and forward azimuth (True North,       
  clockwise). Handles the antipodal edge case via the sinSig < 1e-15 guard.
                                                                                                        
  estimateGeolocationMultiGeo(mtype, z, sigma, sref_deg, srel_deg) → (xhat, P)                          
   
  Core ILS algorithm. Key differences from the Cartesian version:                                       
                  
  ┌────────────────────┬───────────────────────────────────┬───────────────────────────────────────┐    
  │                    │            Cartesian              │              Geographic               │  
  │                    │    (estimateGeolocationMulti)     │     (estimateGeolocationMultiGeo)     │    
  ├────────────────────┼───────────────────────────────────┼───────────────────────────────────────┤  
  │ State              │ [x_m, y_m]                        │ [lat_rad, lon_rad]                    │  
  ├────────────────────┼───────────────────────────────────┼───────────────────────────────────────┤    
  │ DOA reference      │ East (math atan2)                 │ True North (clockwise)                │    
  ├────────────────────┼───────────────────────────────────┼───────────────────────────────────────┤    
  │ Distance           │ Euclidean                         │ Vincenty geodesic                     │    
  ├────────────────────┼───────────────────────────────────┼───────────────────────────────────────┤  
  │ Jacobian ∂h/∂x     │ Flat-earth derivatives            │ WGS84 radii of curvature              │    
  ├────────────────────┼───────────────────────────────────┼───────────────────────────────────────┤
  │ DOA residual       │ No                                │ Yes (% 2π, -π to π)                   │    
  │ wrapping           │                                   │                                       │
  ├────────────────────┼───────────────────────────────────┼───────────────────────────────────────┤    
  │ Convergence check  │ abs(Δx) < 0.1 m                   │ abs(Δ[lat,lon]) × metric < 0.01 m     │
  └────────────────────┴───────────────────────────────────┴───────────────────────────────────────┘    
   
  Initial estimate uses the same Moore-Penrose / DOA+TDOA ray-bisect / centroid hierarchy, implemented  
  in a flat-earth NE frame at the sensor centroid.
                                                                                                        
  geolocate_geo(measurements, containment) → dict                                                       
   
  Wrapper matching the geolocate() interface. After convergence it converts P from (rad_lat, rad_lon)²  
  to (North_m, East_m)² via J = diag(M, N·cos(lat)), then calls calculateEllipse. Returns lat, lon
  (degrees), semimajor/semiminor (metres), orientation (CCW from North, radians), and shape as (N, 2)   
  North/East metre offsets from centre.


### Added
Too slow now, redraws for sensor/target locations and measurements is several seconds slow.

● The culprit is gs_snap. At most zoom levels it hits the 250 cap, giving 250×250 = 62,500 invisible    
  points vs the original 100×100 = 10,000 — a 6× blowup in JSON payload size on every redraw (which
  fires on every store change and every pan/zoom).                                                      
   
  The resolution improvement actually comes entirely from centering the grid on the viewport, not from  
  having more points. With the same GS=100, at zoom 18 the spacing drops from ~2km to ~4.6m
  automatically because the extent shrinks. The max-size formula is unnecessary. 

● Same 10,000 invisible snap points as before, but now centered on the viewport and scaled to the       
  visible extent. Resolution by zoom:                                                                   
                                                                                                        
  ┌──────┬─────────────┬──────────────┐                                                                 
  │ Zoom │ Half-extent │ Snap spacing │                                                                 
  ├──────┼─────────────┼──────────────┤                                                               
  │ 10   │ ~58 km      │ ~1.2 km      │                                                                 
  ├──────┼─────────────┼──────────────┤
  │ 14   │ ~3.7 km     │ ~74 m        │                                                                 
  ├──────┼─────────────┼──────────────┤                                                                 
  │ 18   │ ~228 m      │ ~4.6 m       │                                                                 
  ├──────┼─────────────┼──────────────┤                                                                 
  │ 19   │ ~114 m      │ ~2.3 m       │
  ├──────┼─────────────┼──────────────┤                                                                 
  │ 20   │ ~57 m       │ ~1.1 m       │
  └──────┴─────────────┴──────────────┘                                                                 
                  
  The 1m-class resolution kicks in around zoom 20 (street/building level) without any performance cost  
  over the original.

### Added

  Fix 1 — DOA fill (lines 139-147)                                                                      
                                                                                                        
  The old code built a 3-point triangle [sensor, upper_end, lower_end] for the fill polygon. Plotly's   
  Scattermap (Mapbox-backed) renders fill='toself' based on the winding order of the polygon, which     
  alternates between CW and CCW depending on the DOA angle. This caused the fill to appear or disappear
  inconsistently.                                                                                       
                                                                                                      
  The fix builds a proper wedge/pie-slice polygon: sensor → arc of points sweeping from upper to lower  
  bound → back to sensor. The arc gives the polygon a consistent non-degenerate area regardless of
  angle, and Mapbox fills it reliably.                                                                  
                                                                                                      
  ---
  Fix 2 — Dynamic snap grid (lines 67-86)
                                         
  The old grid always spanned ±100km from origin at fixed GS=100, giving ~2km snap resolution. When
  zoomed in, there might be zero or one snap point visible.                                             
  
  The new grid:                                                                                         
  - Centers on the current viewport (center_lat/lon from view_range)                                  
  - Scales extent to zoom level: half_extent ≈ 40Mm / 2^zoom × 1.5  — at zoom 10 this is ~58km, at zoom 
  18 it's ~228m, at zoom 19 it's ~114m                                                                 
  - Dynamically picks gs_snap targeting 1m spacing, capped at 250 points per axis (62,500 total)        
                                                                                                      
  Resolution by zoom level:                                                                             
                                                                                                      
  ┌──────┬─────────────────────┬──────────────┐                                                         
  │ Zoom │ Visible half-extent │ Grid spacing │
  ├──────┼─────────────────────┼──────────────┤                                                         
  │ 10   │ ~58 km              │ ~468 m       │
  ├──────┼─────────────────────┼──────────────┤
  │ 14   │ ~3.7 km             │ ~29 m        │                                                         
  ├──────┼─────────────────────┼──────────────┤
  │ 18   │ ~228 m              │ ~1.8 m       │                                                         
  ├──────┼─────────────────────┼──────────────┤                                                         
  │ 19   │ ~114 m              │ ~1 m         │
  ├──────┼─────────────────────┼──────────────┤                                                         
  │ 20   │ ~57 m               │ ~1 m         │
  └──────┴─────────────────────┴──────────────┘                                                         

### Added
- Lots of undocumented things....

### Added
- Sensor and target placement via graph clicks, with unlimited points supported
- Toggle buttons to switch between sensor placement mode and target placement mode, with visual active/inactive state (green = active)
- Clear All button to remove all placed sensors and targets at once
- Live tables below the graph showing all placed sensors and targets with their X/Y coordinates
- Individual Remove buttons on each table row to delete a single sensor or target
- State managed via `dcc.Store` (replacing global variables) for correctness across callbacks
- Calculate Geolocation button that invokes a stub function `calculate_geolocation(sensors, targets)`, ready to be wired to an external algorithm
- Result display area below the tables showing the estimated location or an error/warning message

