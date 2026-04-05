import numpy as np
import plotly.express as px
import dash
from dash import Input, Output, State, ALL, ctx, dcc, html
import dash_bootstrap_components as dbc
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import random as rand
import itertools

import coloredlogs, logging
logger = logging.getLogger()#__name__)
coloredlogs.install(level="DEBUG", fmt='%(asctime)s %(name)s %(levelname)s %(message)s')

from geolocation import *

GS = 100                        # snap grid resolution (points per axis)
GRID_MIN, GRID_MAX = 0, 100_000  # m — default / initial view range (100 km)

# Ray length for DOA lines: 1.5× the initial grid diagonal guarantees overshoot
# within the default view; Plotly clips lines that extend past the axis range.
DOA_LINE_LENGTH = (GRID_MAX - GRID_MIN) * 1.5

# NOTE: hint that red values in the UI mean that they are outside the bounds set in this code
#   e.g. if you try to set y=100001 it will show as red and probably will revert to whatever it was before.
# NOTE: It's also very finicky the way it works...
#   e.g. min=0, max=99999, step=1 means you can only have integers between 0 and 99999 inclusive.
#   e.g. min=0, max=99999, step=0.1 means you can now have fractions with one decimal place.
# NOTE: Also when x,y locations are first populated in the table they may start out red b/c the resolution of the table grid is higher than the
#       restrictions in the table... so you get a weird situation where you can place a point at x=35123 but you couldn't edit it to be 35123.5

def make_fig(sensors=None, targets=None, ellipses=None, doa_data=None, hyperbola_data=None, view_range=None):
    if view_range is None:
        view_range = {'x': [GRID_MIN, GRID_MAX], 'y': [GRID_MIN, GRID_MAX]}
    vx, vy = view_range['x'], view_range['y']

    fig = make_subplots(rows=1, cols=1)

    # Invisible snap grid — regenerated to cover the current view so clicks
    # always land on a point within the visible area.
    # hoverinfo='none' suppresses the hover label without disabling click/zoom detection.
    # (hoverinfo='skip' would also disable clicks, breaking sensor/target placement.)
    fig.add_traces(
        px.scatter(x=np.repeat(np.linspace(vx[0], vx[1], GS), GS),
                   y=np.tile(np.linspace(vy[0], vy[1], GS), GS))
        .update_traces(marker_color="rgba(0,0,0,0)", hoverinfo='none')
        .data)

    if sensors:
        fig.add_trace(go.Scatter(
            x=[s['x'] for s in sensors], y=[s['y'] for s in sensors],
            mode='markers', marker=dict(color='rgba(0,0,0,1)', size=10, symbol='circle'),
            name='Sensors'))
    if targets:
        fig.add_trace(go.Scatter(
            x=[t['x'] for t in targets], y=[t['y'] for t in targets],
            mode='markers', marker=dict(color='rgba(152,0,0,0.8)', size=10, symbol='star'),
            name='Targets'))
    if doa_data and sensors:
        main_xs, main_ys = [], []
        bound_xs, bound_ys = [], []
        fill_xs,  fill_ys  = [], []
        ray_len = max(vx[1] - vx[0], vy[1] - vy[0]) * 1.5
        for t_idx, target_doas in enumerate(doa_data):
            for s_idx, doa_deg in enumerate(target_doas):
                sx, sy = sensors[s_idx]['x'], sensors[s_idx]['y']
                sigma_rad = np.radians(sensors[s_idx].get('sigma', 1.0))
                angle_rad = np.radians(doa_deg)

                ex       = sx + ray_len * np.cos(angle_rad)
                ey       = sy + ray_len * np.sin(angle_rad)
                ex_upper = sx + ray_len * np.cos(angle_rad + sigma_rad)
                ey_upper = sy + ray_len * np.sin(angle_rad + sigma_rad)
                ex_lower = sx + ray_len * np.cos(angle_rad - sigma_rad)
                ey_lower = sy + ray_len * np.sin(angle_rad - sigma_rad)

                main_xs  += [sx, ex, None]
                main_ys  += [sy, ey, None]
                bound_xs += [sx, ex_upper, None, sx, ex_lower, None]
                bound_ys += [sy, ey_upper, None, sy, ey_lower, None]
                # Wedge polygon: sensor -> upper tip -> lower tip -> close
                fill_xs  += [sx, ex_upper, ex_lower, None]
                fill_ys  += [sy, ey_upper, ey_lower, None]

        # Draw fill first so lines render on top
        fig.add_trace(go.Scatter(
            x=fill_xs, y=fill_ys,
            mode='lines', line=dict(width=0),
            fill='toself', fillcolor='rgba(255, 140, 0, 0.12)',
            name='DOA ±σ Region', showlegend=True))
        fig.add_trace(go.Scatter(
            x=bound_xs, y=bound_ys,
            mode='lines', line=dict(color='rgba(255, 140, 0, 0.6)', width=1, dash='dash'),
            name='DOA ±σ Bounds', showlegend=True))
        fig.add_trace(go.Scatter(
            x=main_xs, y=main_ys,
            mode='lines', line=dict(color='rgba(255, 140, 0, 0.85)', width=1.5),
            name='DOA Lines', showlegend=True))
    if hyperbola_data:
        hyp_colors = [
            ('rgba(0,150,255,0.75)',  'rgba(0,150,255,0.35)',  'rgba(0,150,255,0.10)'),
            ('rgba(0,200,100,0.75)',  'rgba(0,200,100,0.35)',  'rgba(0,200,100,0.10)'),
            ('rgba(200,0,200,0.75)',  'rgba(200,0,200,0.35)',  'rgba(200,0,200,0.10)'),
            ('rgba(255,100,0,0.75)',  'rgba(255,100,0,0.35)',  'rgba(255,100,0,0.10)'),
        ]
        for t_idx, target_hyperbolas in enumerate(hyperbola_data):
            col_solid, col_dash, col_fill = hyp_colors[t_idx % len(hyp_colors)]
            for h_idx, hyp in enumerate(target_hyperbolas):
                pts_c = np.array(hyp['center'])
                pts_u = np.array(hyp['upper'])
                pts_l = np.array(hyp['lower'])
                first = (h_idx == 0)

                # Fill between upper and lower bounds
                if pts_u.shape[0] > 0 and pts_l.shape[0] > 0:
                    fill_x = pts_u[:, 0].tolist() + pts_l[::-1, 0].tolist()
                    fill_y = pts_u[:, 1].tolist() + pts_l[::-1, 1].tolist()
                    fig.add_trace(go.Scatter(
                        x=fill_x, y=fill_y,
                        mode='lines', line=dict(width=0),
                        fill='toself', fillcolor=col_fill,
                        legendgroup=f'tdoa_t{t_idx}', showlegend=False,
                    ))

                # ±σ bound hyperbolas — dashed
                for pts_b in (pts_u, pts_l):
                    if pts_b.shape[0] == 0:
                        continue
                    fig.add_trace(go.Scatter(
                        x=pts_b[:, 0].tolist(), y=pts_b[:, 1].tolist(),
                        mode='lines',
                        line=dict(color=col_dash, width=1, dash='dash'),
                        name=f'TDOA ±σ Bounds (T{t_idx + 1})',
                        legendgroup=f'tdoa_t{t_idx}',
                        showlegend=False,
                    ))

                # Central hyperbola — solid
                if pts_c.shape[0] == 0:
                    continue
                fig.add_trace(go.Scatter(
                    x=pts_c[:, 0].tolist(), y=pts_c[:, 1].tolist(),
                    mode='lines',
                    line=dict(color=col_solid, width=1.5),
                    name=f'TDOA Hyperbola (T{t_idx + 1})',
                    legendgroup=f'tdoa_t{t_idx}',
                    showlegend=first,
                ))
    if ellipses:
        t = np.linspace(0, 2 * np.pi, 200)
        for e in ellipses:
            label = f"Target {e['target_idx'] + 1}"
            # Ellipse from shape points returned by geolocate
            fig.add_trace(go.Scatter(
                x=e['shape_x'], y=e['shape_y'],
                mode='lines', line=dict(color='rgba(152,0,0,0.8)', width=2),
                fill='toself', fillcolor='rgba(152,0,0,0.1)',
                name=f"Ellipse shape ({label})",
                showlegend=True))
            # Parametric ellipse computed from center + axes + orientation
            a, b, theta = e['semimajor'], e['semiminor'], e['orientation']
            px_vals = e['cx'] + a * np.cos(t) * np.cos(theta) - b * np.sin(t) * np.sin(theta)
            py_vals = e['cy'] + a * np.cos(t) * np.sin(theta) + b * np.sin(t) * np.cos(theta)
            fig.add_trace(go.Scatter(
                x=px_vals.tolist(), y=py_vals.tolist(),
                mode='lines', line=dict(color='rgba(0,100,200,0.8)', width=2, dash='dash'),
                fill='toself', fillcolor='rgba(0,100,200,0.1)',
                name=f"Ellipse parametric ({label})",
                showlegend=True))
    fig.update_layout(
        xaxis_title="X (m)",
        yaxis_title="Y (m)",
        xaxis=dict(range=vx),
        yaxis=dict(range=vy, scaleanchor='x', scaleratio=1),
        uirevision='constant',
    )
    return fig


def point_in_ellipse(px, py, cx, cy, semimajor, semiminor, orientation):
    """Return True if (px, py) lies inside the rotated ellipse."""
    dx, dy = px - cx, py - cy
    x_local =  dx * np.cos(orientation) + dy * np.sin(orientation)
    y_local = -dx * np.sin(orientation) + dy * np.cos(orientation)
    return (x_local / semimajor) ** 2 + (y_local / semiminor) ** 2 <= 1


def make_table(points, point_type, ellipses=None):
    if not points:
        return html.P(f"No {point_type}s placed.", className="text-muted small")
    is_sensor = (point_type == 'Sensor')
    is_target = (point_type == 'Target')
    # Build ellipse lookup: target_idx -> ellipse
    ellipse_map = {e['target_idx']: e for e in (ellipses or [])}

    header_cells = [html.Th("#"), html.Th("X"), html.Th("Y")]
    if is_sensor:
        header_cells.append(html.Th("σ"))
    if is_target:
        header_cells += [html.Th("Area (m²)"), html.Th("Containment"), html.Th("Inside?")]
    header_cells.append(html.Th(""))

    pt = point_type.lower()
    rows = []
    for i, p in enumerate(points):
        cells = [
            html.Td(f"{point_type} {i + 1}"),
            html.Td(dcc.Input(
                id={'type': f'x-{pt}', 'index': i},
                type='number', value=round(p['x']),
                step=1, debounce=True,
                style={'width': '80px'}
            )),
            html.Td(dcc.Input(
                id={'type': f'y-{pt}', 'index': i},
                type='number', value=round(p['y']),
                step=1, debounce=True,
                style={'width': '80px'}
            )),
        ]
        if is_sensor:
            cells.append(html.Td(dcc.Input(
                id={'type': 'sigma-sensor', 'index': i},
                type='number', value=p.get('sigma', 1.0),
                min=0, max=99, step=0.001, debounce=True,
                style={'width': '70px'}
            )))
        if is_target:
            e = ellipse_map.get(i)
            if e:
                area = np.pi * e['semimajor'] * e['semiminor']
                inside = point_in_ellipse(
                    p['x'], p['y'],
                    e['cx'], e['cy'],
                    e['semimajor'], e['semiminor'], e['orientation']
                )
                cells += [
                    html.Td(f"{area:.2f}"),
                    html.Td(f"{e['containment']:.0%}"),
                    html.Td(
                        html.Span("Yes", style={'color': 'green', 'fontWeight': 'bold'})
                        if inside else
                        html.Span("No", style={'color': 'red', 'fontWeight': 'bold'})
                    ),
                ]
            else:
                cells += [html.Td("—"), html.Td("—"), html.Td("—")]
        cells.append(html.Td(dbc.Button(
            "Remove",
            id={'type': f'del-{point_type.lower()}', 'index': i},
            size='sm', color='danger', outline=True, n_clicks=0
        )))
        rows.append(html.Tr(cells))

    return dbc.Table(
        [html.Thead(html.Tr(header_cells)), html.Tbody(rows)],
        bordered=True, hover=True, size='sm'
    )



def get_tdoa_hyperbola(locA, toaA, locB, toaB, hyperbola_length):
    """
    Compute the TDOA hyperbola curve for a pair of sensors A and B.

    The locus of points where |dAE - dBE| = c|TDOA| is one branch of a hyperbola
    with foci at A and B.  This function uses the standard cosh/sinh parametric
    form, which is always valid regardless of sensor/target geometry.

    locA            = sensor A position [x, y]
    toaA            = time of arrival at sensor A (s)
    locB            = sensor B position [x, y]
    toaB            = time of arrival at sensor B (s)
    hyperbola_length = controls how far along the curve to draw; the curve
                       extends to cosh(t_max) = hyperbola_length multiples of
                       the semi-transverse axis from the hyperbola vertex.

    Returns (TDOA, locE) where locE is an Nx2 array of [x, y] points tracing
    the relevant branch, or an empty (0x2) array if the geometry is degenerate.
    """
    c = 299792458   # speed of RF wave through vacuum

    TDOA = toaB - toaA
    D    = c * TDOA         # signed range difference: dBE - dAE = D

    # --- Baseline geometry ---
    diff = locB - locA
    dAB  = np.sqrt(diff[0]**2 + diff[1]**2)

    a = abs(D) / 2          # semi-transverse axis  (half the range difference)
    f = dAB / 2             # focal half-distance

    print(f"\n>> get_tdoa_hyperbola  TDOA={TDOA:.3e}s  D={D:.3f}  dAB={dAB:.3f}  a={a:.3f}  f={f:.3f}")

    # Degenerate: range difference ≥ baseline → emitter at or beyond a sensor on
    # the line connecting them.  No valid hyperbola branch exists.
    if a >= f:
        print("  WARNING: |D| >= dAB — degenerate geometry, returning empty locE")
        return TDOA, np.empty((0, 2))

    b = np.sqrt(f**2 - a**2)   # semi-conjugate axis

    # --- Coordinate frame centred on the midpoint of AB ---
    M = (locA + locB) / 2
    u = diff / dAB              # unit vector A → B
    v = np.array([-u[1], u[0]]) # perpendicular unit vector

    # --- Parametric sweep ---
    # t_max chosen so the tip-to-end distance in the local x direction is
    # a * cosh(t_max) = a * hyperbola_length, i.e. t_max = arccosh(length).
    t_max = np.arccosh(max(float(hyperbola_length), 1.001))
    N     = max(int(100 * hyperbola_length), 200)
    t     = np.linspace(-t_max, t_max, N)

    # Branch selection: D > 0 → toaB > toaA → signal reached A first → E is
    # on the A-side branch (negative local-x, i.e. x_local = -a·cosh(t)).
    branch_sign = -np.sign(D) if D != 0 else 1.0
    x_local = branch_sign * a * np.cosh(t)
    y_local = b * np.sinh(t)

    # Transform back to global coordinates
    locE = M + np.outer(x_local, u) + np.outer(y_local, v)

    print(f"  branch_sign={branch_sign}  a={a:.3f}  b={b:.3f}  t_max={t_max:.3f}  N={N}")
    print(f"  locE[0]={locE[0]}  locE[-1]={locE[-1]}")
    return TDOA, locE


def sim_tdoa():
    # TO DRAW IT
    c = 299792458   # speed of RF wave through vacuum (close enough)
    sigma_time = 0.000001
    length = 10  # how long the line should be, in terms of the # of toaE's lengths..... hard to quantify but 10 gives a good picture for 2 points, would need to calculate this based on the overall scenario map limits

    ########################################
    loc1 = np.array([-18520, 37040])
    loc2 = np.array([0, 37040])
    e1 = np.array([ 5556, 74080])
    sensors = [{'x': loc1[0], 'y': loc1[1], 'sigma': sigma_time}, {'x': loc2[0], 'y': loc2[1], 'sigma': sigma_time}]
    targets = [{'x': e1[0], 'y': e1[1]}]

    for target in targets:
        loc_array = []
        sigma_array = []
        doa_array = []

        d1e = np.sqrt(np.power(loc1[0]-e1[0],2)+np.power(loc1[1]-e1[1],2))
        d2e = np.sqrt(np.power(loc2[0]-e1[0],2)+np.power(loc2[1]-e1[1],2))
        ########################################

        t1 = 0              # TOA at sensor 1 is our reference point
        te = t1 - d1e/c     # emission time will be negative since we're using TOA at sensor 1 for the zero reference point
        t1 = te + d1e/c# + rand.normalvariate(mu=0.0, sigma=sigma_time)

        ########################################
        t2 = te + d2e/c# + rand.normalvariate(mu=0.0, sigma=sigma_time)
        tdoa,hyperbola = get_tdoa_hyperbola(loc1, t1, loc2, t2, length)
        print("tdoa: "+str(tdoa))
        print("")
        ########################################

        tns = [t1]
        loc_array.append(np.vstack([[sensors[0]['x']],[sensors[0]['y']]]))
        sigma_array.append(sigma_time)
        for sensor in sensors[1:]:
            loc_array.append(np.vstack([[sensor['x']],[sensor['y']]]))
            sigma_array.append(sigma_time)

            dne = np.sqrt(np.power(sensor['x']-target['x'],2)+np.power(sensor['y']-target['y'],2))
            tn = te + dne/c# + rand.normalvariate(mu=0.0, sigma=sigma_time)
            tns.append(tn)

        # Get TDOA calculation and hyperbola for every combination of sensors
        tdoas = np.array([])
        indices_combinations = itertools.combinations(range(len(sensors)), 2)
        for indices_tuple in indices_combinations:
            items_tuple = (sensors[i] for i in indices_tuple)
            items = tuple(items_tuple)
            print(f"Indices: {indices_tuple}, Items: {items}")

        for s1, s2 in itertools.combinations(sensors, 2):
            print("s1: "+str(items[0])+", t1: "+str(tns[indices_tuple[0]]))
            print("s2: "+str(items[1])+", t2: "+str(tns[indices_tuple[1]]))
            loc1 = np.array([items[0]['x'], items[0]['y']])
            loc2 = np.array([items[1]['x'], items[1]['y']])
            tdoa,hyperbola = get_tdoa_hyperbola(loc1, tns[indices_tuple[0]], loc2, tns[indices_tuple[1]], length)
            print("tdoa: "+str(tdoa))
            tdoas = np.append(tdoas,tdoa)
        print("tdoas: "+str(tdoas))

        # z = theta basically but for non-DOA 
        # z = np.array([(c*np.abs(tdoa21)), (c*np.abs(tdoa31)), (c*np.abs(tdoa32))])
        doa_array = c*np.abs(tdoas)
        print("doa_array: "+str(doa_array))


def calculate_measurements(sensors, targets, mode='tdoa'):
    """Simulate DOA or TDOA measurements for all targets.

    Args:
        sensors: list of {'x', 'y', 'sigma'} dicts.
        targets: list of {'x', 'y'} dicts.
        mode:    'doa' or 'tdoa'.

    Returns:
        measurements:    list (one per target) of dicts with keys
                         'loc_array', 'doa_array', 'sigma_array' — all JSON-serializable.
        hyperbola_data:  hyperbola_data[target_idx][pair_idx] = [[x, y], ...]
                         (empty lists for DOA mode)
        doa_data:        doa_data[target_idx][sensor_idx] = DOA angle in degrees
                         (empty lists for TDOA mode)
    """
    logger.debug(f">calculate_measurements mode={mode}")

    measurements = []
    hyperbola_data = []
    doa_data = []

    for target in targets:
        logger.debug("  |-Target: "+str(target))

        loc_array = []
        sigma_array = []
        doa_array = []
        doa_per_sensor = []
        hyperbolas_this_target = []

        if mode == 'doa':
            ###################################################################
            # DOA simulation
            # Each sensor produces one angle measurement with Gaussian noise.
            ###################################################################
            for sensor in sensors:
                loc_array.append(np.vstack([[sensor['x']], [sensor['y']]]))

                # sigma from the sensor table is in degrees; convert to radians
                sigma_array.append(sensor['sigma'] * np.pi / 180)

                # True angle from sensor to target (zero = due East, CCW positive)
                theta = np.arctan2(target['y'] - sensor['y'],
                                   target['x'] - sensor['x'])

                # Add Gaussian noise scaled by sensor sigma (degrees)
                error = np.random.normal(0, sensor['sigma'])
                doa = theta + error * np.pi / 180
                doa_array.append(float(doa))
                doa_per_sensor.append(float(np.degrees(doa)))

                logger.debug("    |-Sensor: "+str(sensor))
                logger.debug(f"    |--DOA:   {np.degrees(doa):.3f}°")

        else:
            ###################################################################
            # TDOA simulation
            # Sensor 0 is the time reference; all others are measured relative.
            ###################################################################
            c = 299792458   # speed of RF wave through vacuum
            # Sensor sigma is entered in the UI in microseconds for TDOA mode
            length = 10     # hyperbola sweep length (multiples of initial toaE)

            # --- Phase 1: simulate a noisy TOA for every sensor ---
            d1e = np.sqrt(np.power(sensors[0]['x'] - target['x'], 2) +
                          np.power(sensors[0]['y'] - target['y'], 2))
            te = -d1e / c   # true emission time (sensor 0 TOA = 0 reference)
            tns = []
            for sensor in sensors:
                sigma_n = sensor['sigma'] * 1e-6        # Treat sigma value set by user as microseconds
                dne = np.sqrt(np.power(sensor['x'] - target['x'], 2) +
                              np.power(sensor['y'] - target['y'], 2))
                tns.append(te + dne / c + np.random.normal(0, sigma_n))

            # --- Phase 2: one measurement entry per sensor pair ---
            # loc_array, sigma_array, and doa_array are all built here so their
            # lengths stay in sync: one entry per C(N,2) pair.
            tdoas = np.array([])
            for i, j in itertools.combinations(range(len(sensors)), 2):
                loc1 = np.array([sensors[i]['x'], sensors[i]['y']])
                loc2 = np.array([sensors[j]['x'], sensors[j]['y']])
                sigma_i = sensors[i]['sigma'] * 1e-6
                sigma_j = sensors[j]['sigma'] * 1e-6
                print(f"s{i}: {sensors[i]}, t{i}: {tns[i]}")
                print(f"s{j}: {sensors[j]}, t{j}: {tns[j]}")
                tdoa, hyp_center = get_tdoa_hyperbola(loc1, tns[i], loc2, tns[j], length)
                sigma_combined = np.sqrt(sigma_i**2 + sigma_j**2)
                _,    hyp_upper  = get_tdoa_hyperbola(loc1, tns[i], loc2, tns[j] + sigma_combined, length)
                _,    hyp_lower  = get_tdoa_hyperbola(loc1, tns[i], loc2, tns[j] - sigma_combined, length)
                print("tdoa: "+str(tdoa))
                tdoas = np.append(tdoas, tdoa)
                hyperbolas_this_target.append({
                    'center': hyp_center.tolist(),
                    'upper':  hyp_upper.tolist(),
                    'lower':  hyp_lower.tolist(),
                })
                # midpoint of the pair as the representative sensor position,
                # combined timing sigma via root-sum-square
                midpoint = np.vstack([[(sensors[i]['x'] + sensors[j]['x']) / 2],
                                      [(sensors[i]['y'] + sensors[j]['y']) / 2]])
                loc_array.append(midpoint)
                sigma_array.append(np.sqrt(sigma_i**2 + sigma_j**2))
            print("tdoas: "+str(tdoas))

            doa_array = (c * np.abs(tdoas)).tolist()
            print("doa_array: "+str(doa_array))

        measurements.append({
            'loc_array': [a.tolist() for a in loc_array],
            'doa_array': doa_array,
            'sigma_array': sigma_array,
        })
        hyperbola_data.append(hyperbolas_this_target)
        doa_data.append(doa_per_sensor)

    return measurements, hyperbola_data, doa_data


def calculate_geolocation(measurements, targets, containment=0.95):
    """Compute geolocation ellipses from stored measurements.

    Args:
        measurements: list of per-target dicts produced by calculate_measurements.
        targets:      list of {'x', 'y'} dicts (used only to tag ellipses).
        containment:  probability containment fraction (e.g. 0.95).

    Returns:
        list of ellipse dicts from geolocate(), one per target.
    """
    logger.debug(">calculate_geolocation")

    ellipses = []
    for m, target in zip(measurements, targets):
        loc_array   = [np.array(a) for a in m['loc_array']]
        doa_array   = np.array(m['doa_array'])
        sigma_array = m['sigma_array']

        ellipse = geolocate(loc_array, doa_array, sigma_array, containment)
        logger.debug("  |-x,y:         "+str(ellipse['x'][0])+", "+str(ellipse['y'][0]))
        logger.debug("  |-semimajor:   "+str(ellipse['semimajor']))
        logger.debug("  |-semiminor:   "+str(ellipse['semiminor']))
        logger.debug("  |-orientation: "+str(ellipse['orientation']*180/np.pi))
        ellipse['target_id'] = target
        ellipses.append(ellipse)

    return ellipses


app = dash.Dash(__name__, title='Geo', external_stylesheets=[dbc.themes.BOOTSTRAP])

# Default x range is pre-widened to match a 3:2 aspect ratio (y stays 0–100 km,
# x extends 25% on each side). relayoutData will correct this if the actual
# rendered aspect ratio differs.
_span = GRID_MAX - GRID_MIN
_default_view = {'x': [GRID_MIN - _span * 0.25, GRID_MAX + _span * 0.25],
                 'y': [GRID_MIN, GRID_MAX]}

app.layout = html.Div([
    dcc.Store(id='store', data={'sensors': [], 'targets': [], 'measurements': [], 'ellipses': [], 'doa_data': [], 'hyperbola_data': [], 'show_doa': False}),
    dcc.Store(id='mode-store', data='sensor'),
    dcc.Store(id='view-range', data=_default_view),

    dcc.Graph(id='graph', figure=make_fig(), style={'width': '100%', 'aspectRatio': '3 / 1'},
              config={'scrollZoom': True}),

    html.Div([
        dbc.Row([
            dbc.Col(dbc.Button('Clear All', id='btn-clear', n_clicks=0,
                               color='secondary', outline=True), width='auto'),
            dbc.Col(dbc.Button('Place Sensors', id='btn-place-sensor', n_clicks=0,
                               color='success'), width='auto'),
            dbc.Col([
                dbc.InputGroup([
                    dcc.Input(id='sigma-input', type='number', value=3,
                              min=0, max=99, step=0.001, style={'width': '60px'},
                              className='form-control form-control-sm'),
                    dbc.InputGroupText("σ"),
                ], size='sm'),
            ], width='auto', className='align-self-center'),
            dbc.Col(dbc.Button('Place Targets', id='btn-place-target', n_clicks=0,
                               color='secondary'), width='auto'),
            dbc.Col([
                dbc.InputGroup([
                    dcc.Input(id='containment-input', type='number', value=95,
                              min=1, max=99, step=1, style={'width': '60px'},
                              className='form-control form-control-sm'),
                    dbc.InputGroupText("%"),
                ], size='sm'),
            ], width='auto', className='align-self-center'),
            dbc.Col(
                dbc.RadioItems(
                    id='measure-mode',
                    options=[{'label': 'DOA', 'value': 'doa'},
                             {'label': 'TDOA', 'value': 'tdoa'}],
                    value='tdoa',
                    inline=True,
                ),
                width='auto', className='align-self-center'
            ),
            dbc.Col(dbc.Button('Calculate Measurements', id='btn-measure', n_clicks=0,
                               color='warning'), width='auto'),
            dbc.Col(dbc.Button('Calculate Geolocation', id='btn-calculate', n_clicks=0,
                               color='primary', disabled=True), width='auto'),
            dbc.Col(dbc.Button('Show Measurements', id='btn-show-doa', n_clicks=0,
                               color='info', disabled=True), width='auto'),
        ], className='mt-2 mb-3 g-2'),

        dbc.Row([
            dbc.Col([html.H5('Sensors'), html.Div(id='sensors-table')], width=6),
            dbc.Col([html.H5('Targets'), html.Div(id='targets-table')], width=6),
        ]),

        dbc.Row([
            dbc.Col(html.Div(id='geo-result'), className='mt-3'),
        ]),
    ], className='px-3 pb-3'),
])


@app.callback(
    Output('view-range', 'data'),
    Input('graph', 'relayoutData'),
    State('view-range', 'data'),
    prevent_initial_call=True,
)
def update_view_range(relayout_data, current_range):
    if not relayout_data:
        raise dash.exceptions.PreventUpdate
    if 'xaxis.range[0]' in relayout_data:
        x0, x1 = relayout_data['xaxis.range[0]'], relayout_data['xaxis.range[1]']
        y0 = relayout_data.get('yaxis.range[0]', x0)
        y1 = relayout_data.get('yaxis.range[1]', x1)
        return {'x': [x0, x1], 'y': [y0, y1]}
    if 'xaxis.autorange' in relayout_data or 'autosize' in relayout_data:
        return _default_view
    raise dash.exceptions.PreventUpdate


@app.callback(
    Output('store', 'data'),
    Output('geo-result', 'children'),
    Input('graph', 'clickData'),
    Input('btn-clear', 'n_clicks'),
    Input('btn-measure', 'n_clicks'),
    Input('btn-calculate', 'n_clicks'),
    Input('btn-show-doa', 'n_clicks'),
    Input({'type': 'del-sensor', 'index': ALL}, 'n_clicks'),
    Input({'type': 'del-target', 'index': ALL}, 'n_clicks'),
    Input({'type': 'sigma-sensor', 'index': ALL}, 'value'),
    Input({'type': 'x-sensor', 'index': ALL}, 'value'),
    Input({'type': 'y-sensor', 'index': ALL}, 'value'),
    Input({'type': 'x-target', 'index': ALL}, 'value'),
    Input({'type': 'y-target', 'index': ALL}, 'value'),
    State('store', 'data'),
    State('mode-store', 'data'),
    State('sigma-input', 'value'),
    State('containment-input', 'value'),
    State('measure-mode', 'value'),
)
def update_store(click_data, _clear, _measure, _calc, _show_doa, _del_s, _del_t,
                 sigma_sensor_vals, x_sensor_vals, y_sensor_vals, x_target_vals, y_target_vals,
                 store, mode, sigma_val, containment_pct, measure_mode):
    triggered = ctx.triggered_id
    no_alert = dash.no_update

    if triggered == 'btn-clear':
        return {'sensors': [], 'targets': [], 'measurements': [], 'ellipses': [], 'doa_data': [], 'hyperbola_data': [], 'show_doa': False}, no_alert

    if triggered == 'btn-show-doa':
        store['show_doa'] = not store.get('show_doa', False)
        return store, no_alert

    if isinstance(triggered, dict) and triggered['type'] == 'del-sensor':
        idx = triggered['index']
        pruned_doa = [
            [doa for s_idx, doa in enumerate(target_doas) if s_idx != idx]
            for target_doas in store.get('doa_data', [])
        ]
        return {'sensors': [s for i, s in enumerate(store['sensors']) if i != idx],
                'targets': store['targets'],
                'ellipses': store.get('ellipses', []),
                'doa_data': pruned_doa,
                'show_doa': store.get('show_doa', False)}, no_alert

    if isinstance(triggered, dict) and triggered['type'] == 'del-target':
        idx = triggered['index']
        # Remove ellipses for this target; shift target_idx for remaining ones
        remaining = []
        for e in store.get('ellipses', []):
            if e['target_idx'] == idx:
                continue
            e = dict(e)
            if e['target_idx'] > idx:
                e['target_idx'] -= 1
            remaining.append(e)
        remaining_doa = [d for i, d in enumerate(store.get('doa_data', [])) if i != idx]
        return {'sensors': store['sensors'],
                'targets': [t for i, t in enumerate(store['targets']) if i != idx],
                'ellipses': remaining,
                'doa_data': remaining_doa,
                'show_doa': store.get('show_doa', False)}, no_alert

    if isinstance(triggered, dict) and triggered['type'] == 'sigma-sensor':
        idx = triggered['index']
        val = sigma_sensor_vals[idx]
        if val is not None:
            store['sensors'][idx]['sigma'] = float(val)
        return store, no_alert

    if isinstance(triggered, dict) and triggered['type'] == 'x-sensor':
        idx = triggered['index']
        val = x_sensor_vals[idx]
        if val is not None:
            store['sensors'][idx]['x'] = float(val)
        return store, no_alert

    if isinstance(triggered, dict) and triggered['type'] == 'y-sensor':
        idx = triggered['index']
        val = y_sensor_vals[idx]
        if val is not None:
            store['sensors'][idx]['y'] = float(val)
        return store, no_alert

    if isinstance(triggered, dict) and triggered['type'] == 'x-target':
        idx = triggered['index']
        val = x_target_vals[idx]
        if val is not None:
            store['targets'][idx]['x'] = float(val)
        return store, no_alert

    if isinstance(triggered, dict) and triggered['type'] == 'y-target':
        idx = triggered['index']
        val = y_target_vals[idx]
        if val is not None:
            store['targets'][idx]['y'] = float(val)
        return store, no_alert

    if triggered == 'graph' and click_data:
        pt = click_data['points'][0]
        entry = {'x': pt['x'], 'y': pt['y']}
        if mode == 'sensor':
            entry['sigma'] = float(sigma_val) if sigma_val and sigma_val > 0 else 1.0
            store['sensors'].append(entry)
        else:
            store['targets'].append(entry)
        return store, no_alert

    if triggered == 'btn-measure':
        sensors = store['sensors']
        targets = store['targets']
        if len(sensors) < 2:
            return store, dbc.Alert("Add at least two sensors before calculating measurements.", color='warning')
        if not targets:
            return store, dbc.Alert("Add at least one target before calculating measurements.", color='warning')
        measurements, hyperbola_data, doa_data = calculate_measurements(sensors, targets, mode=measure_mode or 'tdoa')
        store['measurements']    = measurements
        store['hyperbola_data']  = hyperbola_data
        store['doa_data']        = doa_data
        store['ellipses']        = []   # clear stale ellipses when measurements are refreshed
        alert = dbc.Alert(f"Calculated measurements for {len(measurements)} target(s).", color='info')
        return store, alert

    if triggered == 'btn-calculate':
        measurements = store.get('measurements', [])
        targets = store['targets']
        if not measurements:
            return store, dbc.Alert("Calculate measurements first.", color='warning')
        containment = (containment_pct or 95) / 100.0
        results = calculate_geolocation(measurements, targets, containment)
        ellipses = []
        for target_idx, ellipse in enumerate(results):
            shape = np.rot90(ellipse['shape'])      # rotate 90 degrees, original is [[x1 y1], [x2 y2], ...], result is [[y1 y2 ...], [x1 x2 ...]]
            ellipses.append({
                'shape_x': shape[1].tolist(),
                'shape_y': shape[0].tolist(),
                'cx': float(ellipse['x'][0]),
                'cy': float(ellipse['y'][0]),
                'semimajor': float(ellipse['semimajor']),
                'semiminor': float(ellipse['semiminor']),
                'orientation': float(ellipse['orientation']),
                'containment': containment,
                'target_idx': target_idx,
            })
        store['ellipses'] = ellipses
        alert = dbc.Alert(f"Calculated {len(ellipses)} ellipse(s).", color='success')
        return store, alert

    raise dash.exceptions.PreventUpdate


@app.callback(
    Output('graph', 'figure'),
    Output('sensors-table', 'children'),
    Output('targets-table', 'children'),
    Output('btn-show-doa', 'disabled'),
    Output('btn-show-doa', 'children'),
    Input('store', 'data'),
    Input('view-range', 'data'),
)
def update_view(store, view_range):
    sensors = store['sensors']
    targets = store['targets']
    ellipses = store.get('ellipses', [])
    show_doa = store.get('show_doa', False)
    doa_data = store.get('doa_data', []) if show_doa else None
    hyperbola_data = store.get('hyperbola_data', []) if show_doa else None
    btn_disabled = not bool(store.get('measurements'))
    btn_label = 'Hide Measurements' if show_doa else 'Show Measurements'
    return (make_fig(sensors, targets, ellipses, doa_data, hyperbola_data, view_range),
            make_table(sensors, 'Sensor'),
            make_table(targets, 'Target', ellipses),
            btn_disabled,
            btn_label)


@app.callback(
    Output('btn-calculate', 'disabled'),
    Input('store', 'data'),
)
def update_calculate_button(store):
    return not bool(store.get('measurements'))


@app.callback(
    Output('mode-store', 'data'),
    Output('btn-place-sensor', 'color'),
    Output('btn-place-target', 'color'),
    Input('btn-place-sensor', 'n_clicks'),
    Input('btn-place-target', 'n_clicks'),
)
def update_mode(_s, _t):
    if ctx.triggered_id == 'btn-place-target':
        return 'target', 'secondary', 'success'
    return 'sensor', 'success', 'secondary'


if __name__ == '__main__':
    app.run(debug=True)
