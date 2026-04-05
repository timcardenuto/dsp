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
GRID_MIN, GRID_MAX = 0, 100     # km — default / initial view range

# Ray length for DOA lines: 1.5× the initial grid diagonal guarantees overshoot
# within the default view; Plotly clips lines that extend past the axis range.
DOA_LINE_LENGTH = (GRID_MAX - GRID_MIN) * 1.5

# NOTE: hint that red values in the UI mean that they are outside the bounds set in this code
#   e.g. if you try to set y=101 it will show as red and probably will revert to whatever it was before.
# NOTE: It's also very finicky the way it works...
#   e.g. min=0, max=99, step=1 means you can only have integers between 0 and 99 inclusive. 100 won't work. 43.5 won't work.
#   e.g. min=0, max=99, step=0.1 means you can now have fractions but only one decimal place, 43.5 works but 43.55 does not.
#   e.g. min=0, max=99, step=0.01 means you can now have fractions with two decimal place, 43.55 works but 43.555 does not.
# NOTE: Also when x,y locations are first populated in the table they may start out red b/c the resolution of the table grid is higher than the
#       restrictions in the table... so you get a weird situation where you can place a point at x=35.1234 but you couldn't edit it to be 35.1235

def make_fig(sensors=None, targets=None, ellipses=None, doa_data=None, view_range=None):
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
        xaxis_title="X (km)",
        yaxis_title="Y (km)",
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
        header_cells += [html.Th("Area (km²)"), html.Th("Containment"), html.Th("Inside?")]
    header_cells.append(html.Th(""))

    pt = point_type.lower()
    rows = []
    for i, p in enumerate(points):
        cells = [
            html.Td(f"{point_type} {i + 1}"),
            html.Td(dcc.Input(
                id={'type': f'x-{pt}', 'index': i},
                type='number', value=round(p['x'], 4),
                step=0.0001, debounce=True,
                style={'width': '80px'}
            )),
            html.Td(dcc.Input(
                id={'type': f'y-{pt}', 'index': i},
                type='number', value=round(p['y'], 4),
                step=0.0001, debounce=True,
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

    c = 299792458   # speed of RF wave through vacuum (close enough)
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


def calculate_geolocation(sensors, targets, containment=0.95):
    """Stub: estimate target location from sensor positions.

    Args:
        sensors: list of {'x': float, 'y': float} sensor positions
        targets: list of {'x': float, 'y': float} true/initial target positions

    Returns:
        tuple of (ellipses, doa_data) where doa_data[target_idx][sensor_idx]
        is the DOA measurement in degrees.
    """
    logger.debug(">calculate_geolocation")

    ellipses = []
    doa_data = []   # doa_data[target_idx][sensor_idx] = DOA in degrees
    for target in targets:
        logger.debug("  |-Target: "+str(target))
        loc_array = []
        sigma_array = []
        doa_array = []
        doa_per_sensor = []
        if False:
            for sensor in sensors:
                loc_array.append(np.vstack([[sensor['x']],[sensor['y']]]))

                ###################################
                # Simulation for DOA measurements
                # Convert sigma from degrees to radians
                sigma_array.append(sensor['sigma']*np.pi/180 )
                # Calculate relative angle between target and measurement location (zero angle is when target is due East from sensor, negative is clockwise and positive is counter-clockwise)
                theta = np.arctan2(target['y']-sensor['y'], target['x']-sensor['x'])

                # Add measurement error based on sensor sigma
                # TODO: check the logic here....
                err = 1             # max measurement error in degrees, +-
                error = -err + (err+err)*np.random.rand(1)
                doa = theta + error * np.pi/180
                doa_array.append(doa)
                doa_per_sensor.append(float(np.degrees(doa[0])))

                logger.debug("    |-Sensor: "+str(sensor))
                logger.debug("    |--DOA:   "+str(doa))
                ###################################

        if True:
            # TDOA version
            # for every 2 sensor pairs.... but you have to decide how many sensors you want to be part of the final calculation, I think they all might need to be referenced to a common time=zero sensor
                # sensor and target have to be in np.array([x,y])

            # TO DRAW IT
            c = 299792458   # speed of RF wave through vacuum (close enough)
            sigma_time = 0.000001
            length = 10  # how long the line should be, in terms of the # of toaE's lengths..... hard to quantify but 10 gives a good picture for 2 points, would need to calculate this based on the overall scenario map limits
            
            d1e = np.sqrt(np.power(sensors[0]['x']-target['x'],2)+np.power(sensors[0]['y']-target['y'],2))
            t1 = 0              # TOA at sensor 1 is our reference point
            te = t1 - d1e/c     # emission time will be negative since we're using TOA at sensor 1 for the zero reference point
            t1 = te + d1e/c + rand.normalvariate(mu=0.0, sigma=sigma_time)
            tns = [t1]

            loc_array.append(np.vstack([[sensors[0]['x']],[sensors[0]['y']]]))
            sigma_array.append(sigma_time)
            for sensor in sensors[1:]:
                loc_array.append(np.vstack([[sensor['x']],[sensor['y']]]))
                sigma_array.append(sigma_time)

                dne = np.sqrt(np.power(sensor['x']-target['x'],2)+np.power(sensor['y']-target['y'],2))
                tn = te + dne/c + rand.normalvariate(mu=0.0, sigma=sigma_time)
                tns.append(tn)

            # Get TDOA calculation and hyperbola for every combination of sensors
            tdoas = np.array([])
            indices_combinations = itertools.combinations(range(len(sensors)), 2)
            for indices_tuple in indices_combinations:
                items_tuple = (sensors[i] for i in indices_tuple)
                items = tuple(items_tuple)
                print(f"Indices: {indices_tuple}, Items: {items}")

            # for s1, s2 in itertools.combinations(sensors, 2):
                print("s1: "+str(items[0])+", t1: "+str(tns[indices_tuple[0]]))
                print("s1: "+str(items[1])+", t1: "+str(tns[indices_tuple[1]]))
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


        # Calculate geolocation for the full set of sensors and this specific target
        ellipse = geolocate(loc_array, doa_array, sigma_array, containment)
        logger.debug("  |-Geolocation: ")
        logger.debug("    |-x,y:          "+str(ellipse['x'][0])+", "+str(ellipse['y'][0]))
        logger.debug("    |-semimajor:    "+str(ellipse['semimajor']))
        logger.debug("    |-semiminor:    "+str(ellipse['semiminor']))
        logger.debug("    |-orientation:  "+str(ellipse['orientation']*180/np.pi))
        ellipse['target_id'] = target
        ellipses.append(ellipse)
        doa_data.append(doa_per_sensor)

    return ellipses, doa_data


app = dash.Dash(__name__, title='Geo', external_stylesheets=[dbc.themes.BOOTSTRAP])

# Default x range is pre-widened to match a 3:2 aspect ratio (y stays 0–100,
# x extends 25% on each side). relayoutData will correct this if the actual
# rendered aspect ratio differs.
_span = GRID_MAX - GRID_MIN
_default_view = {'x': [GRID_MIN - _span * 0.25, GRID_MAX + _span * 0.25],
                 'y': [GRID_MIN, GRID_MAX]}

app.layout = html.Div([
    dcc.Store(id='store', data={'sensors': [], 'targets': [], 'ellipses': [], 'doa_data': [], 'show_doa': False}),
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
            dbc.Col(dbc.Button('Calculate Geolocation', id='btn-calculate', n_clicks=0,
                               color='primary'), width='auto'),
            dbc.Col(dbc.Button('Show DOA Lines', id='btn-show-doa', n_clicks=0,
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
)
def update_store(click_data, _clear, _calc, _show_doa, _del_s, _del_t,
                 sigma_sensor_vals, x_sensor_vals, y_sensor_vals, x_target_vals, y_target_vals,
                 store, mode, sigma_val, containment_pct):
    triggered = ctx.triggered_id
    no_alert = dash.no_update

    if triggered == 'btn-clear':
        return {'sensors': [], 'targets': [], 'ellipses': [], 'doa_data': [], 'show_doa': False}, no_alert

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

    if triggered == 'btn-calculate':
        sensors = store['sensors']
        targets = store['targets']
        if len(sensors) < 2:
            return store, dbc.Alert("Add at least two sensors before calculating.", color='warning')
        if not targets:
            return store, dbc.Alert("Add at least one target before calculating.", color='warning')
        containment = (containment_pct or 95) / 100.0
        results, doa_data = calculate_geolocation(sensors, targets, containment)
        ellipses = []
        for target_idx, ellipse in enumerate(results):
            shape = np.rot90(ellipse['shape'])      # rotate 90 degress, original is [[x1 y1], [x2 y2], ...], result is [[y1 y2 ...], [x1 x2 ...]]
            # shape = np.array(ellipse['shape'])    # expected 2xN array: row 0 = x, row 1 = y
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
        store['doa_data'] = doa_data
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
    btn_disabled = not bool(ellipses)
    btn_label = 'Hide DOA Lines' if show_doa else 'Show DOA Lines'
    return (make_fig(sensors, targets, ellipses, doa_data, view_range),
            make_table(sensors, 'Sensor'),
            make_table(targets, 'Target', ellipses),
            btn_disabled,
            btn_label)


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
