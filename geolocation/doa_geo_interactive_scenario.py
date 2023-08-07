import numpy as np
import plotly.express as px
import math, json
import dash
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# Add invisible position grid
GS = 100
_fig = make_subplots(rows=1, cols=1)
_fig.add_traces(
    px.scatter(x=np.repeat(np.linspace(0, 1, GS), GS), y=np.tile(np.linspace(0, 1, GS), GS))
    .update_traces(marker_color="rgba(0,0,0,0)")
    .data)


_place_sensors = True
_sensor_xpoints = np.array([])
_sensor_ypoints = np.array([])
_target_xpoints = np.array([])
_target_ypoints = np.array([])


active_button_style = {'background-color': 'red',
                    'color': 'white',
                    'height': '50px',
                    'width': '100px',
                    'margin-top': '50px',
                    'margin-left': '50px'}

inactive_button_style = {'background-color': 'white',
                      'color': 'black',
                      'height': '50px',
                      'width': '100px',
                      'margin-top': '50px',
                      'margin-left': '50px'}


# Build App
#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
external_stylesheets = [dbc.themes.BOOTSTRAP]

app = dash.Dash(__name__, title='Geo', external_stylesheets=external_stylesheets)
app.layout = dash.html.Div([
    dash.dcc.Graph(id="graph", figure=_fig), dash.html.Div(id="where"),
    dbc.Button('Clear', id='btn-clear', n_clicks=0, className="btn btn-outline-secondary"),
    dbc.Button('Place Sensors', id='btn-place-sensor', n_clicks=0, className="btn btn-outline-secondary"),
    dbc.Button('Place Targets', id='btn-place-target', n_clicks=0, className="btn btn-outline-secondary"),
    dash.html.Br(),
    dash.dcc.Input(id="point-x", type="text", placeholder="", debounce=True),
    dash.dcc.Input(id="point-y", type="text", placeholder="", debounce=True),
    dash.html.Br(),
    dash.html.Div(id='sink')
    ])



@app.callback(
    Output("point-x", "value"),
    Input("graph", "clickData"),
)
def click(clickData):
    if not clickData:
        raise dash.exceptions.PreventUpdate
    return clickData["points"][0]["x"]


@app.callback(
    Output("point-y", "value"),
    Input("graph", "clickData"),
)
def click(clickData):
    if not clickData:
        raise dash.exceptions.PreventUpdate
    return clickData["points"][0]["y"]


@app.callback(
    Output("graph", "figure"),
    Input("graph", "clickData"),
    Input('btn-clear', 'n_clicks')
)
def click(clickData, btn1):
    global _sensor_xpoints, _sensor_ypoints, _target_xpoints, _target_ypoints

    if "btn-clear" == dash.ctx.triggered_id:
        fig = make_subplots(rows=1, cols=1)
        fig.add_traces(
            px.scatter(x=np.repeat(np.linspace(0, 1, GS), GS), y=np.tile(np.linspace(0, 1, GS), GS))
            .update_traces(marker_color="rgba(0,0,0,0)")
            .data)
        return fig

    if not clickData:
        raise dash.exceptions.PreventUpdate

    fig = make_subplots(rows=1, cols=1)
    fig.add_traces(
        px.scatter(x=np.repeat(np.linspace(0, 1, GS), GS), y=np.tile(np.linspace(0, 1, GS), GS))
        .update_traces(marker_color="rgba(0,0,0,0)")
        .data) 

    if _place_sensors:
        _sensor_xpoints = np.append(_sensor_xpoints, clickData["points"][0]["x"])
        _sensor_ypoints = np.append(_sensor_ypoints, clickData["points"][0]["y"])

    else:
        _target_xpoints = np.append(_target_xpoints, clickData["points"][0]["x"])
        _target_ypoints = np.append(_target_ypoints, clickData["points"][0]["y"])

    fig.add_trace(go.Scatter(x=_sensor_xpoints, y=_sensor_ypoints, mode='markers', marker_color='rgba(0, 0, 0, 1)'))
    fig.add_trace(go.Scatter(x=_target_xpoints, y=_target_ypoints, mode='markers', marker_color='rgba(152, 0, 0, .8)'))
    return fig





@app.callback(
    Output('btn-place-sensor', 'className'),
    Input('btn-place-sensor', 'n_clicks'),
    Input('btn-place-target', 'n_clicks')
)
def displayClick(btn1, btn2):
    global _place_sensors

    if "btn-place-sensor" == dash.ctx.triggered_id:
        _place_sensors = True
        return "btn btn-success"

    elif "btn-place-target" == dash.ctx.triggered_id:
        _place_sensors = False
        return "btn btn-secondary"


@app.callback(
    Output('btn-place-target', 'className'),
    Input('btn-place-sensor', 'n_clicks'),
    Input('btn-place-target', 'n_clicks')
)
def displayClick(btn1, btn2):
    global _place_sensors

    if "btn-place-sensor" == dash.ctx.triggered_id:
        _place_sensors = True
        return "btn btn-secondary"

    elif "btn-place-target" == dash.ctx.triggered_id:
        _place_sensors = False
        return "btn btn-success"




if __name__ == '__main__':
    app.run_server(debug=True)