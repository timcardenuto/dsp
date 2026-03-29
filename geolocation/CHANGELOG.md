# Changelog

## [Unreleased]

### Updated
- Migrated `doa_geo_interactive_scenario.py` to Dash 4.1.0 syntax:
  - Replaced `from dash.dependencies import Input, Output, State` with `from dash import Input, Output, State`
  - Replaced `app.run_server()` with `app.run()`

### Added
- Sensor and target placement via graph clicks, with unlimited points supported
- Toggle buttons to switch between sensor placement mode and target placement mode, with visual active/inactive state (green = active)
- Clear All button to remove all placed sensors and targets at once
- Live tables below the graph showing all placed sensors and targets with their X/Y coordinates
- Individual Remove buttons on each table row to delete a single sensor or target
- State managed via `dcc.Store` (replacing global variables) for correctness across callbacks
- Calculate Geolocation button that invokes a stub function `calculate_geolocation(sensors, targets)`, ready to be wired to an external algorithm
- Result display area below the tables showing the estimated location or an error/warning message
