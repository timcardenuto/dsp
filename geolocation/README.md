# Geolocation
This software was originally made based on a class I took, and added to over time for research purposes.
I've started using AI to assist me with things I've wanted to do and haven't had time for, like creating an interactive tool for simulating geolocation measurements.

## TODO
* Darkmode support, there's a dark map but colors of sensors/targets/measurements need to change. Also should just have a full dark mode for the rest of the page as well.
* Want the DOA measurement shown on the angles drawn. Want them relative to True North, I think they all might be relative to positive x axis (East). We want to change the math for the carteasian geo algorithm to use North as the zero point if possible.... 
* Probably want a dropdown menu for setting lat/lon and also other settings
* Want a mode for uploading csv or entering a table of measurements for sensor(s) position and measurements and want to plot them all then calculate convergence points WITHOUT a known target location (so basically NOT a sim but when you have the normal set of data)
* Need to add 3D coordinates and elevation data?


## MATLAB to Python Cheatsheet
Coming from MATLAB? Here are some differences in the Python code (mainly about the Numpy library). A bigger list is [here](https://numpy.org/doc/stable/user/numpy-for-matlab-users.html). Also note that I use the Numpy array data type. There *is* a matrix data type in Numpy which has different syntax but they plan to deprecate it so I don't use it. 


| What I want to do                  | MATLAB     | Python/Numpy |
|------------------------------------|------------|---------------------|
| 1 x 3 matrix (1 row, 3 columns)    | [0,1,2]    | numpy.hstack([0,1,2]) |
| 3 x 1 matrix (3 rows, 1 column)    | [0;1;2]    | numpy.vstack([0,1,2]) |
| Matrix multiplication              | a*b        | a@b |
| Matrix element-wise multiplication | a.*b       | a*b |
| Matrix transpose                   | a'         | a.conj().transpose() |
| Matrix inverse                     | inv(a)     | numpy.linalg.inv(a) |
| power                              | a^b        | numpy.power(a,b) |
| arc tangent                        | atan2(a,b) | numpy.arctan2(a,b) |
| 2x2 identity matrix                | eye(2)     | numpy.identity(2) |
