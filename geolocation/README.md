

## MATLAB to Python Cheatsheet
Coming from MATLAB? Here are some differences in the Python code (mainly about the Numpy library). A bigger list is [here](https://numpy.org/doc/stable/user/numpy-for-matlab-users.html). Also note that I use the Numpy array data type. There *is* a matrix data type in Numpy which has different syntax but they plan to deprecate it so I don't use it. 


|--- What I want to do ---|--- MATLAB syntax ---|--- Python Numpy syntax ---|
| 1 x 3 matrix (1 row, 3 columns) | [0,1,2] | numpy.hstack([0,1,2]) |
| 3 x 1 matrix (3 rows, 1 column) | [0;1;2] | numpy.vstack([0,1,2]) |
| Matrix multiplication | a*b | a@b |
| Matrix element-wise multiplication | a.*b | a*b |
| Matrix transpose | a' | a.conj().transpose() |
| Matrix inverse | inv(a) | numpy.linalg.inv(a) |
| power | a^b | numpy.power(a,b) |
| arc tangent | atan2(a,b) | numpy.arctan2(a,b) |
| 2x2 identity matrix | eye(2) | numpy.identity(2) |
