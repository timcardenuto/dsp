
import numpy as np
import matplotlib.pyplot as plt

num_targets = 1     # number of targets to generate
area = 10           # grid radius to search for targets
num_sensors = 1
num_mlocs = 100     # number of measurement locations (sample size)
path = np.pi           # length of flight path, in radians
err = 1             # random +- error in degrees of measurements

# determine sigma (standard deviation)
rng = np.random.default_rng(12345)
e = np.zeros((100000,1))
for i in range(100000):
    e[0]= -err + (err+err)*np.random.rand(1)
sigma = np.std(e)

plt.plot(0,0,'.b')

# generate targets within x,y coordinate grid
targets = np.zeros((num_targets,2))
for num in range(num_targets):
    x = (-area + (area+area)*np.random.rand(1))
    y = (-area + (area+area)*np.random.rand(1))
    targets[num,:] = np.hstack([x,y])
    plt.plot(x,y,'*r')


# generate measurement locations on flight path
radius = area * np.sqrt(2)
theta = np.linspace(0,path,num_mlocs)   # measurement spacing
mlocs = np.vstack([radius*np.cos(theta), radius*np.sin(theta)]).conj().transpose()
plt.plot(mlocs[:,0], mlocs[:,1],'ob')

# note - mlocs should really be a 3D matrix of [[x,y] x path x sensor]
#      - doa should be 2D matrix of [doa(s) x sensor], after association add
#       target dimension
#      - 

# flat vector b/c you don't know which measurements are which targets
doa = np.zeros((num_mlocs*num_targets,1))

# time step, update target/platform locations
for mloc in range(num_mlocs):
    # generate lobs per sensor at current mloc
    for sensor in range(num_sensors):

        for target in range(num_targets):
            print("mloc: ",mloc)
            print("sensor: ", sensor)
            print("target: ", target)   
            # calculate true angle between target and mloc 
            theta = np.arctan2(mlocs[mloc,1]-targets[target,1], mlocs[mloc,0]-targets[target,0])
            
            # plot true line to target
            #plot([targets(target,1) mlocs(mloc,1)],[targets(target,2) mlocs(mloc,2)], 'r')
            # get true length to target just for plotting to make sure we
            # can see crossover!!! this would not be known
            distance = np.sqrt(np.power(targets[target,0]-mlocs[mloc,0],2) + np.power(targets[target,1]-mlocs[mloc,1],2))
            # convert to true angle from mloc relative to due East
            theta = -(np.pi - theta)
            # add random error - not sigma, that's a statistical measure not
            # each individual measurement error
            error = -err + (err+err)*np.random.rand(1)
            # doa index to save measurement since we're using a flat vector
            index = (mloc-1)*len(targets[:,0]) + target
            doa[index] = theta + error * np.pi/180
            # don't forget to add the sensor xy offset (mlocs(mloc,1 and mlocs(mloc,2)) for the plot
            # uses the range as the magnitude for drawing a line to target
            plt.plot([mlocs[mloc,0], (distance*np.cos(doa[index])+mlocs[mloc,0])],[mlocs[mloc,1], (distance*np.sin(doa[index])+mlocs[mloc,1])], 'k')

# save doa and measurement data to cvs, last row is sigma
# measurement data x,y can be used as nautical miles
doa = np.vstack([doa, sigma])
mlocs = np.vstack([mlocs, np.hstack([sigma,sigma])]) 
data = np.hstack([doa, mlocs])             # matrix with columns [doa x y]
#csvwrite('DOAgeneration.dat',data);
#saveas(figure,'DOAgeneration.png');
plt.show()

