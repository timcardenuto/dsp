
import numpy as np
import matplotlib.pyplot as plt

plt.ion()
fig,axs = plt.subplots(1,1)
axs.set_xlabel('west/east (nm)')
axs.set_ylabel('south/north (nm)')
axs.set_title('DOA Scenario Generation')

num_targets = 5     # number of targets to generate
area = 10           # grid radius to place targets
num_sensors = 1
num_mlocs = 100     # number of measurement locations (sample size)
path = 0.5*np.pi    # length of flight path, in radians
err = 1             # max measurement error in degrees, +-

# TODO make seed a configuration, turn on/off
np.random.seed(12345)

# determine sigma (standard deviation)
e = np.zeros((100000,1))
for i in range(100000):
    e[0] = -err + (err+err)*np.random.rand(1)
sigma = np.std(e)
print("sigma: ",sigma)

# generate targets within x,y coordinate grid
targets = np.zeros((num_targets,2))
for num in range(num_targets):
    # x = (-area + (area+area)*np.random.rand(1))  # not sure what shape I was going for here
    # y = (-area + (area+area)*np.random.rand(1))
    x = np.random.uniform(-area,area)  # using a box instead
    y = np.random.uniform(-area,area)
    targets[num,:] = np.hstack([x,y])
    print("Target",num,": ",x,",",y)
    plt.plot(x,y,'*r')


# generate measurement locations on flight path
radius = area * np.sqrt(2)
theta = np.linspace(0,path,num_mlocs)   # measurement spacing
mlocs = np.vstack([radius*np.cos(theta), radius*np.sin(theta)]).conj().transpose() # measurement locations x,y (mlocs)
plt.plot(mlocs[:,0], mlocs[:,1],'ob')

# note - mlocs should really be a 3D matrix of [[x,y] x path x sensor]
#      - doa should be 2D matrix of [doa(s) x sensor], after association add
#       target dimension
#      - 

# flat vector b/c you don't know which measurements are which targets
doa = np.zeros((num_mlocs*num_targets,1))

# time step, update target/platform locations
for timestep in range(num_mlocs):
    # generate doa per sensor at current mloc
    for sensor in range(num_sensors):

        for target in range(num_targets):
            #print("timestep: ",timestep)
            #print("sensor: ", sensor)
            #print("target: ", target)   

            # Calculate relative angle between target and measurement location 
            theta = np.arctan2(mlocs[timestep,1]-targets[target,1], mlocs[timestep,0]-targets[target,0])
            
            # Calculate true angle to target, from measurement location relative to due East
            theta = -(np.pi - theta)

            # Add measurement error based on sensor sigma
            error = -err + (err+err)*np.random.rand(1)
            theta = theta + error * np.pi/180
            # print("error: ",error)
            # print("theta: ",theta)

            # plot true line to target
            # plt.plot([targets[target,0], mlocs[timestep,0]],[targets[target,1], mlocs[timestep,1]], 'r')

            # Get distance to target just for plotting to make sure we can see crossover. This would not be known in real life. Increase the length by 1.5 so it's more realistic.
            distance = 1.5 * np.sqrt(np.power(targets[target,0]-mlocs[timestep,0],2) + np.power(targets[target,1]-mlocs[timestep,1],2))
            # print(distance)

            # Save the measurement to the doa array 
            # TODO why use flat array and not 2 dimensions? seems like that would be more straightforward
            index = (timestep-1)*len(targets[:,0]) + target
            doa[index] = theta
            # print(index)
            # print(doa[index])

            # Draw the measured DOA (with error) from measurement location to target, using the distance as the magnitude
            plt.plot([mlocs[timestep,0], (distance*np.cos(doa[index][0])+mlocs[timestep,0])],[mlocs[timestep,1], (distance*np.sin(doa[index][0])+mlocs[timestep,1])], 'k')


# print("doa: ",doa)
# print("sigma: ",sigma)
# print("mlocs: ",mlocs)
# save doa and measurement data to csv, last row is sigma
# measurement data x,y can be used as nautical miles
# doa = np.vstack([doa, sigma])
# mlocs = np.vstack([mlocs, np.hstack([sigma,sigma])]) 
# data = np.hstack([doa, mlocs])             # matrix with columns [doa x y]
#csvwrite('DOAgeneration.dat',data);
#saveas(figure,'DOAgeneration.png');


plt.draw()
input("Press Enter to continue...")

