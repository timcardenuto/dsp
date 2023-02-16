# Linear Estimation Problem 5 - Fit the following model to measured data z
# t = [0:.01:3]'
# v = normpdf(t,0,(sigma*sigma))
# z = x1 + x2*exp(-t)*sin(4*pi*t) + v
 
 
import numpy as np
import matplotlib.pyplot as plt

def normpdf(x, mu, sigma):
    u = (x-mu)/abs(sigma)
    y = (1/(np.sqrt(2*np.pi)*abs(sigma)))*np.exp(-u*u/2)
    return y


# x1 and x2 estimate
x1 = 1
x2 = 1
sigma = 0.2
t = np.transpose(np.arange(0,3,0.01))
v = normpdf(t, 0, (sigma*sigma))
z = x1 + x2*np.exp(-t)*np.sin(4*np.pi*t) + v
H = np.vstack((np.ones(len(t)), np.exp(-t)*np.sin(4*np.pi*t)))
H = np.transpose(H)
R = (sigma*sigma)*np.identity(len(t))

try:
    # xhat = inv(H' * inv(R) * H) * H' * inv(R) * z
    xhat = np.linalg.inv((np.transpose(H).dot(np.linalg.inv(R))).dot(H)).dot(np.transpose(H)).dot(np.linalg.inv(R)).dot(z)
except np.linalg.LinAlgError:
    print("Not invertible. Skip this one.")
    pass
else:
    print("continue with what you were doing")

print("xhat = ",xhat)

# Plot measured data and estimates
y = xhat[0] + np.multiply(xhat[1]*np.exp(-t),np.sin(4*np.pi*t))
plt.subplot(2,2,1)
plt.plot(t, z, '.', label='data z')
plt.plot(t, y, label='estimated xhat')
plt.xlabel('Time t')
plt.ylabel('Data z or xhat')
plt.title('Measured data z and fitted curve')

 
# Estimate error covariance matrix P
P = np.linalg.inv(np.transpose(H).dot(np.linalg.inv(R)).dot(H))
print("Covariance Matrix = ",P)

# 95% confidence intervals
x1min = xhat[0]-np.sqrt(3.8415*P[0][0])
x1max = xhat[0]+np.sqrt(3.8415*P[0][0])
x2min = xhat[1]-np.sqrt(3.8415*P[1][1])
x2max = xhat[1]+np.sqrt(3.8415*P[1][1])
 
# 95% EEP
w,v = np.linalg.eig(P)
print("Eigenvalues = ",w)
k = 5.9915
semimajor = np.sqrt(k*max(w))
semiminor = np.sqrt(k*min(w))
 
## Plot 95% confidence intervals and EEP
 
# Plot estimate location & confidence intervals
xview=np.array([[1.6,2.6],[1.6,2.6]])
yview=np.array([[2.5,3.5],[2.5,3.5]])
x1axis=np.array([[x1min,x1max], [x1min,x1max]])
x2axis=np.array([[x2min,x2max], [x2min,x2max]])

plt.subplot(2,2,2)
plt.plot(xhat[0], xhat[1],'*', label='xhat')
plt.plot(x1axis, xview, label='x1')
plt.plot(yview, x2axis, label='x2')
 
# Create and scale ellipse
theta = np.linspace(0,2*np.pi)
ellipse = np.transpose(np.vstack([semimajor*np.cos(theta), semiminor*np.sin(theta)]))
#plot(ellipse(:,1),ellipse(:,2),'-')
 
# Rotate ellipse
eigenvector = v[:,1]
print("Eigenvector = ",eigenvector)
angle = np.arctan2(eigenvector[1], eigenvector[0])
if(angle < 0):
    angle = angle + 2*np.pi

R = np.vstack([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])
ellipse = ellipse.dot(R)
#plot(ellipse(:,1),ellipse(:,2),'-r')
 
# Shift ellipse and plot 
plt.plot(ellipse[:,0]+xhat[0], ellipse[:,1]+xhat[1], '-', label='Error Ellipse')
plt.xlabel('x1')
plt.ylabel('x2')
plt.title('95% Confidence Intervals and EEP')
plt.show()


