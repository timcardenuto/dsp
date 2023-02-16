% Passive Geolocation Homework 7 Problem 3 - Tim Cardenuto

clear
a1=[0;0];
theta1=50/180*pi;
sigma1=3*pi/180;
u1=[cos(theta1);sin(theta1)];
M1 = eye(2)-u1*u1';
plot(a1(1), a1(2), 'ok')
hold on
plot([a1(1) 100*cos(theta1)+a1(1)],[a1(2) 100*sin(theta1)+a1(2)], 'k')

a2=[20;0];
theta2=70/180*pi;
sigma2=2*pi/180;
u2=[cos(theta2);sin(theta2)];
M2 = eye(2)-u2*u2';
plot(a2(1), a2(2), 'ok')
plot([a2(1) 100*cos(theta2)+a2(1)],[a2(2) 100*sin(theta2)+a2(2)], 'k')

a3=[60;0];
theta3=110/180*pi;
sigma3=1*pi/180;
u3=[cos(theta3);sin(theta3)];
M3 = eye(2)-u3*u3';
plot(a3(1), a3(2), 'ok')
plot([a3(1) 100*cos(theta3)+a3(1)],[a3(2) 100*sin(theta3)+a3(2)], 'k')

% Initial position estimate given by Moore-Penrose solution
xhat = inv(M1+M2+M3)*(M1*a1+M2*a2+M3*a3)
plot(xhat(1), xhat(2), '*r')

% Final Estimation by Iterated Least Squares 
% h(x) = theta = arctan((x2-a2)/(x1-a1))
% H = (1/r)*[-sin(h) cos(h)]
count = 1;
maxcount = 10;
error = 100;
maxerror = .1;
xhatold = xhat;
theta = [theta1; theta2; theta3]
R = [sigma1*sigma1, 0, 0; 0, sigma2*sigma2, 0; 0, 0, sigma3*sigma3]
while (count < maxcount & error > maxerror)
    %thetahat = atan2((xhatold(2)-a3(2)),(xhatold(1)-a3(1)))
    thetahat = [atan2((xhatold(2)-a1(2)),(xhatold(1)-a1(1))); ...
        atan2((xhatold(2)-a2(2)),(xhatold(1)-a2(1))); ...
        atan2((xhatold(2)-a3(2)),(xhatold(1)-a3(1)))];
    %H = (1/(((xhatold-a3)'*(xhatold-a3))^.5))*[-sin(thetahat) cos(thetahat)]
    H = [(1/(((xhatold-a1)'*(xhatold-a1))^.5))*[-sin(thetahat(1)) cos(thetahat(1))]; ...
        (1/(((xhatold-a2)'*(xhatold-a2))^.5))*[-sin(thetahat(2)) cos(thetahat(2))]; ...
        (1/(((xhatold-a3)'*(xhatold-a3))^.5))*[-sin(thetahat(3)) cos(thetahat(3))]];
    P = inv(H'*inv(R)*H);
    xhatnew = xhatold + P*H'*inv(R)*(theta-thetahat)
    error = abs(xhatnew - xhatold)
    xhatold = xhatnew;
    count = count+1;
end

% 95% EEP
eigenvalues = eig(P);
k = 5.9915;
semimajor = sqrt(k*max(eigenvalues))
semiminor = sqrt(k*min(eigenvalues))

plot(xhatnew(1), xhatnew(2), '*g')
hold on
% Create and scale ellipse
theta = linspace(0,2*pi);
ellipse = [semimajor*cos(theta);semiminor*sin(theta)]';
% Rotate ellipse
[U, D] = eig(P);
eigenvector = U(:,2);
angle = atan2(eigenvector(2), eigenvector(1));
if(angle < 0)
    angle = angle + 2*pi;
end
R = [cos(angle),sin(angle);-sin(angle),cos(angle)];
ellipse = ellipse * R;
% Shift ellipse and plot 
plot((ellipse(:,1)+xhatnew(1)),(ellipse(:,2)+xhatnew(2)),'-')
% axis([-1 70,-5 1])
