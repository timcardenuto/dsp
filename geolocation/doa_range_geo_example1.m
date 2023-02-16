% Passive Geolocation Homework 8 Problem 2 - Tim Cardenuto

clear
% True Target location
target = plot(0,0,'or')
hold on

% Buoy 1 - DOA
a1=[(sin(60/180*pi)*15*1852);(cos(60/180*pi)*15*1852)];
range1=15*1852;
theta1=-(90+60)/180*pi;
sigma1=2/180*pi;
u1=[cos(theta1);sin(theta1)];
M1 = eye(2)-u1*u1';
buoys = plot(a1(1)/1852, a1(2)/1852, 'ok')
plot([a1(1)/1852 (100*cos(theta1)/1852)],[a1(2)/1852 (100*sin(theta1)/1852)], 'k')

% Buoy 2 - range
a2=[-(sin(60/180*pi)*10*1852);(cos(60/180*pi)*10*1852)];
range2=10*1852;
theta2=-(90-60)/180*pi;
sigma2=100;
u2=[cos(theta2);sin(theta2)];
M2 = eye(2)-u2*u2';
plot(a2(1)/1852, a2(2)/1852, 'ok')
th = 0:pi/50:2*pi;
xunit = range2/1852 * cos(th) + a2(1)/1852;
yunit = range2/1852 * sin(th) + a2(2)/1852;
plot(xunit, yunit, 'k');

% Initial position estimate given by Moore-Penrose solution
xhat = inv(M1+M2)*(M1*a1+M2*a2);
guess = plot(xhat(1)/1852, xhat(2)/1852, '*k')

% Final Estimation by Iterated Least Squares 
% DOA h and H
% h(x) = theta = arctan((x2-a2),(x1-a1))
% H = (1/r)*[-sin(h) cos(h)]
% Range h and H
% h = r = ((x-a)'*(x-a))^.5
% H = (1/r)*(x-a)'
count = 1;
maxcount = 10;
error = 100;
maxerror = .1;
xhatold = xhat;
% Using two measurements - 1 DOA and 1 range
z = [theta1;range2];
R = [sigma1*sigma1, 0; 0, sigma2*sigma2];
while (count < maxcount & error > maxerror)
    h = [atan2((xhatold(2)-a1(2)),(xhatold(1)-a1(1))); ...
        ((xhatold-a2)'*(xhatold-a2))^.5];
    H = [(1/(((xhatold-a1)'*(xhatold-a1))^.5))*[-sin(h(1)) cos(h(1))]; ...
        (1/(((xhatold-a2)'*(xhatold-a2))^.5))*(xhatold-a2)'];
    P = inv(H'*inv(R)*H);
    xhatnew = xhatold + P*H'*inv(R)*(z-h);
    error = abs(xhatnew - xhatold);
    xhatold = xhatnew;
    count = count+1
end

% 95% EEP
eigenvalues = eig(P);
k = 5.9915;
semimajor = sqrt(k*max(eigenvalues))
semiminor = sqrt(k*min(eigenvalues))

final = plot(xhatnew(1)/1852, xhatnew(2)/1852, '*b')
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
eep = plot((ellipse(:,1)+xhatnew(1))/1852,(ellipse(:,2)+xhatnew(2))/1852,'-b')
xlabel('west/east (nm)')
ylabel('south/north (nm)')
title('Geolocation Estimate')
legend([target,buoys,guess,final,eep],'True Target','Buoys','Initial Estimate','Final Estimate','95% Ellipse')