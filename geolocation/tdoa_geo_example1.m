% Passive Geolocation Homework 7 Problem 5 - Tim Cardenuto
% TDOA h(x) = r1 - r2 = ((x-a1)'*(x-a1))^.5 - ((x-a2)'*(x-a2))^.5
% H = (1/r1)*(x-a1)' - (1/r2)*(x-a2)'

clear
sigma = 50
R = sigma*sigma
xhat=[0;0]
z = [13000;-5600]
% 1 nautical mile = 1852 meters
a1 = [-10*1852;20*1852]
a2 = [0;20*1852]
a3 = [10*1852;20*1852]
count = 1
maxcount = 10
error = 100
maxerror = .1

while (count < maxcount & error > maxerror)
    h = [((xhat(:,count)-a1)'*(xhat(:,count)-a1))^.5 - ((xhat(:,count)-a2)' ...
        *(xhat(:,count)-a2))^.5; ((xhat(:,count)-a3)'*(xhat(:,count)-a3))^.5 ...
        - ((xhat(:,count)-a2)'*(xhat(:,count)-a2))^.5];
    H = [(1/((xhat(:,count)-a1)'*(xhat(:,count)-a1))^.5)*(xhat(:,count)-a1)' ...
        - (1/((xhat(:,count)-a2)'*(xhat(:,count)-a2))^.5)*(xhat(:,count)-a2)'; ...
        (1/((xhat(:,count)-a3)'*(xhat(:,count)-a3))^.5)*(xhat(:,count)-a3)' ... 
        - (1/((xhat(:,count)-a2)'*(xhat(:,count)-a2))^.5)*(xhat(:,count)-a2)'];
    P = inv(H'*inv(R)*H)
    xhat = [xhat, xhat(:,count) + P*H'*inv(R)*(z-h)]
    error = abs(xhat(:,count+1) - xhat(:,count));
    count = count+1;
end

% 95% EEP
eigenvalues = eig(P);
k = 5.9915;
semimajor = sqrt(k*max(eigenvalues))
semiminor = sqrt(k*min(eigenvalues))


plot(a1(1)/1852, a1(2)/1852, 'ob')
hold on
plot(a2(1)/1852, a2(2)/1852, 'ob')
plot(a3(1)/1852, a3(2)/1852, 'ob')
%plot(xhat(1,:)/1852, xhat(2,:)/1852, '.b')
plot(xhat(1,count)/1852,xhat(2,count)/1852, '*g')

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
plot((ellipse(:,1)+xhat(1,count))/1852,(ellipse(:,2)+xhat(2,count))/1852,'-')
axis([-15 15 0 25])
xlabel('south/north (nm)')
ylabel('west/east (nm)')
title('Geolocation Estimate')
legend('Aircraft 1', 'Aircraft 2', 'Aircraft 3', 'Estimated location', 'Error Ellipse')