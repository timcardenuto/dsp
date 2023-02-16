% Passive Geolocation Homework9 Problem 1 - Tim Cardenuto
clear
m = matfile('data_association_hw')
xhats = [m.xhat1, m.xhat2, m.xhat3];
Ps = [m.P1, m.P2, m.P3];
ac = m.ac;
z = m.z/180*pi;
sigma = m.sigma/180*pi;

B=.05;
a=.95;
k=3.8415;
largestLL=0;
i=1;
j=1;
for (count=1:length(xhats))
    xhat = xhats([i;i+1]);
    P=Ps([j,j+2;j+1,j+3]);
    zhat = atan2(xhat(2)-ac(2),xhat(1)-ac(1));
    H = (1/(((xhat-ac)'*(xhat-ac))^.5))*[-sin(zhat) cos(zhat)];
    md = ((z-zhat)^2) / (H*P*H'+(sigma^2));
    if md < k
        LL = -.5*(log(2*pi)+log(H*P*H'+sigma^2)+md);
        if LL > largestLL
            largestLL = LL
            targetxhat = xhat
            targetP = P
            targetzhat = zhat;
            targetH = H;
        end
    end
    i=i+2;
    j=j+4;
end

% Plot measurement and most likely target 
sensor = plot(ac(1), ac(2), 'ok')
hold on
measurement = plot([ac(1) 10*cos(z)],[ac(2) 10*sin(z)], 'k')
error = plot([ac(1) 10*cos(z+sigma)],[ac(2) 10*sin(z+sigma)], '--k')
plot([ac(1) 10*cos(z-sigma)],[ac(2) 10*sin(z-sigma)], '--k')
targetellipse = drawEllipse(targetxhat, targetP)
xlabel('west/east (nm)')
ylabel('south/north (nm)')
title('DOA Target Association')
legend([sensor,measurement,error,targetellipse],'Sensor Location', ...
    'DOA Measurement', 'DOA Sigma', 'Associated Target')

% Add new measurement to most likely target
% xhat = targetxhat;
% h = targetzhat;
% H = targetH;
% P = targetP;
% R = sigma^2;    
% K = P*H'*inv(H*P*H'+R);
% xhat = xhat+K*(z-h);
% P = (1-K*H)*P;
% drawEllipse(xhat, P)
