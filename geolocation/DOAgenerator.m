clear
hold on
axis([-20,20,-20,20])

num_targets = 1     % number of targets to generate
area = 10           % grid radius to search for targets
num_sensors = 1
num_mlocs = 100     % number of measurement locations (sample size)
path = pi           % length of flight path, in radians
err = 1             % random +- error in degrees of measurements

% determine sigma (standard deviation)
rng('shuffle', 'twister')
e = zeros(100000,1);
for i = 1:100000 
    e(i)= -err + (err+err)*rand(1);
end
sigma = std(e);

figure = plot(0,0,'.b')

% generate targets within x,y coordinate grid
targets = zeros(num_targets,2);
for num = 1:num_targets
    x = (-area + (area+area)*rand(1));
    y = (-area + (area+area)*rand(1));
    targets(num,:) = [x,y];
    plot(x,y,'*r')
end

% generate measurement locations on flight path
radius = area * sqrt(2);
theta = linspace(0,path,num_mlocs);   % measurement spacing
mlocs = [radius*cos(theta);radius*sin(theta)]';
plot(mlocs(:,1),mlocs(:,2),'ob');

% note - mlocs should really be a 3D matrix of [[x,y] x path x sensor]
%      - doa should be 2D matrix of [doa(s) x sensor], after association add
%       target dimension
%      - 

% flat vector b/c you don't know which measurements are which targets
doa = zeros(num_mlocs*num_targets,1); 

% time step, update target/platform locations
for mloc = 1:num_mlocs
    % generate lobs per sensor at current mloc
    for sensor = 1:num_sensors
        for target = 1:num_targets
           % calculate true angle between target and mloc 
           theta = atan2(mlocs(mloc,2)-targets(target,2),mlocs(mloc,1)-targets(target,1));
           % plot true line to target
           %plot([targets(target,1) mlocs(mloc,1)],[targets(target,2) mlocs(mloc,2)], 'r')
           % get true length to target just for plotting to make sure we
           % can see crossover!!! this would not be known
           range = sqrt((targets(target,1)-mlocs(mloc,1))^2 + (targets(target,2)-mlocs(mloc,2))^2);
           % convert to true angle from mloc relative to due East
           theta = -(pi - theta);
           % add random error - not sigma, that's a statistical measure not
           % each individual measurement error
           error = -err + (err+err)*rand(1);
           % doa index to save measurement since we're using a flat vector
           index = (mloc-1)*length(targets(:,1)) + target; 
           doa(index) = theta + error * pi/180;
           % don't forget to add the sensor xy offset (mlocs(mloc,1 and mlocs(mloc,2)) for the plot
           % uses the range as the magnitude for drawing a line to target
           plot([mlocs(mloc,1) (range*cos(doa(index))+mlocs(mloc,1))],[mlocs(mloc,2) (range*sin(doa(index))+mlocs(mloc,2))], 'k')
       end
    end 
end

% save doa and measurement data to cvs, last row is sigma
% measurement data x,y can be used as nautical miles
doa = [doa; sigma];
mlocs = [mlocs; [sigma,sigma]]; 
data = [doa mlocs];             % matrix with columns [doa x y]
csvwrite('DOAgeneration.dat',data);
saveas(figure,'DOAgeneration.png');

