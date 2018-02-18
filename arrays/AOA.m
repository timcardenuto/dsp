clear

Fc = 1e9;                          % carrier frequency
wavelength = (3*10^8 / Fc);

d = (3*10^8 / Fc)/2;     % half wavelength in meters @ center frequency Fc in Hz
N = 9;                  % number of array elements
array = [-4*d, -3*d, -2*d, -1*d, 0, d, 2*d, 3*d, 4*d]; % linear element array
figure()

for angle = 1:180
    clf
    plot(array,0,'or');
    axis([-1,1,-1,1])
    hold on 
    
    aoa = angle;                    % 45 degree incoming signal
    aoaradians = aoa*pi/180;        % in radians

    % Calculate phase delay using equation here (http://www.radartutorial.eu/06.antennas/Phased%20Array%20Antenna.en.html)
    phasedelay = 2*pi*array*sin(aoaradians)/wavelength;
    phasedelay = mod(abs(phasedelay),(2*pi)).*(phasedelay./abs(phasedelay));
    % phasedelayradians = pd*pi/180;

    % This is more likely calculation based on (https://web.wpi.edu/Pubs/E-project/Available/E-project-101012-211424/unrestricted/DirectionFindingPaper.pdf)    
    %pi_aoa = asin(wavelength*(phasedelay/(2*pi))./array);
    % Or more simply from wikipedia (same result)
    pi_aoa = asin((wavelength*phasedelay)./(2*pi*array));
    pi_aoa = pi_aoa*180/pi
    
    
    for i = 1:N
        mag = 1;
        ang = pi_aoa(i);
        x = array(i);
        y = 0;
        plot([x mag*cos(ang*pi/180)+x],[y mag*sin(ang*pi/180)+y], 'r'); 
    end
    pause
end