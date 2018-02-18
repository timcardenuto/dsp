%% FM

% Sample time scale
Fs = 2000; % should be > 2 times the highest frequency component
T = 1/Fc;
t = 0:T/999:30*T;

% Carrier
Ac = 5;
Fc = 1000;
carrier = Ac*cos(2*pi*Fc*t);

% Message
Fm = 100;
Am = 5;
message = Am*cos(2*pi*Fm*t);

kf = 100;
fm = Ac*cos(2*pi*Fc*t + (kf/Fm)*Am*sin(2*pi*Fm*t));

subplot(5,1,1)
plot(t,carrier)
title ('Carrier');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,2)
plot(t,message)
title ('FM Message');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,3)
plot(t,fm)
title ('Frequency Modulation');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');


%% FSK
% message = square(2*pi*Fm*t);
% fsk = cos(2*pi*Fc*t).*cumsum(message);
% subplot(5,1,4)
% plot(t,message)
% title ('FSK Message');
% xlabel ('time(sec)');
% ylabel ('Amplitude(volt)');
% subplot(5,1,5)
% plot(t,fsk)
% title ('Frequency Shift Keying');
% xlabel ('time(sec)');
% ylabel ('Amplitude(volt)');

N = 8;

% Generate a random bit stream
bit_stream = [1 0 1 0 1 0 1 0]

% Enter the two frequencies 
% Frequency component for 0 bit
f1 = 5; 

% Frequency component for 1 bit
f2 = 10;

% Sampling rate - This will define the resoultion
fs = 2000;

% Time for one bit
t = 0: 1/fs : 1;

% This time variable is just for plot
time = [];

FSK_signal = [];
Digital_signal = [];

for ii = 1: 1: length(bit_stream)
    
    % The FSK Signal
    FSK_signal = [FSK_signal (bit_stream(ii)==0)*sin(2*pi*f1*t)+...
        (bit_stream(ii)==1)*sin(2*pi*f2*t)];
    
    % The Original Digital Signal
    Digital_signal = [Digital_signal (bit_stream(ii)==0)*...
        zeros(1,length(t)) + (bit_stream(ii)==1)*ones(1,length(t))];
    
    time = [time t];
    t =  t + 1;
   
end

% Plot the FSK Signal
subplot(5,1,4);
plot(time,Digital_signal,'r','LineWidth',2);
title ('FSK Message');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
% axis([0 time(end) -0.5 1.5]);
grid  on;

% Plot the Original Digital Signal
subplot(5,1,5);
plot(time,FSK_signal);
title ('Frequency Shift Keying');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
% axis([0 time(end) -1.5 1.5]);
grid on;