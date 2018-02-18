clear

%% PM
Fc = 1000;
Fm = 100;
T = 1/Fc;
t = 0:T/999:30*T;
carrier = cos(2*pi*Fc*t);
message = cos(2*pi*Fm*t);
kp = 2*pi;
pm = cos(2*pi*Fc*t + kp*message);
subplot(5,1,1)
plot(t,carrier)
title ('Carrier');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,2)
plot(t,message)
title ('PM Message');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,3)
plot(t,pm)
title ('Phase Modulation');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');

%% PSK
message = square(2*pi*Fm*t);
psk = cos(2*pi*Fc*t).*message;

subplot(5,1,4)
plot(t,message)
title ('PSK Message');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,5)
plot(t,psk)
title ('Phase Shift Keying');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');