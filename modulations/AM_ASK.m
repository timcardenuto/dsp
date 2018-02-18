clf
clear

Fm = 10;
Fc = Fm*10;
T = 1/Fc;
t = 0:T/999:30*T;
Am = 5;
Sm = Am*sin(2*pi*Fm*t);

m = 1;      % modulation index m = Am/Ac
Ac = Am/m;  % carrier amplitude based on modulation index
Sc = Ac*sin(2*pi*Fc*t);

subplot(5,1,1)
plot(t,Sc)
title ('Carrier');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,2)
plot(t,Sm)
title ('Analog Message');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');

subplot(5,1,3)
am1 = Sc.*(1 + m/Am*Sm);
plot(t,am1)
title ('Amplitude Modulated (AM) Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
% return

%% ASK
message = Am*square(2*pi*Fm*t) + Am;
ask = sin(2*pi*Fc*t).*message;

subplot(5,1,4)
plot(t,message)
title ('Digital Message');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,5)
plot(t,ask)
title ('Amplitude Shift Keying (ASK) Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');

