clf
clear

Fm = 10;                             % message frequency
Fc = Fm*10;                          % carrier frequency
T = 1/Fc;                            % carrier period
sample_rate = Fc*2;                  % # of samples per carrier period (T)
t = 0:T/sample_rate:(30*T - T/1000); % same # of samples no matter what Fc
Am = 5;                              % message amplitude
Sm = Am*sin(2*pi*Fm*t);     % message signal

m = 1;                               % modulation index m = Am/Ac
Ac = Am/m;                           % carrier amplitude based on modulation index
Sc = Ac*sin(2*pi*Fc*t);     % carrier signal

amSig = Sc.*(1 + m/Am*Sm);    % modulated signal

rectifiedSig = abs(amSig);      % rectified signal

% low pass filter signal
averaging_length = Fc;    % to get envelope recovered, use average = # of samples per Fc period = sample_rate
filtered_signal = zeros(1,averaging_length);
for i=averaging_length+1:length(rectifiedSig)
    filtered_signal = [filtered_signal sum(rectifiedSig((i-averaging_length):i-1))/averaging_length];
end

%% Plot everything

subplot(5,1,1)
plot(t,Sm)
title ([num2str(Fm),' Hz Analog Message']);
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,2)
plot(t,Sc)
title ([num2str(Fc),' Hz Carrier']);
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,3)
plot(t,amSig)
title ('Amplitude Modulated (AM) Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,4)
plot(t,rectifiedSig)
title ('Rectified Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(5,1,5)
plot(t,filtered_signal)
title ([num2str(averaging_length),' Sample Average Filtered Signal']);
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');