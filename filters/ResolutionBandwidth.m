clear
clf

Fm = 100;                             % message frequency
Fc = Fm*10;                          % carrier frequency
Tc = 1/Fc;                            % carrier period
sample_rate = Fc*2;                  % # of samples per carrier period (T)
t = 0:Tc/2:(400*Tc/2);

Am = 5;                              % message amplitude
Sm = Am*sin(2*pi*Fm*t);     % message signal

m = 1;                               % modulation index m = Am/Ac
Ac = Am/m;                           % carrier amplitude based on modulation index
Sc = Ac*sin(2*pi*Fc*t);     % carrier signal

% t = 1:1:20000;
% 
% Fc = 100;                          % carrier frequency
% Sc = sin(2*pi*t/Fc);
% Fm = 1000;                        % message frequency
% Sm = sin(2*pi*t/Fm);

multSig = Sm .* Sc;         % composite signal, not really AM

% check spectrum
Fs = 1000;            % Sampling frequency                    
Tc = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*Tc;        % Time vector
S = sin(2*pi*50*t) .* sin(2*pi*120*t);
Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fc*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

return
%% Narrowband Filter
interval = 900;     % sample average interval = filter resonant frequency
num_samples = 500;    % num samples to use in average = inverse filter bandwidth (sensitivity)
sig_in = multSig;
sig_out = [];
for i=1:length(t)

    sumSig = sig_in(i);
    for j=1:(num_samples-1)
        % what do you do when run out of samples in signal for the averaging?
        % use zero? stop running experiment cuz this wouldn't happen
        % IRL? % do nothing?
        if (i+(j*interval)) > length(t)   
            sumSig = sumSig + 0;
        else  % average normally
            sumSig = sumSig + sig_in(i+(j*interval));
        end
    end
    avrgSig = sumSig / num_samples;

    sig_out = [sig_out avrgSig];
end
nbSig = sig_out;


%% Wideband Filter
interval2 = 950;      % sample average interval = filter resonant frequency
num_samples2 = 75;    % num samples to use in average = inverse filter bandwidth (sensitivity)
sig_in = multSig;
sig_out = [];
for i=1:length(t)

    sumSig = sig_in(i);
    for j=1:(num_samples2-1)
        % what do you do when run out of samples in signal for the averaging?
        % use zero? stop running experiment cuz this wouldn't happen
        % IRL? % do nothing?
        if (i+(j*interval2)) > length(t)   
            sumSig = sumSig + 0;
        else  % average normally
            sumSig = sumSig + sig_in(i+(j*interval2));
        end
    end
    avrgSig = sumSig / num_samples2;

    sig_out = [sig_out avrgSig];
end
wbSig = sig_out;

%% Plot everything
% limit = T;
subplot(5,1,1)
plot(t,Sm)
title ([num2str(Fm),' Hz Analog Message']);
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
% xlim([1 limit]);
subplot(5,1,2)
plot(t,Sc)
title ([num2str(Fc),' Hz Carrier']);
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
% xlim([1 limit]);
subplot(5,1,3)
plot(t,multSig)
title ('Composite Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
% xlim([1 limit]);
subplot(5,1,4)
plot(t,nbSig)
title ([num2str(num_samples),' Sample Narrowband Filtered Signal']);
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
% xlim([1 limit]);
subplot(5,1,5)
plot(t,wbSig)
title ([num2str(num_samples2),' Sample Wideband Filtered Signal']);
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
% xlim([1 limit]);