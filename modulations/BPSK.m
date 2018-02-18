clc
clear

%% Create test signal
% Generation of bit pattern
s=round(rand(1,25));    % Generating 20 bits
signal=[];  
carrier=[];
t=[0:2*pi/119:2*pi];     % Creating 60 samples for one cosine 
for k=1:25
    if s(1,k)==0
        sig=-ones(1,120);    % 120 minus ones for bit 0
    else
        sig=ones(1,120);     % 120 ones for bit 1
    end
    c=cos(t);   
    carrier=[carrier c];
    signal=[signal sig];
end
subplot(4,1,1);
plot(signal);
axis([-100 3100 -1.5 1.5]);
title('\bf\it Original Bit Sequence');

% BPSK Modulation of the signal
bpsk_sig=signal.*carrier;   % Modulating the signal
subplot(4,1,2);
plot(bpsk_sig)
axis([-100 3100 -1.5 1.5]);
title('\bf\it BPSK Modulated Signal');

%% 


% Frequency data
n = 2^nextpow2(length(bpsk_sig)); % find power of two larger than data length to optimize fft
y = fft(bpsk_sig,n);
f = 1000/2*linspace(0,1,n/2+1);
figure()
plot(f,2*abs(y(1:n/2+1)/2))



moment(signal,4)
var(signal)
kurtosis(signal)
skewness(signal)
