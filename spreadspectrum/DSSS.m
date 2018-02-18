% Lab 06
% WiCom_3
% By Kashif Shahzad 
% 01-ET-31
% 3rd July 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct Sequence Spread Spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

% Generating the bit pattern with each bit 6 samples long
b=round(rand(1,20));
pattern=[];
for k=1:20
    if b(1,k)==0
        sig=zeros(1,6);
    else
        sig=ones(1,6);
    end
    pattern=[pattern sig];
end
plot(pattern);
axis([-1 130 -.5 1.5]);
title('\bf\it Original Bit Sequence');

% Generating the pseudo random bit pattern for spreading
spread_sig=round(rand(1,120));
figure,plot(spread_sig);
axis([-1 130 -.5 1.5]);
title('\bf\it Pseudorandom Bit Sequence');

% XORing the pattern with the spread signal
hopped_sig=xor(pattern,spread_sig);

% Modulating the hopped signal
dsss_sig=[];
t=[0:100];    
fc=.1
c1=cos(2*pi*fc*t);
c2=cos(2*pi*fc*t+pi);
for k=1:120
    if hopped_sig(1,k)==0
        dsss_sig=[dsss_sig c1];
    else
        dsss_sig=[dsss_sig c2];
    end
end
figure,plot([1:12120],dsss_sig);
axis([-1 12220 -1.5 1.5]);
title('\bf\it DSSS Signal');

% Plotting the FFT of DSSS signal
figure,plot([1:12120],abs(fft(dsss_sig)))