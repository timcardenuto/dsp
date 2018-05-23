% Example of GPS signal transmit/receive
% Using Space Vehicle (SV) signal, Psuedo-Random Number (PRN) 1 = taps [2, 6]
clear

%% Transmitter portion

% generate the chipping, aka spreading, sequence
% PRN = generateGPSPRN([2,6])

% generate the chipping, aka spreading, sequence
G1 = generateGenericPRN(10, ones(1,10), [3,10], [10]);
G2_1 = generateGenericPRN(10, ones(1,10), [2,3,6,8,9,10], [2,6]);
PRN1 = bitxor(G1,G2_1);

G2_5 = generateGenericPRN(10, ones(1,10), [2,3,6,8,9,10], [1,9]);
PRN5 = bitxor(G1,G2_5);

G2_16 = generateGenericPRN(10, ones(1,10), [2,3,6,8,9,10], [9,10]);
PRN16 = bitxor(G1,G2_16);

% map unipolar to bipolar (0 => -1)
for i=1:length(PRN1)
    if (PRN1(i) == 0)
        PRN1(i) = -1;
    end
    if (PRN5(i) == 0)
        PRN5(i) = -1;
    end 
    if (PRN16(i) == 0)
        PRN16(i) = -1;
    end
end

% baseband messages
Fs = 10;
T = 1/Fs;
t = 0:T/31:((100*T)-(T/31)); % to get longer signal, more samples, just add more copies (increase 10*T)
message = sin(2*pi*t);

% actually instead of expanding single chipping sequence to be = to
% 1 message period, repeat the chipping sequence over the message length as
% much as needed
sampled_PRN1 = [];
sampled_PRN5 = [];
sampled_PRN16 = [];
while length(sampled_PRN1) < length(message)
    sampled_PRN1 = [sampled_PRN1 PRN1];
    sampled_PRN5 = [sampled_PRN5 PRN5];
    sampled_PRN16 = [sampled_PRN16 PRN16];
end

% discard any runover bits.... shouldn't really matter
if length(sampled_PRN1) > length(message)
    sampled_PRN1 = sampled_PRN1(1:length(message));
    sampled_PRN5 = sampled_PRN5(1:length(message));
    sampled_PRN16 = sampled_PRN16(1:length(message));
end

% spread message with chipping sequence
code1 = sampled_PRN1 .* message;
code5 = sampled_PRN5 .* message;
code16 = sampled_PRN16 .* message;

% TODO combine with modulation
transmit_signal = code1 + code5 + code16;

%% Receiver portion

% despread received signal with expected chipping sequence
% TODO add noise
receive_signal = transmit_signal .* sampled_PRN1;

% low pass filter despread signal
averaging_length = 31;
filtered_signal = zeros(1,averaging_length);
for i=averaging_length+1:length(receive_signal)
    filtered_signal = [filtered_signal sum(receive_signal((i-averaging_length):i-1))/averaging_length];
end

%% Plots

% transmissions
subplot(4,1,1)
plot(t,message)
title ('Message');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,1,2)
plot(t,sampled_PRN1)
title ('PRN 1 Chipping Sequence');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,1,3)
plot(t,code1)
title ('Spread PRN 1 Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,1,4)
plot(t,transmit_signal)
title ('Combined PRN 1,5,16 RF Transmitted Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');

figure()
subplot(3,1,1)
plot(t,transmit_signal)
title ('Combined PRN 1,5,16 RF Received Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(3,1,2)
plot(t,receive_signal)
title ('Despread PRN 1 Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(3,1,3)
plot(t,filtered_signal)
title ('Filtered/Recovered Message from PRN 1');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
hold on
plot(t,message)

