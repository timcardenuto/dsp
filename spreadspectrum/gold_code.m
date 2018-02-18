%% 5 bit Gold code example
clear

L = 5;          % # of bits, or shift registers
N = 2^L - 1;    % # of codes, also sequence repeat length in bits

% initial 5 register values (the 'seed')
A1 = [0 0 0 0 1];
A2 = [0 0 0 0 1];
B1 = [0 0 0 1 0];
B2 = [0 0 0 1 0];
C1 = [0 0 1 0 0];
C2 = [0 0 1 0 0];

% loop (aka clock) for one sequence length
% this actuall creates a map with each row being the register values at
% that clock cycle
for i=1:N-1   
    % 1 + x2 + x3 + x4 + x5
    A1 = [A1; xor(xor(xor(A1(i,5),A1(i,4)),A1(i,3)),A1(i,2)) A1(i,1) A1(i,2) A1(i,3) A1(i,4)];
    % 1 + x4 + x5
    A2 = [A2; xor(A2(i,5),A2(i,3)) A2(i,1) A2(i,2) A2(i,3) A2(i,4)];
    
    % 1 + x2 + x3 + x4 + x5
    B1 = [B1; xor(xor(xor(B1(i,5),B1(i,4)),B1(i,3)),B1(i,2)) B1(i,1) B1(i,2) B1(i,3) B1(i,4)];
    % 1 + x3 + x5
    B2 = [B2; xor(B2(i,5),B2(i,2)) B2(i,1) B2(i,2) B2(i,3) B2(i,4)];
    
    % 1 + x2 + x3 + x4 + x5
    C1 = [C1; xor(xor(xor(C1(i,5),C1(i,4)),C1(i,3)),C1(i,2)) C1(i,1) C1(i,2) C1(i,3) C1(i,4)];
    % 1 + x2 + x5
    C2 = [C2; xor(C2(i,5),C2(i,4)) C2(i,1) C2(i,2) C2(i,3) C2(i,4)];

end

% chipping sequences are actually last column of each LFSR map (x5 value
% for each clock cycle), xor'd between their primary and secondary
% polynomials
% chipA = A(:,5)';
% chipB = B(:,5)';
chipA = [];
chipB = [];
chipC = [];
for i=1:length(A1(:,5))
    chipA = [chipA xor(A1(i,5),A2(i,5))];
    chipB = [chipB xor(B1(i,5),B2(i,5))];
    chipC = [chipC xor(C1(i,5),C2(i,5))];
end

% map unipolar to bipolar (0 => -1)
for i=1:length(chipA)
    if (chipA(i) == 0)
        chipA(i) = -1;
    end
    if (chipB(i) == 0)
        chipB(i) = -1;
    end
    if (chipC(i) == 0)
        chipC(i) = -1;
    end
end

% baseband messages
Fs = 10;
T = 1/Fs;
t = 0:T/31:((10*T)-(T/31));
messageA = sin(2*pi*t);
messageB = sin(4*pi*t);
messageC = sin(8*pi*t);

% % need to increase sample rate of chipping sequence, 1 sample per bit is
% % not good to match the sample rate of the message
% sampled_chipA = [];
% sampled_chipB = [];
% padding = 1;            % the larger this is, the crappier the signal recovery since it decreases chips per baseband signal bit
% for i=1:length(chipA)
%     if (chipA(i) == 1)
%         sampled_chipA = [sampled_chipA ones(1,padding)];
%     else
%         sampled_chipA = [sampled_chipA -1*ones(1,padding)];
%     end
%     if (chipB(i) == 1)
%         sampled_chipB = [sampled_chipB ones(1,padding)];
%     else
%         sampled_chipB = [sampled_chipB -1*ones(1,padding)];
%     end
% end

% actually instead of expanding single chipping sequence to be = to
% 1 message period, repeat the chipping sequence over the message length as
% much as needed
sampled_chipA = [];
sampled_chipB = [];
sampled_chipC = [];
while length(sampled_chipA) < length(messageA)
    sampled_chipA = [sampled_chipA chipA];
    sampled_chipB = [sampled_chipB chipB];
    sampled_chipC = [sampled_chipC chipC];
end

% discard any runover bits.... shouldn't really matter
if length(sampled_chipA) > length(messageA)
    sampled_chipA = sampled_chipA(1:length(messageA))
    sampled_chipB = sampled_chipB(1:length(messageB))
    sampled_chipC = sampled_chipC(1:length(messageC))
end

% spread message with chipping sequence
codeA = sampled_chipA .* messageA;
codeB = sampled_chipB .* messageB;
codeC = sampled_chipC .* messageC;

% combine different encoded signals
signal = codeA + codeB + codeC;

% despread received signal with expected chipping sequence
receive_signalA = signal .* sampled_chipA;
receive_signalB = signal .* sampled_chipB;
receive_signalC = signal .* sampled_chipC;

% low pass filter despread signal
averaging_length = 31;
filtered_signalA = zeros(1,averaging_length);
filtered_signalB = zeros(1,averaging_length);
filtered_signalC = zeros(1,averaging_length);
for i=averaging_length+1:length(receive_signalA)
    filtered_signalA = [filtered_signalA sum(receive_signalA((i-averaging_length):i-1))/averaging_length];
    filtered_signalB = [filtered_signalB sum(receive_signalB((i-averaging_length):i-1))/averaging_length];
    filtered_signalC = [filtered_signalC sum(receive_signalC((i-averaging_length):i-1))/averaging_length];
end

%% Plots

% transmissions
subplot(4,3,1)
plot(t,messageA)
title ('Message A');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,3,2)
plot(t,messageB)
title ('Message B');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,3,3)
plot(t,messageC)
title ('Message C');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,3,4)
plot(t,sampled_chipA)
title ('Chipping Sequence A');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,3,5)
plot(t,sampled_chipB)
title ('Chipping Sequence B');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,3,6)
plot(t,sampled_chipC)
title ('Chipping Sequence C');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,3,7)
plot(t,codeA)
title ('Spread Signal A');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,3,8)
plot(t,codeB)
title ('Spread Signal B');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,3,9)
plot(t,codeC)
title ('Spread Signal C');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(4,3,10:12)
plot(t,signal)
title ('Combined RF Transmitted Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');

% signal A
figure()
subplot(3,3,1:3)
plot(t,signal)
title ('Combined RF Received Signal');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(3,3,4)
plot(t,receive_signalA)
title ('Despread Message A');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(3,3,5)
plot(t,receive_signalB)
title ('Despread Message B');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(3,3,6)
plot(t,receive_signalC)
title ('Despread Message C');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(3,3,7)
plot(t,filtered_signalA)
title ('Filtered/Recovered Message A');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
hold on
plot(t,messageA)
title ('Original Message A');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(3,3,8)
plot(t,filtered_signalB)
title ('Filtered/Recovered Message B');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
hold on
plot(t,messageB)
title ('Original Message B');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
subplot(3,3,9)
plot(t,filtered_signalC)
title ('Filtered/Recovered Message C');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');
hold on
plot(t,messageC)
title ('Original Message C');
xlabel ('time(sec)');
ylabel ('Amplitude(volt)');