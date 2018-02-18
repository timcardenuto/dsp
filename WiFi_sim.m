M = 4;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 100;                  % Number of bits to process
numSamplesPerSymbol = 1;    % Oversampling factor

% Create a binary data stream as a column vector
%rng default                 % Use default random number generator
bits = randi([0 1],n,1);  % Generate vector of binary data

% Plot the first 40 bits in a stem plot
stem(bits(1:40),'filled');
axis([0 40 -1 2])
title('Pseudorandom Bits');
xlabel('Bit Index');
ylabel('Binary Value');

% Perform a bit-to-symbol mapping
dataInMatrix = reshape(bits,length(bits)/k,k);   % Reshape data into binary 4-tuples
dataSymbols = bi2de(dataInMatrix);                 % Convert to integers

% Plot the first 10 symbols in a stem plot
figure; % Create new figure window.
stem(dataSymbols(1:10));
title('Pseudorandom Symbols');
xlabel('Symbol Index');
ylabel('Integer Value');

% Apply modulation/coding to get complex output....
dataMod = pskmod(dataSymbols,M,0);         % Binary coding, phase offset = 0
dataModG = pskmod(dataSymbols,M,0,'gray'); % Gray coding, phase offset = 0

%figure;
scatterplot(dataMod);
title('Pseudorandom PSK modulated data');
xlabel('Inphase');
ylabel('Quadrature');


% How do I get RF time/freq data?
T = 100;
t = [2*pi/T:2*pi/T:2*pi];
data = [];
carrier = [];
for i = 1:length(bits)
    if bits(i)==0
        data = [data -ones(1,T)];
    else
        data = [data ones(1,T)];
    end
    carrier = [carrier cos(t)]; % extend carrier signal for each bit/symbol
end

rf = data.*carrier;
figure();
plot(rf)

break;

%%
% Calculate the SNR when the channel has an Eb/N0 = 10 dB.
EbNo = 0;
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);

% Pass the signal through the AWGN channel for both the binary and Gray coded symbol mappings.
receivedSignal = awgn(dataMod,snr,'measured');
receivedSignalG = awgn(dataModG,snr,'measured');

sPlotFig = scatterplot(receivedSignal,1,0,'g.');
hold on
scatterplot(dataMod,1,0,'k*',sPlotFig)
