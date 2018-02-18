M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 30000;                  % Number of bits to process
numSamplesPerSymbol = 1;    % Oversampling factor

% Create a binary data stream as a column vector
%rng default                 % Use default random number generator
dataIn = randi([0 1],n,1);  % Generate vector of binary data

% Plot the first 40 bits in a stem plot
stem(dataIn(1:40),'filled');
title('Random Bits');
xlabel('Bit Index');
ylabel('Binary Value');

% Perform a bit-to-symbol mapping
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary 4-tuples
dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers

% Plot the first 10 symbols in a stem plot
figure; % Create new figure window.
stem(dataSymbolsIn(1:10));
title('Random Symbols');
xlabel('Symbol Index');
ylabel('Integer Value');

% Apply modulation
dataMod = pskmod(dataSymbolsIn,M,0);         % Binary coding, phase offset = 0
dataModG = pskmod(dataSymbolsIn,M,0,'gray'); % Gray coding, phase offset = 0

% Calculate the SNR when the channel has an Eb/N0 = 10 dB.
EbNo = 20;
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);

% Pass the signal through the AWGN channel for both the binary and Gray coded symbol mappings.
receivedSignal = awgn(dataMod,snr,'measured');
receivedSignalG = awgn(dataModG,snr,'measured');

sPlotFig = scatterplot(receivedSignal,1,0,'g.');
hold on
scatterplot(dataMod,1,0,'k*',sPlotFig)