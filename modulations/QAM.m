M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 1000;                  % Number of bits to process
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
dataMod = qammod(dataSymbolsIn,M,0);         % Binary coding, phase offset = 0
% dataModG = qammod(dataSymbolsIn,M,0,'gray'); % Gray coding, phase offset = 0

% Calculate the SNR when the channel has an Eb/N0 = 10 dB.
EbNo = 5;
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);

% Pass the signal through the AWGN channel for both the binary and Gray coded symbol mappings.
receivedSignal = awgn(dataMod,snr,'measured');
% receivedSignalG = awgn(dataModG,snr,'measured');

sPlotFig = scatterplot(receivedSignal,1,0,'g.');
hold on
scatterplot(dataMod,1,0,'k*',sPlotFig)

% sPlotFig = scatterplot(receivedSignalG,1,0,'b.');
% hold on
% scatterplot(dataModG,1,0,'k*',sPlotFig)

%% Given N copies of the same signal
% How he creates data - generate different randomness to each signal copy
nSym = 1000;
order = 16;
n = wgn(nSym, order, 0, 'complex');
for sigcopy=1:order
    channelgain = raylrnd(1/sqrt(2),nSym,order); % add channel gain
    channelphase = exp(j*2*pi*rand(nSym,order)); % add channel phase
    channelnoise = channelgain.*channelphase; % complex channel noise
    n(:,sigcopy)=n(:,sigcopy)/(sqrt(sum(n(:,sigcopy).*conj(n(:,sigcopy)))/nSym));
end
yn = [sig1, sig2, sig3, sig4];
Th = ; % How do I get "theoretical"?
nSym = size(yn,1);


fft40 = fft(yn(:,1).*yn(:,2).*yn(:,3).*yn(:,4));
fft42 = fft(yn(:,1).*yn(:,2).*conj(yn(:,3)).*conj(yn(:,4)));

fftReg=[fft40];
fftCon=[fft42];

% Calculate moments, find max
[momReg, peakIdxReg] = max(abs(fftReg)/nSym);
[momCon, peakIdxCon] = max(abs(fftCon)/nSym);
momReg=momReg';
momCon=momCon';
[maxPeak, maxPeakIdx] = max(momReg);
outIdx = [maxPeakIdx peakIdxReg(maxPeakIdx)];

% Calculate min distance
vv = abs(momReg./momCon);
for jjj = 1:size(Th,2)
    temp = abs(vv-Th(:,jjj));
    distVec(jjj)=0;
    for iii = 1:length(orderIdx)
        distVec(jjj) = distVec(jjj) + temp(orderIdx(iii));
    end
end
[dum idxm] = min(distVec);

figure(1)
plot(abs(fft40/nSym));
title('m40')
hold on
plot(peakIdxReg(1), abs(fft40(peakIdxReg(1)))/nSym, 's')

figure(2)
plot(abs(fft42/nSym));
title('m42')
hold on
plot(peakIdxCon(1), abs(fft42(peakIdxCon(1)))/nSym, 's')


