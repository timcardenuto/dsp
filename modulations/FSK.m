clear;

M = 2;        % Number of symbols (modulation order)
k = log2(M);  % Number of bits per symbol
n = 1000;    % Number of bits to process
freqsep = 115000;  % Frequency separation (Hz)
numSamplesPerSymbol = 2;    % Number of samples per symbol
Fs = (M-1)*freqsep;      % Sample rate (Hz)

% Generate random M-ary symbols (skipping bit/symbol map)
dataSymbolsIn = randi([0 M-1],n,1);

% Apply FSK modulation
dataMod = fskmod(dataSymbolsIn,M,freqsep,numSamplesPerSymbol,Fs);

% h = dsp.SpectrumAnalyzer('SampleRate',Fs);
% step(h,dataMod)

% Calculate the SNR when the channel has an Eb/N0 = 10 dB.
EbNo = 10;
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);

% Pass the signal through the AWGN channel for both the binary and Gray coded symbol mappings.
receivedSignal = awgn(dataMod,snr,'measured');

sPlotFig = scatterplot(receivedSignal,1,0,'g.');
hold on
scatterplot(dataMod,1,0,'k*',sPlotFig)

