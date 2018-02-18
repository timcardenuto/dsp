clear
clf

% Fmax = 35;
% T = 1/Fmax;                                 % carrier period
% sample_rate = Fmax*2;                       % # of samples per carrier period (T)
% t = 0:T/sample_rate:(30*T - T/sample_rate); % same # of samples no matter what Fc

t = 1:1:200;

for Fc = 32:36;                       % carrier frequency
    
%     sigA = sin(2*pi*t*Fc);          % no idea why this doesn't work when using real time steps
    sigA = sin(2*pi*t/Fc);


% Filter
% DON'T change the filter properties, you test the output response by
% altering the input frequency
% if sample_rate = 20, and Fc = 10 Hz then interval should be 20 * 10 = 200 sample interval 
    
    resonant_freq = 36;
    % interval should really = filter resonant frequency * sample_rate 
    % (in homework sample rate is usually 1 per second so this term was invisible)
    interval = resonant_freq;     % sample average interval 
    % number of samples used in average should really = filter bandwidth * sample_rate
    % (basically the number of 'cycles' to include)
    num_samples = 4;              % num samples to use in average
    sigA_out = [];
    for i=1:length(t)

        sumA = sigA(i);
        for j=1:(num_samples-1)
            % what do you do when run out of samples in signal for the averaging?
            % use zero? stop running experiment cuz this wouldn't happen IRL?
            if (i+(j*interval)) > length(t)   
                sumA = sumA + 0;
            else  % average normally
                sumA = sumA + sigA(i+(j*interval));
            end
        end
        avrgA = sumA / num_samples;

        sigA_out = [sigA_out avrgA];
    end


    % Plot
    clf
    plot(t,sigA)
    hold on
    plot(t,sigA_out)
    xlim([0 resonant_freq]); 
    title ('Bandpass Filter');
    xlabel('Time(sec)');
    ylabel('Amplitude(volt)');
    legend(['Input Signal w/ Frequency = ',num2str(Fc),'Hz'],['Output Signal w/ resonance = ',num2str(resonant_freq),'Hz']);
    drawnow
    pause(1)
    
end

mag = abs(max(sigA_out))/abs(max(sigA))
[c,lag]=xcorr(sigA,sigA_out);
[maxC,I]=max(c);
sample_lag = lag(I)
