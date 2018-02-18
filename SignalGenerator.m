% Time data
time=[.0001:.0001:.1];
u = rand(1000,1);
plot(time,u)

% Frequency data
n = 2^nextpow2(length(u)); % find power of two larger than data length to optimize fft
y = fft(u,n);
f = 1000/2*linspace(0,1,n/2+1);
figure()
plot(f,2*abs(y(1:n/2+1)/2))