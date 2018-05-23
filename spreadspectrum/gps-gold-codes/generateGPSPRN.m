% GPS PRN Generator
% Outputs 1023 bit Gold code created from modulo-2 sum of two 1023 bit
% sequences G1 and G2. These LFSR's are always the same polynomials since
% we're hard coded for GPS application
% The taps argument should be a two value array for the G2 taps which 
% corresponds to a specific satellite PRN, ex. taps = [2, 6] is PRN 1

function PRN = generateGPSPRN(taps)

L = 10;         % # of bits, or shift registers
N = 2^L - 1;    % # of codes, also sequence repeat length in bits

% initialization vector for G1 and G2 are all ones
G1 = ones(1,L);
G2 = ones(1,L);
PRN = zeros(1,N);

% loop (aka clock) for one sequence length
% this actuall creates a map with each row being the register values at that clock cycle
for i=1:N
    % G1 = 1 + x3 + x10
    G1 = [xor(G1(10),G1(3)) G1(1) G1(2) G1(3) G1(4) G1(5) G1(6) G1(7) G1(8) G1(9)];

    % G2 = 1 + x2 + x3 + x6 + x8 + x9 + x10
    G2 = [xor(xor(xor(xor(xor(G2(10),G2(9)),G2(8)),G2(6)),G2(3)),G2(2)) G2(1) G2(2) G2(3) G2(4) G2(5) G2(6) G2(7) G2(8) G2(9)];

    % PRN 1 = xor of taps G1(10) and G2(2,6)
    PRN(i) = xor(G1(10),xor(G2(taps(1)),G2(taps(2))));
end
