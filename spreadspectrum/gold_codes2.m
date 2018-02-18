%% EN525.783 Module 5-1
%  Author: Tim Cardenuto
%  5th order LFSR

L = 5;          % # of bits, or shift registers
N = 2^L - 1;    % # of codes, also sequence repeat length in bits

% initial register values (the 'seed')
C = [1 1 1 1 1];

% loop (aka clock) for one sequence length
for i=1:N-1       
    % 1 + x2 + x5
    % remember, index in matlab starts at 1 so this looks a little wierd
    C = [C; C(i,2) C(i,3) C(i,4) C(i,5) xor(C(i,3),C(i,1)) ];
end

% chipping sequences are actually last column of each LFSR map
chipC = C(:,1)'
