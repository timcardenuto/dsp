%% EN525.783 Module 4-1
%  Author: Tim Cardenuto
%  Description: Linear Feedback Shift Registers (LFSR) 
%  using 3rd, 4th, and 5th order polynomials:
%   1 + x + x3
%   1 + x + x4
%   1 + x2 + x5
clear

%% 3rd order LFSR
L = 3;          % # of bits, or shift registers
N = 2^L - 1;    % # of codes, also sequence repeat length in bits

% initial register values (the 'seed')
A = [0 0 1];

% loop (aka clock) for one sequence length
for i=1:N-1   
    % 1 + x + x3
    A = [A; xor(A(i,3),A(i,1)) A(i,1) A(i,2)];
end

% chipping sequences are actually last column of each LFSR map
chipA = A(:,3)';


%% 4th order LFSR
L = 4;          % # of bits, or shift registers
N = 2^L - 1;    % # of codes, also sequence repeat length in bits

% initial register values (the 'seed')
B = [0 0 0 1];

% loop (aka clock) for one sequence length
for i=1:N-1       
    % 1 + x + x4
    	  
%     B = [B; xor(B(i,4),B(i,1)) B(i,1) B(i,2) B(i,3)];  % my original guess, output on the right
    
%     B = [B; B(i,2) B(i,3) B(i,4) xor(B(i,4),B(i,1))];  % gives the reverse, based on professor output

      B = [B; B(i,2) B(i,3) B(i,4) xor(B(i,2),B(i,1))]; % correct according to professor
end

% chipping sequences are the output column of each LFSR map
% chipB = B(:,4)'; % original guess, output on right, or reverse in the reverse case above
chipB = B(:,1)'; % correct according to professor

%% 5th order LFSR
L = 5;          % # of bits, or shift registers
N = 2^L - 1;    % # of codes, also sequence repeat length in bits

% initial register values (the 'seed')
C = [0 0 0 0 1];

% loop (aka clock) for one sequence length
for i=1:N-1       
    % 1 + x2 + x5
    C = [C; xor(C(i,5),C(i,2)) C(i,1) C(i,2) C(i,3) C(i,4)];
end

% chipping sequences are actually last column of each LFSR map
chipC = C(:,5)';


%% Analysis

B_ordered = [0 0 0 1;
             0 0 1 0;
             0 0 1 1;
             0 1 0 0;
             0 1 0 1;
             0 1 1 0;
             0 1 1 1;
             1 0 0 0;
             1 0 0 1;
             1 0 1 0;
             1 0 1 1;
             1 1 0 0;
             1 1 0 1;
             1 1 1 0;
             1 1 1 1;]
         
B

ones_A = length(find(chipA==1))
zeros_A = length(find(chipA==0))
assert((length(find(chipA==0)) == (length(find(chipA==1))-1)),'number of zeros != number of ones - 1')
ones_B = length(find(chipB==1))
zeros_B = length(find(chipB==0))
assert((length(find(chipB==0)) == (length(find(chipB==1))-1)),'number of zeros != number of ones - 1')
ones_C = length(find(chipC==1))
zeros_C = length(find(chipC==0))
assert((length(find(chipC==0)) == (length(find(chipC==1))-1)),'number of zeros != number of ones - 1')
