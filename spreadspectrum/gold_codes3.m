%% EN525.783 Module 5-3
%  Author: Tim Cardenuto
%  3rd order Gold Code

L = 3;          % # of bits, or shift registers
N = 2^L - 1;    % # of codes, also sequence repeat length in bits

% initial register values (the 'seed')
seq = [0 0 1];

% loop (aka clock) for one sequence length
for i=1:N-1   
    % 1 + x + x3
    seq = [seq; seq(i,2) seq(i,3) xor(seq(i,2),seq(i,1))];
end

% chipping sequences are the output column of each LFSR map
chip_org = seq(:,1)';

% decimate original chipping sequence
chip_decimate = [];
j = 1;  % starting bit
decimate_rate = 5;
for i=1:N
    chip_decimate(i) = chip_org(j);
    j = j + decimate_rate;
    if j > N
        j = j - N;
    end
end

% shift decimated version with offset
chip_shift = chip_decimate;
for i=2:N
    for j=1:N-1
        chip_shift(i,j) = chip_shift(i-1,j+1);
    end
    chip_shift(i,j+1) = chip_shift(i-1,1);
end

% and XOR shifted with original
chips = [chip_org; chip_decimate];
for i=1:N
    for j=1:N
        chips(i+2,j) = xor(chip_shift(i,j),chip_org(j));
    end
end

% correlate
corr_mul = chips(1,:) .* chips(1,:);
corr_sum = sum(corr_mul);
for i=2:N+2
    corr_mul = [corr_mul; (chips(1,:) .* chips(i,:))];
    corr_sum = [corr_sum; sum(corr_mul((i),:))];
end

%% 10th order Gold Code
clear

L = 10;          % # of bits, or shift registers
N = 2^L - 1;    % # of codes, also sequence repeat length in bits

% initial register values (the 'seed')
seq = [0 0 0 0 0 0 0 0 0 1];
% seq = [1 1 1 1 1 1 1 1 1 1];  
seq2 = seq;

% loop (aka clock) for one sequence length
for i=1:N-1   
    % 1 + x3 + x10
    seq = [seq; seq(i,2) seq(i,3) seq(i,4) seq(i,5) seq(i,6) seq(i,7) seq(i,8) seq(i,9) seq(i,10) xor(seq(i,4),seq(i,1))];

    % 1 + x2 + x3 + x6 + x8 + x9 + x10
    seq2 = [seq2; seq2(i,2) seq2(i,3) seq2(i,4) seq2(i,5) seq2(i,6) seq2(i,7) seq2(i,8) seq2(i,9) seq2(i,10) xor(seq2(i,10),xor(seq2(i,9),xor(seq2(i,7),xor(seq2(i,4),xor(seq2(i,3),seq2(i,1))))))];

end

% chipping sequences are the output column of each LFSR map
chip_org = seq(:,1)';
chip_org2 = seq2(:,1)';

% compare both original
if (strcmp(num2str(chip_org),num2str(chip_org2)))
        num2str(chip_org)
        num2str(chip_org2)
end

% decimate original chipping sequence
chip_decimate = [];
j = 1;  % starting bit
decimate_rate = 65;
for i=1:N
    chip_decimate(i) = chip_org(j);
    j = j + decimate_rate;
    if j > N
        j = j - N;
    end
end

% compare decimated with original
if (strcmp(num2str(chip_decimate),num2str(chip_org2)))
        num2str(chip_decimate)
        num2str(chip_org2)
end

% shift decimated version with offset
chip_shift = chip_decimate;
for i=2:N
    for j=1:N-1
        chip_shift(i,j) = chip_shift(i-1,j+1);
    end
    chip_shift(i,j+1) = chip_shift(i-1,1);
    
    % check this shifted version
    if (strcmp(num2str(chip_shift(i,:)),num2str(chip_org2)))    % TODO try subset of range to find something close and see what's wrong
        i-1
        num2str(chip_shift(i,:))
        num2str(chip_org2)
    end
end

