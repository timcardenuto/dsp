%% EN525.783 Module 4-2
%  Author: Tim Cardenuto
clf
clear

L = 4;          % # of bits, or shift registers
N = 2^L - 1;    % # of codes, also sequence repeat length in bits

% initial register values (the 'seed')
B = [1 1 0 0];

% loop (aka clock) for one sequence length
for i=1:N-1       
    % 1 + x + x4
    B = [B; xor(B(i,4),B(i,1)) B(i,1) B(i,2) B(i,3)];
end

% chipping sequences are actually last column of each LFSR map
chipB = B(:,4)';

% map unipolar to bipolar (0 => -1)
for i=1:length(chipB)
    if (chipB(i) == 0)
        chipB(i) = -1;
    end
end

% shift and correlate
corr_mul = chipB .* chipB;
corr_sum = sum(corr_mul);
for i=1:length(chipB)

    % shift right
    temp = zeros(1,length(chipB));
    temp = chipB(length(chipB));
    for j=2:length(chipB)
        temp(j) = chipB(j-1);
    end
    chipB = [chipB; temp];    
    
    % correlate
    corr_mul = [corr_mul; (chipB(1,:) .* temp)];
    corr_sum = [corr_sum; sum(corr_mul((i+1),:))];
end


%% Reversed taps
% original taps were x^4 + x^1 + x^0
% m = number of LFSR stages = 4, to reverse taps, take m - tap #
% reverse taps are x^(4-4) + x^(4-1) + x^(4-0) = 1 + x3 + x4

% initial register values (the 'seed')
B_reverse = [0 0 0 1];

% loop (aka clock) for one sequence length
for i=1:N-1       
    % 1 + x3 + x4
    B_reverse = [B_reverse; xor(B_reverse(i,4),B_reverse(i,3)) B_reverse(i,1) B_reverse(i,2) B_reverse(i,3)];
end

% chipping sequences are actually last column of each LFSR map
chipB_reverse = B_reverse(:,4)';

% map unipolar to bipolar (0 => -1)
for i=1:length(chipB_reverse)
    if (chipB_reverse(i) == 0)
        chipB_reverse(i) = -1;
    end
end

% shift and correlate
corr_mul_rev = chipB(1,:) .* chipB_reverse;      % 1st correlation
corr_sum_rev = sum(corr_mul_rev);

% not sure I need to do this? correlate all new sequence shifts against original chipB from part 1?
for i=1:length(chipB_reverse)

    % shift right
    temp = zeros(1,length(chipB_reverse));
    temp = chipB_reverse(length(chipB_reverse));
    for j=2:length(chipB_reverse)
        temp(j) = chipB_reverse(j-1);
    end
    chipB_reverse = [chipB_reverse; temp];    
    
    % correlate
    corr_mul_rev = [corr_mul_rev; (chipB(1,:) .* temp)];
    corr_sum_rev = [corr_sum_rev; sum(corr_mul_rev((i+1),:))];
end

subplot(3,1,1)
stem(chipB(1,:))
title ('4th Order Chipping Sequence');
xlabel ('Chip');
ylabel ('Value');

subplot(3,1,2)
stem(chipB_reverse(1,:))
title ('Reverse 4th Order Chipping Sequence');
xlabel ('Chip');
ylabel ('Value');

n = 0:N;
subplot(3,1,3)
plot(n,corr_sum)
grid on
hold on
plot(n,corr_sum_rev)
title ('Cross-Correlation of Sequences');
xlabel ('# Bits Shifted Right');
ylabel ('Correlation Value');
legend('Original sequence correlated with itself','Reversed sequence correlated with original');
