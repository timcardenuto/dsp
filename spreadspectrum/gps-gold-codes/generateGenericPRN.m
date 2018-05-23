% Generic Linear Feedback Shift Register (LFSR)
% order = polynomial order, aka num of bits/registers, ex. 10
% fill = initial fill for registers, must be length of the polynomial,
%        ex. [0 0 0 0 0 0 0 0 0 1] could be a fill for polynomial order 10
% feedback_taps = registers to xor on every clock cycle for feedback input to 1st register
% output_taps = registers to xor on every clock cycle for the output bit,
%               to just to use the last register bit as output, only enter
%               the order (last register location) for this argument

function output = generateGenericPRN(order, fill, feedback_taps, output_taps)

N = 2^order - 1;    % # of codes, also sequence repeat length in bits

X = fill;  % initialization vector for X registers

output = zeros(1,N);

% loop (aka clock) for one sequence length
for i=1:N    
    % shift bits based on order, xor the feedback taps
    feedback = 0;
    j = order;
    while j > 1;
        if ismember(j,feedback_taps)
            feedback = xor(feedback,X(j));
        end
        X(j) = X(j-1);
        j = j - 1;
    end
    X(1) = feedback;

    % clock output bit into PRN sequence
    for j=1:length(output_taps)
        output(i) = xor(output(i), X(output_taps(j)));
    end
end
