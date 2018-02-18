%% EN525.783 Module 4-3
%  Author: Tim Cardenuto

clear

L = 4;          % # of bits, or shift registers
N = 2^L - 1;    % # of codes, also sequence repeat length in bits

LFSR_complete_fill = [0 0 0 1;
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
             1 1 1 1;];
            
% record matrix
LFSR_record = [LFSR_complete_fill, zeros(15,1)];

         
% initial register values (the 'seed')
LFSR_init = LFSR_complete_fill;

for k=1:length(LFSR_init(:,1))
    
    LFSR = LFSR_init(k,:);
    % loop (aka clock) for one sequence length
    for i=1:N-1       
        % non-maximal length LFSR
        % 1 + x2 + x4
%         LFSR = [LFSR; xor(xor(xor(LFSR(i,4),LFSR(i,3)),LFSR(i,2)),LFSR(i,1)) LFSR(i,1) LFSR(i,2) LFSR(i,3)];
        LFSR = [LFSR; xor(LFSR(i,4),LFSR(i,2)) LFSR(i,1) LFSR(i,2) LFSR(i,3)];
    end

    return
    
    % loop through each row in LFSR_fill (all possible 4th order register fill combinations)      
    % and see if they're each in the generated LFSR matrix.
    % Easiest way to compare is to turn each row into string and do string comparison
    for i=1:length(LFSR_complete_fill(:,1))

        for j=1:length(LFSR(:,1))

            % if we find a match then record the initial fill number required for that row
            if(strcmp(num2str(LFSR_complete_fill(i,:)), num2str(LFSR(j,:))))
                if(LFSR_record(i,5) == 0)
                   LFSR_record(i,5) = k;
                end
            end
        end
    end
end
    