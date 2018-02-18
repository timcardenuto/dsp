% only supports looping from 65-90 (upper case) and 97-122 (lower case) letters
% other decimal numbers representing other ASCII aren't handled right now

function shifted_decimal = shiftRightDecimalASCII(decimal, shift)
    
    shifted_decimal = zeros(1,length(decimal));

    for i=1:length(decimal)
        % lower case shift
        if (96 < decimal(i) && decimal(i) < 123)
%             fprintf("lower case\n")
            shifted_decimal(i) = decimal(i) + shift;
            if (shifted_decimal(i) > 122)
                shifted_decimal(i) = 96 + mod(shifted_decimal(i),122);
            end
        % upper case shift
        elseif (64 < decimal(i) && decimal(i) < 91)
%             fprintf("upper case\n")
            shifted_decimal(i) = decimal(i) + shift;
            if (shifted_decimal(i) > 90)
                shifted_decimal(i) = 64 + mod(shifted_decimal(i),90);
            end
        else
            shifted_decimal(i) = 0;  % not a alphabet character A-Za-z
        end
    end
end