function decimal = convertBinaryToDecimalASCII(binary)

    if (mod(length(binary),8) ~= 0)
        % this means there are leftover bits, probably error
    end
    numchar = floor(length(binary)/8);

    decimal = [];
    j = 1;
    for i=1:numchar
        decimal = [decimal (binary(j)*128 + binary(j+1)*64 + binary(j+2)*32 + binary(j+3)*16 + binary(j+4)*8 + binary(j+5)*4 + binary(j+6)*2 + binary(j+7)*1)]; 
        j = j + 8;
    end
end