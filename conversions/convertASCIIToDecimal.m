function decimal = convertASCIIToDecimal(ascii)

    decimal = zeros(1,length(ascii));
    for i=1:length(ascii)
        decimal(i) = ascii(i);
    end
end