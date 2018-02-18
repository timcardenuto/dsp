function binary = convertDecimalToBinary(decimal)

    binary = [];
    for i=1:length(decimal)
        
        byte = zeros(1,8);
        
        byte(1) = floor(decimal(i)/128);
        decimal(i) = decimal(i) - 128*byte(1);
        byte(2) = floor(decimal(i)/64);
        decimal(i) = decimal(i) - 64*byte(2);
        byte(3) = floor(decimal(i)/32);
        decimal(i) = decimal(i) - 32*byte(3);
        byte(4) = floor(decimal(i)/16);
        decimal(i) = decimal(i) - 16*byte(4);
        byte(5) = floor(decimal(i)/8);
        decimal(i) = decimal(i) - 8*byte(5);
        byte(6) = floor(decimal(i)/4);
        decimal(i) = decimal(i) - 4*byte(6);
        byte(7) = floor(decimal(i)/2);
        decimal(i) = decimal(i) - 2*byte(7);
        byte(8) = floor(decimal(i)/1);
        decimal(i) = decimal(i) - 1*byte(8);
    
        binary = [binary byte];
    end
end