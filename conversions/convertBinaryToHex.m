function hex = convertBinaryToHex(binary)

    if (mod(length(binary),4) ~= 0)
        % this means there are leftover bits, probably error
    end
    
    numchars = floor(length(binary)/4);
    hex = [];
    j = 1;
    for i=1:numchars
        
        decimal = binary(j)*8 + binary(j+1)*4 + binary(j+2)*2 + binary(j+3)*1;       

        switch decimal
            case 0
                hex = [hex '0'];
            case 1
                hex = [hex '1'];
            case 2
                hex = [hex '2'];
            case 3
                hex = [hex '3'];
            case 4
                hex = [hex '4'];
            case 5
                hex = [hex '5'];
            case 6
                hex = [hex '6'];
            case 7
                hex = [hex '7'];
            case 8
                hex = [hex '8'];
            case 9
                hex = [hex '9'];
            case 10
                hex = [hex 'A'];
            case 11
                hex = [hex 'B'];
            case 12
                hex = [hex 'C'];
            case 13
                hex = [hex 'D'];
            case 14
                hex = [hex 'E'];
            case 15
                hex = [hex 'F'];
        end
        j = j + 4;
    end
end