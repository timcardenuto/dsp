function binary = convertHexToBinary(hex)

    binary = zeros(1,length(hex)*4);
    j = 1;
    for i=1:length(hex)
        switch hex(i)
            case '0'
                binary(j:(j+3)) = [0 0 0 0];
            case '1'
                binary(j:(j+3)) = [0 0 0 1];
            case '2'
                binary(j:(j+3)) = [0 0 1 0];
            case '3'
                binary(j:(j+3)) = [0 0 1 1];
            case '4'
                binary(j:(j+3)) = [0 1 0 0]; 
            case '5'
                binary(j:(j+3)) = [0 1 0 1]; 
            case '6'
                binary(j:(j+3)) = [0 1 1 0]; 
            case '7'
                binary(j:(j+3)) = [0 1 1 1]; 
            case '8'
                binary(j:(j+3)) = [1 0 0 0]; 
            case '9'
                binary(j:(j+3)) = [1 0 0 1]; 
            case 'A'
                binary(j:(j+3)) = [1 0 1 0]; 
            case 'B'
                binary(j:(j+3)) = [1 0 1 1]; 
            case 'C'
                binary(j:(j+3)) = [1 1 0 0]; 
            case 'D'
                binary(j:(j+3)) = [1 1 0 1]; 
            case 'E'
                binary(j:(j+3)) = [1 1 1 0]; 
            case 'F'
                binary(j:(j+3)) = [1 1 1 1]; 
        end
        j = j + 4;
    end
end