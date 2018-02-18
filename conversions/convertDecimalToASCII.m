function ascii = convertDecimalToASCII(decimal)

    ascii = [];
    for i=1:length(decimal)
       
        switch decimal(i)
            case 33
                ascii = [ascii '!'];
            case 34
                ascii = [ascii '"'];
            case 35
                ascii = [ascii '#'];
            case 36
                ascii = [ascii '$'];
            case 37
                ascii = [ascii '%'];
            case 38
                ascii = [ascii '&'];
            case 39
                ascii = [ascii "'"];
            case 40
                ascii = [ascii '('];
            case 41
                ascii = [ascii ')'];
            case 42
                ascii = [ascii '*'];
            case 43
                ascii = [ascii '+'];
            case 44
                ascii = [ascii ','];
            case 45
                ascii = [ascii '-'];
            case 46
                ascii = [ascii '.'];
            case 47
                ascii = [ascii '/'];
            case 48
                ascii = [ascii '0'];
            case 49
                ascii = [ascii '1'];
            case 50
                ascii = [ascii '2'];
            case 51
                ascii = [ascii '3'];
            case 52
                ascii = [ascii '4'];
            case 53
                ascii = [ascii '5'];
            case 54
                ascii = [ascii '6'];
            case 55
                ascii = [ascii '7'];
            case 56
                ascii = [ascii '8'];
            case 57
                ascii = [ascii '9'];
            case 58
                ascii = [ascii ':'];
            case 59
                ascii = [ascii ';'];
            case 60
                ascii = [ascii '<'];
            case 61
                ascii = [ascii '='];
            case 62
                ascii = [ascii '>'];
            case 63
                ascii = [ascii '?'];
            case 64
                ascii = [ascii '@'];
            case 65
                ascii = [ascii 'A'];
            case 66
                ascii = [ascii 'B'];
            case 67
                ascii = [ascii 'C'];
            case 68
                ascii = [ascii 'D'];
            case 69
                ascii = [ascii 'E'];
            case 70
                ascii = [ascii 'F'];
            case 71
                ascii = [ascii 'G'];
            case 72
                ascii = [ascii 'H'];
            case 73
                ascii = [ascii 'I'];
            case 74
                ascii = [ascii 'J'];
            case 75
                ascii = [ascii 'K'];
            case 76
                ascii = [ascii 'L'];
            case 77
                ascii = [ascii 'M'];
            case 78
                ascii = [ascii 'N'];
            case 79
                ascii = [ascii 'O'];
            case 80
                ascii = [ascii 'P'];
            case 81
                ascii = [ascii 'Q'];
            case 82
                ascii = [ascii 'R'];
            case 83
                ascii = [ascii 'S'];
            case 84
                ascii = [ascii 'T'];
            case 85
                ascii = [ascii 'U'];
            case 86
                ascii = [ascii 'V'];
            case 87
                ascii = [ascii 'W'];
            case 88
                ascii = [ascii 'X'];
            case 89
                ascii = [ascii 'Y'];
            case 90
                ascii = [ascii 'Z'];
            case 91
                ascii = [ascii '['];
            case 92
                ascii = [ascii '\'];
            case 93
                ascii = [ascii ']'];
            case 94
                ascii = [ascii '^'];
            case 95
                ascii = [ascii '_'];
            case 96
                ascii = [ascii '`'];
            case 97
                ascii = [ascii 'a'];
            case 98
                ascii = [ascii 'b'];
            case 99
                ascii = [ascii 'c'];
            case 100
                ascii = [ascii 'd'];
            case 101
                ascii = [ascii 'e'];
            case 102
                ascii = [ascii 'f'];
            case 103
                ascii = [ascii 'g'];
            case 104
                ascii = [ascii 'h'];
            case 105
                ascii = [ascii 'i'];
            case 106
                ascii = [ascii 'j'];
            case 107
                ascii = [ascii 'k'];
            case 108
                ascii = [ascii 'l'];
            case 109
                ascii = [ascii 'm'];
            case 110
                ascii = [ascii 'n'];
            case 111
                ascii = [ascii 'o'];
            case 112
                ascii = [ascii 'p'];
            case 113
                ascii = [ascii 'q'];
            case 114
                ascii = [ascii 'r'];
            case 115
                ascii = [ascii 's'];
            case 116
                ascii = [ascii 't'];
            case 117
                ascii = [ascii 'u'];
            case 118
                ascii = [ascii 'v'];
            case 119
                ascii = [ascii 'w'];
            case 120
                ascii = [ascii 'x'];
            case 121
                ascii = [ascii 'y'];
            case 122
                ascii = [ascii 'z'];
            case 123
                ascii = [ascii '{'];
            case 124
                ascii = [ascii '|'];
            case 125
                ascii = [ascii '}'];
            case 126
                ascii = [ascii '~'];
            otherwise
                ascii = [ascii ' '];
%                 fprintf('Uknown ASCII character\n')
        end
    end
end