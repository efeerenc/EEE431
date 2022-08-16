function output = Fy(input)
    if input < -9
        output = 0;
    elseif -9 <= input && input < -4
        output = ((input+9)^2)/90;
    elseif -4 <= input && input < 0
        output = (5/18) + (input+4)/9;
    elseif 0 <= input && input <= 5
        output = (13/8) + (input/18) + (5-input)*(input)/90;
    else
        output = 0;
    end
end