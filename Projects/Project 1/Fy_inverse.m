function output = Fy_inverse(input)
    if input < 0
        output = 0;
    elseif 0 <= input && input < (25/90)
        output = sqrt(90*input) - 9;
    elseif (25/90) <= input && input < (13/18)
        output = (9*input) - (13/2);
    elseif (13/18) <= input && input <= 1
        output = 5 - sqrt(90*(1-input));
    else
        output = 0;
end