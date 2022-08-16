%PART1
%% 2
close all;
clear all;
rand('seed',sum(100*clock));

%% 3
X1 = 5.*rand(1,1000000); % a=5
X2 = -9 + 9.*rand(1,1000000); % b=9
Y = X1 + X2;

figure;
histogram(Y, 100);
xlabel("y");
ylabel("f_Y(y)");
title("Histogram of 1000000 samples of Y");
power = sum(Y.^2)/1000000;
display(power);

%% 4
N = 8;
maxval = max(Y);
minval = min(Y);

delta = (maxval-minval)/N;
a_initial = [minval:delta:maxval];
recon = [minval+(delta/2):delta:maxval-(delta/2)];

mask = ones(1,1000000);

Y_tilde = recon(1)*(mask.*(Y <= a_initial(2))) + recon(2)*(mask.*(Y > a_initial(2) & Y <= a_initial(3))) + recon(3)*(mask.*(Y > a_initial(3) & Y <= a_initial(4))) + recon(4)*(mask.*(Y > a_initial(4) & Y <= a_initial(5))) + recon(5)*(mask.*(Y > a_initial(5) & Y <= a_initial(6))) + recon(6)*(mask.*(Y > a_initial(6) & Y <= a_initial(7))) + recon(7)*(mask.*(Y > a_initial(7) & Y <= a_initial(8))) + recon(8)*(mask.*(Y > a_initial(8)));

e = Y - Y_tilde;

error_power = sum(e.^2)/1000000;
display(error_power);

sqnr = power/error_power;
sqnr_db = 10*log10(sqnr);
display(sqnr_db);

%% 5
values = recon;
probs = [sum(Y_tilde == recon(1)), sum(Y_tilde == recon(2)), sum(Y_tilde == recon(3)), sum(Y_tilde == recon(4)), sum(Y_tilde == recon(5)), sum(Y_tilde == recon(6)), sum(Y_tilde == recon(7)), sum(Y_tilde == recon(8))]./1000000;

figure;
bar(values,probs);
xlabel("y_{hat}");
ylabel("p_Y_{tilde}(y_{hat})");
title("PMF of the possible quantization outputs");

% Encoder scheme

[sorted_probs, indices] = sort(probs);
[huffman_table, huffman_codes] = HuffmanDeneme(sorted_probs);

codewords = huffman_codes(:, 1);
avg_code_len_huffman = 0;

for i = 1:size(codewords,1)
    disp([sorted_probs(i), codewords(i)]);
    avg_code_len_huffman = avg_code_len_huffman + strlength(codewords(i))*sorted_probs(i);
end

disp("Average code length with Huffman coding: " + string(avg_code_len_huffman));

avg_code_len_vanilla = sum(3*sorted_probs);
disp("Average code length without Huffman coding: " + string(avg_code_len_vanilla));

recon_sorted = recon(indices);
codewords = flip(codewords);

for index = 1:size(Y_tilde, 2)
    temp_codeword = codewords(recon_sorted == Y_tilde(index));
    encoded(index) = temp_codeword;
end

encoded = strjoin(encoded);
encoded = char(regexprep(encoded, ' ', '')); % Encoded string

%% 6
% Decoder scheme

temp_str = "";
decoder_index = 1;

for index = 1:strlength(encoded)

    temp_str = strcat(temp_str, encoded(index));
    
    if ismember(temp_str, codewords)      
        temp_recon = recon_sorted(codewords == temp_str);
        decoded(decoder_index) = temp_recon;
        disp(decoder_index);
        temp_str = "";
        decoder_index = decoder_index + 1;
    end
    
end

if sum(decoded == Y_tilde)/1000000 == 1
    disp("Original input sequence and the output of the decoder are the same.");
else
    disp("Mismatch between the input sequence and the output of the decoder.")
end

%% 8
%Implement both Fy and Fy_inverse as functions. Check files.
%But I will use array manipulation, since it is faster to do it this way,
%instead of using for loops. Also check the report for the functions.

%% 9
compressor = ((-9 <= Y) & (Y < -4)).*(((Y+9).^2)/90) + ((-4 <= Y) & (Y < 0)).*((5/18) + ((Y+4)/9)) + ((0 <= Y) & (Y <= 5)).*((13/18) + (Y/18) + (((5-Y).*Y)/90));

N = 8;
maxval = max(compressor);
minval = min(compressor);

delta = (maxval-minval)/N;
a = [minval:delta:maxval];
recon = [minval+(delta/2):delta:maxval-(delta/2)];

%% 10

y_example = linspace(-9,5,300);
qy = (y_example <= 5 & y_example >= (5- sqrt(90*(1-a(8))))).*recon(8) + (y_example < (5- sqrt(90*(1-a(8)))) & y_example >= (5- sqrt(90*(1-a(7))))).*recon(7) + (y_example < (5- sqrt(90*(1-a(7)))) & y_example >= (9*a(6))-(13/2)).*recon(6) + (y_example < (9*a(6))-(13/2) & y_example >= (9*a(5))-(13/2)).*recon(5) + (y_example < (9*a(5))-(13/2) & y_example >= (9*a(4))-(13/2)).*recon(4) + (y_example < (9*a(4))-(13/2) & y_example >= (sqrt(90*a(3))-9)).*recon(3) + (y_example < (sqrt(90*a(3))-9) & y_example >= (sqrt(90*a(2))-9)).*recon(2) + (y_example < (sqrt(90*a(2))-9) & y_example >= (sqrt(90*a(1))-9)).*recon(1);
plot(y_example, qy);
xlabel("y");
ylabel("Q(y)");

%% 11

mask = ones(1,1000000);
Y_tilde_compressed = recon(1)*(mask.*(compressor <= a(2))) + recon(2)*(mask.*(compressor > a(2) & compressor <= a(3))) + recon(3)*(mask.*(compressor > a(3) & compressor <= a(4))) + recon(4)*(mask.*(compressor > a(4) & compressor <= a(5))) + recon(5)*(mask.*(compressor > a(5) & compressor <= a(6))) + recon(6)*(mask.*(compressor > a(6) & compressor <= a(7))) + recon(7)*(mask.*(compressor > a(7) & compressor <= a(8))) + recon(8)*(mask.*(compressor > a(8)));

expander = ((0 <= Y_tilde_compressed) & (Y_tilde_compressed < (25/90))).*(sqrt(90*Y_tilde_compressed)-9) + (((25/90) <= Y_tilde_compressed) & (Y_tilde_compressed < (13/18))).*(9*Y_tilde_compressed - (13/2)) + (((13/18) <= Y_tilde_compressed) & (Y_tilde_compressed <= 1)).*(5-sqrt(90*(1-Y_tilde_compressed)));

e_compressed = Y - expander;
error_avg_pow_compressed = sum(e_compressed.^2)/1000000;
display(error_avg_pow_compressed);

sqnr_compressed = power/error_avg_pow_compressed;
sqnr_compressed_db = 10*log10(sqnr_compressed);
display(sqnr_compressed_db);

%% 12

values = unique(expander);

probs = [sum(expander == values(1)), sum(expander == values(2)), sum(expander == values(3)), sum(expander == values(4)), sum(expander == values(5)), sum(expander == values(6)), sum(expander == values(7)), sum(expander == values(8))]./1000000;

figure;
bar(values,probs);
xlabel("y_{hat}");
ylabel("p_Y_{tilde}(y_{hat})");
title("PMF of the possible quantization outputs");

% Encoder scheme

[sorted_probs, indices] = sort(probs);
[huffman_table2, huffman_codes2] = HuffmanDeneme(sorted_probs);

codewords = huffman_codes2(:, 1);
avg_code_len_huffman2 = 0;

for i = 1:size(codewords,1)
    disp([sorted_probs(i), codewords(i)]);
    avg_code_len_huffman2 = avg_code_len_huffman2 + strlength(codewords(i))*sorted_probs(i);
end

disp("Average code length with Huffman coding: " + string(avg_code_len_huffman2));

avg_code_len_vanilla2 = sum(3*sorted_probs);
disp("Average code length without Huffman coding: " + string(avg_code_len_vanilla2));

values_sorted = values(indices);
codewords = flip(codewords);

encoded = "";

for index = 1:size(Y_tilde_compressed, 2)
    temp_codeword = codewords(values_sorted == expander(index));
    encoded(index) = temp_codeword;
end

encoded = strjoin(encoded);
encoded = char(regexprep(encoded, ' ', '')); % Encoded string


%% 13

% Decoder scheme

temp_str = "";
decoder_index = 1;

for index = 1:strlength(encoded)

    temp_str = strcat(temp_str, encoded(index));
    
    if ismember(temp_str, codewords)      
        temp_recon = recon_sorted(codewords == temp_str);
        decoded(decoder_index) = temp_recon;
        disp(decoder_index);
        temp_str = "";
        decoder_index = decoder_index + 1;
    end
    
end

if sum(decoded == Y_tilde_compressed)/1000000 == 1
    disp("Original input sequence and the output of the decoder are the same.");
else
    disp("Mismatch between the input sequence and the output of the decoder.")
end

%% 14

[a_final, recon_final] = LloydMaxFunc(Y, a_initial);

%% 15

y_example = linspace(-9,5,300);
q_lloyd = (y_example <= 5 & y_example >= a_final(8)).*recon_final(8) + ((y_example < a_final(8)) & (y_example >= a_final(7))).*recon_final(7) + ((y_example < a_final(7)) & y_example >= a_final(6)).*recon_final(6) + ((y_example < a_final(6)) & y_example >= a_final(5)).*recon_final(5) + (y_example < a_final(5) & y_example >= a_final(4)).*recon_final(4) + (y_example < a_final(4) & y_example >= a_final(3)).*recon_final(3) + (y_example < a_final(3) & y_example >= a_final(2)).*recon_final(2) + (y_example < a_final(2)).*recon_final(1);
plot(y_example, q_lloyd);

%% 16

Y_tilde_lloyd = recon_final(1)*(mask.*(Y <= a_final(2))) + recon_final(2)*(mask.*(Y > a_final(2) & Y <= a_final(3))) + recon_final(3)*(mask.*(Y > a_final(3) & Y <= a_final(4))) + recon_final(4)*(mask.*(Y > a_final(4) & Y <= a_final(5))) + recon_final(5)*(mask.*(Y > a_final(5) & Y <= a_final(6))) + recon_final(6)*(mask.*(Y > a_final(6) & Y <= a_final(7))) + recon_final(7)*(mask.*(Y > a_final(7) & Y <= a_final(8))) + recon_final(8)*(mask.*(Y > a_final(8)));

e_lloyd = Y - Y_tilde_lloyd;

epower_lloyd = sum(e_lloyd.^2)/1000000;
display(epower_lloyd);
sqnr_lloyd = power/epower_lloyd;
sqnr_lloyd_db = 10*log10(sqnr_lloyd);
display(sqnr_lloyd_db);