%PART1
% READ:
%I appended the extra functions at the end of the code, so to run the
%script correctly please create a seperate file for each of them.
%I uploaded all of them in a single file because the project
%assignment implies that way.
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

% Part 2
%% 1
corpus = 'FEYNMAN WAS BORN ON MAY IN QUEENS NEW YORK CITY TO LUCILLE NE PHILLIPS A HOMEMAKER AND MELVILLE ARTHUR FEYNMAN A SALES MANAGER ORIGINALLY FROM MINSK IN BELARUS THEN PART OF THE RUSSIAN EMPIRE FEYNMAN WAS A LATE TALKER AND DID NOT SPEAK UNTIL AFTER HIS THIRD BIRTHDAY AS AN ADULT HE SPOKE WITH A NEW YORK ACCENT STRONG ENOUGH TO BE PERCEIVED AS AN AFFECTATION OR EXAGGERATION SO MUCH SO THAT HIS FRIENDS WOLFGANG PAULI AND HANS BETHE ONCE COMMENTED THAT FEYNMAN SPOKE LIKE A BUM THE YOUNG FEYNMAN WAS HEAVILY INFLUENCED BY HIS FATHER WHO ENCOURAGED HIM TO ASK QUESTIONS TO CHALLENGE ORTHODOX THINKING AND WHO WAS ALWAYS READY TO TEACH FEYNMAN SOMETHING NEW FROM HIS MOTHER HE GAINED THE SENSE OF HUMOR THAT HE HAD THROUGHOUT HIS LIFE AS A CHILD HE HAD A TALENT FOR ENGINEERING MAINTAINED AN EXPERIMENTAL LABORATORY IN HIS HOME AND DELIGHTED IN REPAIRING RADIOS THIS RADIO REPAIRING WAS PROBABLY THE FIRST JOB FEYNMAN HAD AND DURING THIS TIME HE SHOWED EARLY SIGNS OF AN APTITUDE FOR HIS LATER CAREER IN THEORETICAL PHYSICS WHEN HE WOULD ANALYZE THE ISSUES THEORETICALLY AND ARRIVE AT THE SOLUTIONS WHEN HE WAS IN GRADE SCHOOL HE CREATED A HOME BURGLAR ALARM SYSTEM WHILE HIS PARENTS WERE OUT FOR THE DAY RUNNING ERRANDS WHEN RICHARD WAS FIVE HIS MOTHER GAVE BIRTH TO A YOUNGER BROTHER HENRY PHILLIPS WHO DIED AT AGE FOUR WEEKS FOUR YEARS LATER RICHARDS SISTER JOAN WAS BORN AND THE FAMILY MOVED TO FAR ROCKAWAY QUEENS THOUGH SEPARATED BY NINE YEARS JOAN AND RICHARD WERE CLOSE AND THEY BOTH SHARED A CURIOSITY ABOUT THE WORLD THOUGH THEIR MOTHER THOUGHT WOMEN LACKED THE CAPACITY TO UNDERSTAND SUCH THINGS RICHARD ENCOURAGED JOANS INTEREST IN ASTRONOMY AND JOAN EVENTUALLY BECAME AN ASTROPHYSICIST FEYNMANS PARENTS WERE BOTH FROM JEWISH FAMILIES BUT NOT RELIGIOUS AND BY HIS YOUTH FEYNMAN DESCRIBED HIMSELF AS AN AVOWED ATHEIST MANY YEARS LATER IN A LETTER TO TINA LEVITAN DECLINING A REQUEST FOR INFORMATION FOR HER BOOK ON JEWISH NOBEL PRIZE WINNERS HE STATED TO SELECT FOR APPROBATION THE PECULIAR ELEMENTS THAT COME FROM SOME SUPPOSEDLY JEWISH HEREDITY IS TO OPEN THE DOOR TO ALL KINDS OF NONSENSE ON RACIAL THEORY ADDING AT THIRTEEN I WAS NOT ONLY CONVERTED TO OTHER RELIGIOUS VIEWS BUT I ALSO STOPPED BELIEVING THAT THE JEWISH PEOPLE ARE IN ANY WAY THE CHOSEN PEOPLE LATER IN HIS LIFE DURING A VISIT TO THE JEWISH THEOLOGICAL SEMINARY HE ENCOUNTERED THE TALMUD FOR THE FIRST TIME HE SAW THAT IT CONTAINED THE ORIGINAL TEXT IN A LITTLE SQUARE ON THE PAGE AND SURROUNDING IT WERE COMMENTARIES WRITTEN OVER TIME BY DIFFERENT PEOPLE IN THIS WAY THE TALMUD HAD EVOLVED AND EVERYTHING THAT WAS DISCUSSED WAS CAREFULLY RECORDED DESPITE BEING IMPRESSED FEYNMAN WAS DISAPPOINTED WITH THE LACK OF INTEREST FOR NATURE AND THE OUTSIDE WORLD EXPRESSED BY THE RABBIS WHO CARED ABOUT ONLY THOSE QUESTIONS WHICH ARISE FROM THE TALMUD';
corpus = strcat(corpus, '#');

keyset = {'#' 'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z' ' '};
values = linspace(0,27,28);

dictionary = containers.Map(keyset, values);

first_char = corpus(1);
word_index = 2;
dict_index = 28;

encoded = "";

counts_char = containers.Map();
counts_len = containers.Map();

while word_index <= size(corpus, 2) + 1
    
    if word_index <= size(corpus, 2)
        next_char = corpus(word_index);
    end
    
    if first_char == "#"
        code = dec2bin(dictionary(first_char), ceil(log2(size(dictionary.keys(), 2))));
        encoded = strcat(encoded, string(code));
        if ~isKey(counts_char, char(first_char))
            counts_char(char(first_char)) = 1;
            counts_len(char(first_char)) = size(code, 2);
        else
            counts_char(char(first_char)) = counts_char(char(first_char)) + 1;
            counts_len(char(first_char)) = counts_len(char(first_char)) + size(code, 2);
        end
        word_index = word_index + 1;
    elseif isKey(dictionary, strcat(string(first_char), string(next_char)))
        first_char = strcat(string(first_char), string(next_char));
        word_index = word_index + 1;
    else
        code = dec2bin(dictionary(first_char), ceil(log2(size(dictionary.keys(), 2))));
        encoded = strcat(encoded, string(code));
        
        if ~isKey(counts_char, char(first_char))
            counts_char(char(first_char)) = 1;
            counts_len(char(first_char)) = size(code, 2);
        else
            counts_char(char(first_char)) = counts_char(char(first_char)) + 1;
            counts_len(char(first_char)) = counts_len(char(first_char)) + size(code, 2);
        end
        
        disp(strcat(string(first_char), string(next_char)));
        
        dictionary(strcat(string(first_char), string(next_char))) = dict_index;
        dict_index = dict_index + 1;
        
        first_char = next_char;
        word_index = word_index + 1; 
    end       
end

%% 2 
keys = dictionary.keys();
num = 0;
for i = 1:size(dictionary.keys(), 2)
    num = num + strlength(keys(i));
end

bits = 0;
occur = 0;
code_keys = string(counts_char.keys());
for i = 1:size(code_keys, 2)
    current_key = code_keys(i);
    current_occurrence = counts_char(current_key);
    current_length = counts_len(current_key);
    bits = current_length + bits;
    occur = occur + current_occurrence;
end
    
avg_codeword_length_lzw = bits/occur;
avg_codeword_length = 5;

disp(avg_codeword_length_lzw);

%% 3
keyset = {'#' 'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z' ' '};
values = linspace(0,27,28);
dictionary = containers.Map(keyset, values);

values = dec2bin(values, ceil(log2(size(dictionary.keys(), 2))));
inverse_dict = containers.Map(string(values), keyset);

bitstream = char(encoded);
decoded = "";

first_code = bitstream(1:5);
decoded = strcat(decoded, string(inverse_dict(first_code)));

character = first_code;

code_index = 2;
code_size = 5;

dict_index = 28;

start_index = 6;
last_index = 10;

finished = 0;

while (code_index < size(corpus, 2)) && finished == 0

    %new_code = bitstream((code_index-1)*code_size+1:code_index*code_size);     
    
    new_code = bitstream(start_index:last_index);
    if ~isKey(inverse_dict, new_code)
        text = inverse_dict(first_code);
        text = strcat(string(text), string(character));
        %code_size = code_size + 1;
    else
        text = inverse_dict(new_code);      
    end
    
    decoded = strcat(decoded, string(text));
    incase = char(text);
    character = incase(1);
    
    if rem(log2(size(inverse_dict.keys(),2)), 1) == 0
        first_code = strcat("0", first_code);
        inverse_dict((dec2bin(dict_index, (ceil(log2(size(inverse_dict.keys(), 2))))+1))) = strcat(string(inverse_dict(first_code)), string(character));
    else
        inverse_dict((dec2bin(dict_index, ceil(log2(size(inverse_dict.keys(), 2)))))) = strcat(string(inverse_dict(first_code)), string(character));
    end
    
    
    dict_index = dict_index + 1;
    
    if rem(log2(size(inverse_dict.keys(),2)), 1) == 0
        code_size = log2(size(inverse_dict.keys(),2)) + 1;
        
        keys = inverse_dict.keys();
        
        for i = 1:size(keys,2)
            value = inverse_dict(char(keys(i)));
            new_key = char(strcat("0",string(keys(i))));
            remove(inverse_dict, keys(i));
            inverse_dict(new_key) = value;
        end
        
    else
        code_size = ceil(log2(size(inverse_dict.keys(),2)));
    end    
    
    first_code = new_code;
    code_index = code_index + 1;
    
    start_index = last_index + 1;
    last_index = start_index + code_size - 1;
    
    checking = char(decoded);
    
    if checking(end) == "#"
        finished = 1;
    end
    
end

if decoded == corpus
    disp("Input sequence and the decoder output are the same.");
else
    disp("Error.");
end

%% 4
for i = 1:size(keyset, 2)
    occurrence(i) = length(find(corpus == char(keyset(i))));
end

xaxis = linspace(1,28,28);

prob = occurrence/sum(occurrence);

figure;
bar(xaxis, prob);

entropy = sum(prob.*log2(1./prob));
disp(entropy);


%% Functions
function [table, codes] = HuffmanDeneme(probs)

    probs = sort(probs);
    
    huffman = transpose(probs);
    
    num_columns = 1;
    
    while num_columns < size(probs, 2) - 1 % Huffman tree loop
    
        huffman(1, num_columns + 1) = huffman(1, num_columns) + huffman(2, num_columns); % New Huffman Probability
        unchanged_huffmans = huffman(3:size(probs,2), num_columns);
        
%         if size(unchanged_huffmans, 1) <= size(probs, 2) - num_columns
%             for i = 1:size(probs,2)-num_columns
%                 unchanged_huffmans = vertcat(unchanged_huffmans, 0)
%             end
%         end
        unchanged_huffmans = vertcat(unchanged_huffmans, 0);
        
        huffman(2:size(probs, 2), num_columns + 1) = unchanged_huffmans;
        
        huffman = huffman + (huffman == 0).*ones(size(huffman));
        
        huffman = sort(huffman);
        
        num_columns = num_columns + 1;
        
    end
    
    code_array = strings(size(probs, 2), size(probs, 2) - 1);
    
    code_array(1, num_columns) = "0";
    code_array(2, num_columns) = "1";
    
    while num_columns > 1 % Huffman coding loop        
          
        flag = 1;
        index = 1;
        
        while (flag == 1)
            
            if huffman(1, num_columns - 1) + huffman(2, num_columns - 1) == huffman(index, num_columns)
                code_array(1, num_columns - 1) = strcat(code_array(index, num_columns), "0");
                code_array(2, num_columns - 1) = strcat(code_array(index, num_columns), "1");
                flag = 0;
            else
                index = index + 1;
            end                   
        
        end
        
        temp_codes = code_array(:, num_columns);
        temp_codes = temp_codes(1:end ~= index);
        placement = 3;
        
        for i = 1:size(temp_codes,1)
            if temp_codes(i) ~= ""
                code_array(placement, num_columns - 1) = temp_codes(i);
                placement = placement + 1;
            end
        end
    
        num_columns = num_columns - 1;
        
    end
    
    table = huffman;
    codes = code_array;

end

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

function [quant_final, recon_final] = LloydMaxFunc(Y, a)

    for i = 1:size(a, 2) - 1
        recon(i) = (a(i) + a(i + 1))/2;
    end

    mask = ones(1,1000000);
    Y_tilde = recon(1)*(mask.*(Y <= a(2))) + recon(2)*(mask.*(Y > a(2) & Y <= a(3))) + recon(3)*(mask.*(Y > a(3) & Y <= a(4))) + recon(4)*(mask.*(Y > a(4) & Y <= a(5))) + recon(5)*(mask.*(Y > a(5) & Y <= a(6))) + recon(6)*(mask.*(Y > a(6) & Y <= a(7))) + recon(7)*(mask.*(Y > a(7) & Y <= a(8))) + recon(8)*(mask.*(Y > a(8)));
    e = Y - Y_tilde;
    
    D_prev = sum(e.^2)/1000000;
    
    delta = 1;
    
    while delta > 0.05
        recon(1) = sum(Y.*(mask.*(Y <= a(2))))/sum(mask.*(Y <= a(2)));
        recon(2) = sum(Y.*(mask.*(Y > a(2) & Y <= a(3))))/sum(mask.*(Y > a(2) & Y <= a(3)));
        recon(3) = sum(Y.*(mask.*(Y > a(3) & Y <= a(4))))/sum(mask.*(Y > a(3) & Y <= a(4)));
        recon(4) = sum(Y.*(mask.*(Y > a(4) & Y <= a(5))))/sum(mask.*(Y > a(4) & Y <= a(5)));
        recon(5) = sum(Y.*(mask.*(Y > a(5) & Y <= a(6))))/sum(mask.*(Y > a(5) & Y <= a(6)));
        recon(6) = sum(Y.*(mask.*(Y > a(6) & Y <= a(7))))/sum(mask.*(Y > a(6) & Y <= a(7)));
        recon(7) = sum(Y.*(mask.*(Y > a(7) & Y <= a(8))))/sum(mask.*(Y > a(7) & Y <= a(8)));
        recon(8) = sum(Y.*(mask.*(Y > a(8))))/sum(mask.*(Y > a(8)));
        
        for i = 2:size(a, 2) - 1          
            a(i) = (recon(i) + recon(i-1))/2;            
        end
        
        Y_tilde = recon(1)*(mask.*(Y <= a(2))) + recon(2)*(mask.*(Y > a(2) & Y <= a(3))) + recon(3)*(mask.*(Y > a(3) & Y <= a(4))) + recon(4)*(mask.*(Y > a(4) & Y <= a(5))) + recon(5)*(mask.*(Y > a(5) & Y <= a(6))) + recon(6)*(mask.*(Y > a(6) & Y <= a(7))) + recon(7)*(mask.*(Y > a(7) & Y <= a(8))) + recon(8)*(mask.*(Y > a(8)));
        
        e = Y - Y_tilde;
        
        D_new = sum(e.^2)/1000000;
        
        delta = abs(D_new - D_prev);
        
    end
    
    quant_final = a;
    recon_final = recon;
    
end