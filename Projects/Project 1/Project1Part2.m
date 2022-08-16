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