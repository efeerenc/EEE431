%Part 2
%initializing the sentence
sentence = 'LOREM IPSUM DOLOR SIT AMET CONSECTETUER ADIPISCING ELIT';

%%
strlength(sentence); % length = 55
x = repmat(sentence, 1, 182); % 10000/55 = 181.8182

dictionary = string.empty; % initialize the character dictionary

for i = 1:strlength(sentence) % this for loop fills the dictionary
    if ~ismember(sentence(i), dictionary)
        dictionary(end + 1) = sentence(i);
    end
end

countTable = double.empty; % initialize the numbers of each character

for i = 1:size(dictionary, 2)
    countTable(end + 1) = count(x, dictionary(i));
end

    