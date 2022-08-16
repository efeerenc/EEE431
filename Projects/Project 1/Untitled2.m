function codewords = HuffmanMatrix(probs)
    probs = sort(probs);
    
    symbols = cell(1, size(probs, 2));
    huffman = transpose(probs); % Sorted probabilities are listed in columns in ascending order
    
    num_columns = 1
    
    while num_columns < size(probs, 2)
    
        huffman(1, num_columns + 1) = huffman(1, num_columns) + huffman(2, num_columns); % New Huffman Probability
        unchanged_huffmans = huffman(3:size(probs,2), num_columns);
        
        if size(unchanged_huffmans, 1) < size(probs, 2) - num_columns
            for i = 1:size(probs,2)-num_columns
                unchanged_huffmans = vertcat(unchanged_huffmans, 0)
            end
        end
        
        huffman(2:size(probs, 2), num_columns + 1) = unchanged_huffmans;
        
        huffman = sort(huffman);
        
    end
    
    codewords = huffman;
        
end