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