function [table, codess] = HuffmanMatrix(probs)
    probs = sort(probs);
    
    huffman = transpose(probs); % Sorted probabilities are listed in columns in ascending order
    
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
    
    codes = cell(size(huffman));
    
    codes{1, num_columns} = '0';
    codes{2, num_columns} = '1';
    
    while num_columns > 1 % Huffman coding loop
    
        if huffman(1, num_columns) == huffman(1, num_columns - 1) + huffman(2, num_columns - 1)
        
            codes{1, num_columns - 1} = strcat(codes{1, num_columns}, '0');
            codes{2, num_columns - 1} = strcat(codes{1, num_columns}, '1');
            
            temp_codes = codes{:, num_columns};
            
            temp_codes{1} = [];           
            
        else 
            
            codes{1, num_columns - 1} = strcat(codes{2, num_columns}, '0');
            codes{2, num_columns - 1} = strcat(codes{2, num_columns}, '1');
        
            temp_codes = codes(:, num_columns);
            
            temp_codes{2} = [];           
            
        end
        
        index = 1;
        nonzero_index = 1;
        nonzero_temp = cell(1, 1);
        
        while index <= size(probs, 2)
        
            if ~isempty(temp_codes{index})
                nonzero_temp{nonzero_index} = temp_codes{index}; 
                nonzero_index = nonzero_index + 1;
            end
            
            index = index + 1;
            
        end
        
        index = 1;
        nonzero_index = 1;
        loop_count = size(nonzero_temp, 2);

        while loop_count > 0
            if isempty(codes{index, num_columns - 1})
                loop_count = loop_count - 1;
                codes{index, num_columns - 1} = nonzero_temp{nonzero_index};
                nonzero_index = nonzero_index + 1;
            end
            
            index = index + 1;
            
            if index > size(probs, 2)
                loop_count = 0;
            end
            
        end
        
        num_columns = num_columns - 1;
        
    end
        
    table = huffman;
    codess = codes;
        
end