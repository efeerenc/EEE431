function updated = create_Huffman_node(nodes)

    if nodes(1).Prob > nodes(2).Prob
        low = nodes(1).Prob;
        low_index = 1;
        lower = nodes(2).Prob;
        lower_index = 2;
    else
        lower = nodes(1).Prob;
        lower_index = 1;
        low = nodes(2).Prob;
        low_index = 2;
    end
    
    for i = 3:size(nodes,2) %Find the lowest two probabilities
        if (nodes(i).Prob) <= low && (nodes(i).Prob <= lower)
            low = lower;
            low_index = lower_index;
            lower = nodes(i).Prob;
            lower_index = i;
        elseif (nodes(i).Prob <= low) && (nodes(i).Prob > lower)
            low = nodes(i).Prob;
            low_index = i;
        end
    end
    
    new_prob = low + lower;
    
%     nodes(low_index).Code = '1';
%     nodes(lower_index).Code= '0';

    nodes(end+1).Prob = new_prob;
    nodes(end).Code = "";
    nodes(end).LeftSub = nodes(lower_index);
    nodes(end).RightSub = nodes(low_index);
    
    if low_index < lower_index
        nodes(lower_index) = [];
        nodes(low_index) = [];
    else
        nodes(low_index) = [];
        nodes(lower_index) = [];
    end

    updated = nodes;
    
end