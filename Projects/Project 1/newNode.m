function updated_tree = newNodee(tree) %Creates a new HuffmanNode in tree using two nodes with the lowest probabilities 
    if tree.Nodes(1).Prob > tree.Nodes(2).Prob
        low = tree.Nodes(1).Prob;
        low_index = 1;
        lower = tree.Nodes(2).Prob;
        lower_index = 2;
    else
        lower = tree.Nodes(1).Prob;
        lower_index = 1;
        low = tree.Nodes(2).Prob;
        low_index = 2;
    end            

    for i = 3:size(tree.Nodes,2) %Find the lowest two probabilities
        if (tree.Nodes(i).Prob) <= low && (tree.Nodes(i).Prob <= lower)
            low = lower;
            low_index = lower_index;
            lower = tree.Nodes(i).Prob;
            lower_index = i;
        elseif (tree.Nodes(i).Prob <= low) && (tree.Nodes(i).Prob > lower)
            low = tree.Nodes(i).Prob;
            low_index = i;
        end
    end

    new_prob = low + lower;

    new_node = HuffmanNode;
    new_node.Prob = new_prob;
    new_node.LeftSub = tree.Nodes(lower_index);
    new_node.RightSub = tree.Nodes(low_index);

    if low_index < lower_index
        tree.Nodes(lower_index) = [];
        tree.Nodes(low_index) = [];
    else
        tree.Nodes(low_index) = [];
        tree.Nodes(lower_index) = [];
    end

    tree.Nodes(end + 1) = new_node;

    updated_tree = tree;
end