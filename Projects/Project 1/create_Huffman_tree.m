function tree = create_Huffman_tree(probs)
        
    nodes = struct.empty(0, size(probs,2));

    for i = 1:size(probs,2)
        nodes(i).Prob = probs(i);
        nodes(i).Code = "";
        nodes(i).LeftSub = "";
        nodes(i).RightSub = "";
    end

    while size(nodes,2) > 1
        nodes = create_Huffman_node(nodes);
    end
    
    tree = nodes;
end