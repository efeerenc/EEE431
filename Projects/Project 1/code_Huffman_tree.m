function result = code_Huffman_tree(tree)
    
    if isa(tree.LeftSub, 'struct') && isa(tree.RightSub, 'struct')
        tree.LeftSub.Code = strcat(tree.Code, '0');
        tree.RightSub.Code = strcat(tree.Code, '1');      
    end

    flag = 1;
    
    while flag == 1
    
        
    
    end
end