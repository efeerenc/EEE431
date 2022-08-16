function codewords = denemedeneme(probs)

init_probs = probs;

codewords = strings([1,size(probs,2)]);

while size(probs,2) > 1
   probs = sort(probs, "ascend");
   
   combined_prob = sum(probs(1:2));
   
   
end

end