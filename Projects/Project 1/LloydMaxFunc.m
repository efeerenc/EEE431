function [quant_final, recon_final] = LloydMaxFunc(Y, a)

    for i = 1:size(a, 2) - 1
        recon(i) = (a(i) + a(i + 1))/2;
    end

    mask = ones(1,1000000);
    Y_tilde = recon(1)*(mask.*(Y <= a(2))) + recon(2)*(mask.*(Y > a(2) & Y <= a(3))) + recon(3)*(mask.*(Y > a(3) & Y <= a(4))) + recon(4)*(mask.*(Y > a(4) & Y <= a(5))) + recon(5)*(mask.*(Y > a(5) & Y <= a(6))) + recon(6)*(mask.*(Y > a(6) & Y <= a(7))) + recon(7)*(mask.*(Y > a(7) & Y <= a(8))) + recon(8)*(mask.*(Y > a(8)));
    e = Y - Y_tilde;
    
    D_prev = sum(e.^2)/1000000;
    
    delta = 1;
    
    while delta > 0.05
        recon(1) = sum(Y.*(mask.*(Y <= a(2))))/sum(mask.*(Y <= a(2)));
        recon(2) = sum(Y.*(mask.*(Y > a(2) & Y <= a(3))))/sum(mask.*(Y > a(2) & Y <= a(3)));
        recon(3) = sum(Y.*(mask.*(Y > a(3) & Y <= a(4))))/sum(mask.*(Y > a(3) & Y <= a(4)));
        recon(4) = sum(Y.*(mask.*(Y > a(4) & Y <= a(5))))/sum(mask.*(Y > a(4) & Y <= a(5)));
        recon(5) = sum(Y.*(mask.*(Y > a(5) & Y <= a(6))))/sum(mask.*(Y > a(5) & Y <= a(6)));
        recon(6) = sum(Y.*(mask.*(Y > a(6) & Y <= a(7))))/sum(mask.*(Y > a(6) & Y <= a(7)));
        recon(7) = sum(Y.*(mask.*(Y > a(7) & Y <= a(8))))/sum(mask.*(Y > a(7) & Y <= a(8)));
        recon(8) = sum(Y.*(mask.*(Y > a(8))))/sum(mask.*(Y > a(8)));
        
        for i = 2:size(a, 2) - 1          
            a(i) = (recon(i) + recon(i-1))/2;            
        end
        
        Y_tilde = recon(1)*(mask.*(Y <= a(2))) + recon(2)*(mask.*(Y > a(2) & Y <= a(3))) + recon(3)*(mask.*(Y > a(3) & Y <= a(4))) + recon(4)*(mask.*(Y > a(4) & Y <= a(5))) + recon(5)*(mask.*(Y > a(5) & Y <= a(6))) + recon(6)*(mask.*(Y > a(6) & Y <= a(7))) + recon(7)*(mask.*(Y > a(7) & Y <= a(8))) + recon(8)*(mask.*(Y > a(8)));
        
        e = Y - Y_tilde;
        
        D_new = sum(e.^2)/1000000;
        
        delta = abs(D_new - D_prev);
        
    end
    
    quant_final = a;
    recon_final = recon;
    
end