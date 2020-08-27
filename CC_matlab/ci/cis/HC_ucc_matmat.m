function [SIGMA] = HC_ucc_matmat(C,sys)

    SIGMA = zeros(size(C));
    
    for i = 1:size(C,2)
        sigma = HC_ucc_matvec(C(:,i),sys);
        SIGMA(:,i) = sigma;
    end

end
