function [SIGMA] = HR_ucc_matmat(R,HBar_t,cc_t,sys)

    SIGMA = zeros(size(R));
    
    for i = 1:size(R,2)
        sigma = HR_ucc_matvec(R(:,i),HBar_t,cc_t,sys);
        SIGMA(:,i) = sigma;
    end

end
