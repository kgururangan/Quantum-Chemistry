function [SIGMA] = HR_matmat(R,HBar)

    [Nocc, Nunocc] = size(HBar{1}{1,2});
    Nov = Nocc*Nunocc;

    SIGMA = zeros(size(R));
    
    for i = 1:size(R,2)
        r1 = reshape(R(1:Nov,i),Nunocc,Nocc);
        r2 = reshape(R(Nov+1:end,i),Nunocc,Nunocc,Nocc,Nocc);
        [sigma_r1, sigma_r2] = build_HBar_R(r1,r2,HBar);
        SIGMA(1:Nov,i) = sigma_r1(:);
        SIGMA(Nov+1:end,i) = sigma_r2(:);
    end

end

