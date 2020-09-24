function [MM23] = build_MM23(t1,t2,HBar)

    fprintf('MMCC(2,3) construction... ')
    
    tic

    h_ov = HBar{1}{1,2};
    h_vooo = HBar{2}{2,1,1,1};
    h_vvov = HBar{2}{2,2,1,2};
    I_vvov = h_vvov+einsum_kg(h_ov,t2,'me,abim->abie');
    
    % = -A(k/ij)A(a/bc) h(amij)*t(bcmk)
    h_vooo_t2 =   -einsum_kg(h_vooo,t2,'amij,bcmk->abcijk');
    h_vooo_t2 = h_vooo_t2 - permute(h_vooo_t2,[2,1,3,4,5,6]) - permute(h_vooo_t2,[3,2,1,4,5,6]);
    h_vooo_t2 = h_vooo_t2 - permute(h_vooo_t2,[1,2,3,6,5,4]) - permute(h_vooo_t2,[1,2,3,4,6,5]);
    
    % A(i/jk)A(c/ab)(h(abie)-h(me)*t(abim))*t(ecjk)
    I_vvov_t2 =   +einsum_kg(I_vvov,t2,'abie,ecjk->abcijk');
    I_vvov_t2 = I_vvov_t2 - permute(I_vvov_t2,[3,2,1,4,5,6]) - permute(I_vvov_t2,[1,3,2,4,5,6]);
    I_vvov_t2 = I_vvov_t2 - permute(I_vvov_t2,[1,2,3,5,4,6]) - permute(I_vvov_t2,[1,2,3,6,5,4]);
              
    MM23 = h_vooo_t2 + I_vvov_t2;

    fprintf(' finished in %4.2f s\n',toc)
end

