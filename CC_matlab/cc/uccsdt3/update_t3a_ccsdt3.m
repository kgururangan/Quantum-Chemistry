function [t3a] = update_t3a_ccsdt3(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, HBar_t, VT3_t, sys, shift)

    H1A = HBar_t.H1A;
    H2A = HBar_t.H2A;
    H2B = HBar_t.H2B;
    VTA = VT3_t.A;

    I2A_vvov = H2A.vvov-einsum_kg(H1A.ov,t2a,'me,abim->abie');
    
    % MM23A
    M23_D1 =   -einsum_kg(H2A.vooo + VTA.vooo,t2a,'amij,bcmk->abcijk'); 
    M23_D1 = M23_D1 - permute(M23_D1,[1,2,3,6,5,4]) - permute(M23_D1,[1,2,3,4,6,5]) - permute(M23_D1,[3,2,1,4,5,6])...
                    - permute(M23_D1,[2,1,3,4,5,6]) + permute(M23_D1,[2,1,3,6,5,4]) + permute(M23_D1,[3,2,1,6,5,4]) ...
                    + permute(M23_D1,[2,1,3,4,6,5]) + permute(M23_D1,[3,2,1,4,6,5]);

         
    M23_D2 =   +einsum_kg(I2A_vvov + VTA.vvov,t2a,'abie,ecjk->abcijk');
    M23_D2 = M23_D2 - permute(M23_D2,[1,2,3,5,4,6]) - permute(M23_D2,[1,2,3,6,5,4]) ...
                    - permute(M23_D2,[3,2,1,4,5,6]) - permute(M23_D2,[1,3,2,4,5,6]) + permute(M23_D2,[3,2,1,5,4,6]) ...
                    + permute(M23_D2,[1,3,2,5,4,6]) + permute(M23_D2,[3,2,1,6,5,4]) + permute(M23_D2,[1,3,2,6,5,4]);
        
    MM23A = M23_D1 + M23_D2;

    % (HBar*T3)_C
    
    D1 = -einsum_kg(H1A.oo-diag(diag(sys.fa_oo)),t3a,'mk,abcijm->abcijk');
    
    D2 = einsum_kg(H1A.vv-diag(diag(sys.fa_vv)),t3a,'ce,abeijk->abcijk');
    
    D3 = 0.5*einsum_kg(H2A.oooo,t3a,'mnij,abcmnk->abcijk');
    D13 = D1 + D3;
    D13 = D13 - permute(D13,[1,2,3,6,5,4]) - permute(D13,[1,2,3,4,6,5]);
    
    D4 = 0.5*einsum_kg(H2A.vvvv,t3a,'abef,efcijk->abcijk');
    D24 = D2 + D4;
    D24 = D24 - permute(D24,[3,2,1,4,5,6]) - permute(D24,[1,3,2,4,5,6]);
    
    D5 = einsum_kg(H2A.voov,t3a,'cmke,abeijm->abcijk');
    D6 = einsum_kg(H2B.voov,t3b,'cmke,abeijm->abcijk');
    D56 = D5 + D6;

    D56 = D56 - permute(D56,[1,2,3,6,5,4]) - permute(D56,[1,2,3,4,6,5])...
              - permute(D56,[3,2,1,4,5,6]) - permute(D56,[1,3,2,4,5,6]) ...
              + permute(D56,[3,2,1,6,5,4]) + permute(D56,[3,2,1,4,6,5]) ...
              + permute(D56,[1,3,2,6,5,4]) + permute(D56,[1,3,2,4,6,5]);    
     
     X3A_abcijk = MM23A + D13 + D24 + D56;
     
     t3a = zeros(size(t3a));
     for a = 1:sys.Nvir_alpha
        for b = a+1:sys.Nvir_alpha
            for c = b+1:sys.Nvir_alpha
                for i = 1:sys.Nocc_alpha
                    for j = i+1:sys.Nocc_alpha
                        for k = j+1:sys.Nocc_alpha
                            
                            % (1)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(a,b,c,i,j,k) = X3A_abcijk(a,b,c,i,j,k)/...
                                (sys.fa_oo(i,i)+sys.fa_oo(j,j)+sys.fa_oo(k,k)-sys.fa_vv(a,a)-sys.fa_vv(b,b)-sys.fa_vv(c,c)-shift);    
                            t3a(a,b,c,k,i,j) = t3a(a,b,c,i,j,k);
                            t3a(a,b,c,j,k,i) = t3a(a,b,c,i,j,k);
                            t3a(a,b,c,i,k,j) = -t3a(a,b,c,i,j,k);
                            t3a(a,b,c,j,i,k) = -t3a(a,b,c,i,j,k);
                            t3a(a,b,c,k,j,i) = -t3a(a,b,c,i,j,k);
                            
                            % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(b,a,c,i,j,k) = -t3a(a,b,c,i,j,k);
                            t3a(b,a,c,k,i,j) = -t3a(a,b,c,i,j,k);
                            t3a(b,a,c,j,k,i) = -t3a(a,b,c,i,j,k);
                            t3a(b,a,c,i,k,j) = t3a(a,b,c,i,j,k);
                            t3a(b,a,c,j,i,k) = t3a(a,b,c,i,j,k);
                            t3a(b,a,c,k,j,i) = t3a(a,b,c,i,j,k);
                            
                            % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(a,c,b,i,j,k) = -t3a(a,b,c,i,j,k);
                            t3a(a,c,b,k,i,j) = -t3a(a,b,c,i,j,k);
                            t3a(a,c,b,j,k,i) = -t3a(a,b,c,i,j,k);
                            t3a(a,c,b,i,k,j) = t3a(a,b,c,i,j,k);
                            t3a(a,c,b,j,i,k) = t3a(a,b,c,i,j,k);
                            t3a(a,c,b,k,j,i) = t3a(a,b,c,i,j,k);
                            
                            % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(c,b,a,i,j,k) = -t3a(a,b,c,i,j,k);
                            t3a(c,b,a,k,i,j) = -t3a(a,b,c,i,j,k);
                            t3a(c,b,a,j,k,i) = -t3a(a,b,c,i,j,k);
                            t3a(c,b,a,i,k,j) = t3a(a,b,c,i,j,k);
                            t3a(c,b,a,j,i,k) = t3a(a,b,c,i,j,k);
                            t3a(c,b,a,k,j,i) = t3a(a,b,c,i,j,k);
                            
                            % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(b,c,a,i,j,k) = t3a(a,b,c,i,j,k);
                            t3a(b,c,a,k,i,j) = t3a(a,b,c,i,j,k);
                            t3a(b,c,a,j,k,i) = t3a(a,b,c,i,j,k);
                            t3a(b,c,a,i,k,j) = -t3a(a,b,c,i,j,k);
                            t3a(b,c,a,j,i,k) = -t3a(a,b,c,i,j,k);
                            t3a(b,c,a,k,j,i) = -t3a(a,b,c,i,j,k);
                            
                            % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(c,a,b,i,j,k) = t3a(a,b,c,i,j,k);
                            t3a(c,a,b,k,i,j) = t3a(a,b,c,i,j,k);
                            t3a(c,a,b,j,k,i) = t3a(a,b,c,i,j,k);
                            t3a(c,a,b,i,k,j) = -t3a(a,b,c,i,j,k);
                            t3a(c,a,b,j,i,k) = -t3a(a,b,c,i,j,k);
                            t3a(c,a,b,k,j,i) = -t3a(a,b,c,i,j,k);
                        end
                    end
                end
            end
        end
    end
       

end

