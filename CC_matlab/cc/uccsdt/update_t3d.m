function [t3d] = update_t3d(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, HBar_t, VT3_t, sys, shift)

    H1B = HBar_t.H1B;
    H2B = HBar_t.H2B;
    H2C = HBar_t.H2C;
    VTC = VT3_t.C;

    I2C_vvov = H2C.vvov-einsum_kg(H1B.ov,t2c,'me,abim->abie');
    
    % MM23A
    M23_D1 =   -einsum_kg(H2C.vooo + VTC.vooo,t2c,'amij,bcmk->abcijk'); 
    M23_D1 = M23_D1 - permute(M23_D1,[1,2,3,6,5,4]) - permute(M23_D1,[1,2,3,4,6,5]) - permute(M23_D1,[3,2,1,4,5,6])...
                    - permute(M23_D1,[2,1,3,4,5,6]) + permute(M23_D1,[2,1,3,6,5,4]) + permute(M23_D1,[3,2,1,6,5,4]) ...
                    + permute(M23_D1,[2,1,3,4,6,5]) + permute(M23_D1,[3,2,1,4,6,5]);

         
    M23_D2 =   +einsum_kg(I2C_vvov + VTC.vvov,t2c,'abie,ecjk->abcijk');
    M23_D2 = M23_D2 - permute(M23_D2,[1,2,3,5,4,6]) - permute(M23_D2,[1,2,3,6,5,4]) ...
                    - permute(M23_D2,[3,2,1,4,5,6]) - permute(M23_D2,[1,3,2,4,5,6]) + permute(M23_D2,[3,2,1,5,4,6]) ...
                    + permute(M23_D2,[1,3,2,5,4,6]) + permute(M23_D2,[3,2,1,6,5,4]) + permute(M23_D2,[1,3,2,6,5,4]);
        
    MM23D = M23_D1 + M23_D2;

    % (HBar*T3)_C

    D1 = -einsum_kg(H1B.oo-diag(diag(sys.fb_oo)),t3d,'mk,abcijm->abcijk');

    D2 = einsum_kg(H1B.vv-diag(diag(sys.fb_vv)),t3d,'ce,abeijk->abcijk');

    D3 = 0.5*einsum_kg(H2C.oooo,t3d,'mnij,abcmnk->abcijk');

    D4 = 0.5*einsum_kg(H2C.vvvv,t3d,'abef,efcijk->abcijk');

    D5 = einsum_kg(H2B.ovvo,t3c,'maei,ebcmjk->abcijk');

    D6 = einsum_kg(H2C.voov,t3d,'amie,ebcmjk->abcijk');

    D13 = D1 + D3;
    D13 = D13 - permute(D13,[1,2,3,6,4,5]) - permute(D13,[1,2,3,4,6,5]);

    D24 = D2 + D4;
    D24 = D24 - permute(D24,[3,2,1,4,5,6]) - permute(D24,[1,3,2,4,5,6]);

    D56 = D5 + D6;
    D56 = D56 - permute(D56,[1,2,3,6,4,5]) - permute(D56,[1,2,3,5,4,6])...
              -permute(D56,[2,1,3,4,5,6]) - permute(D56,[3,2,1,4,5,6]) ...
              +permute(D56,[2,1,3,6,4,5]) + permute(D56,[2,1,3,5,4,6]) ...
              +permute(D56,[3,2,1,6,4,5]) + permute(D56,[3,2,1,5,4,6]); 
     
     X3D_abcijk = MM23D + D13 + D24 + D56;
     
     t3d = zeros(size(t3d));
     for a = 1:sys.Nvir_beta
        for b = a+1:sys.Nvir_beta
            for c = b+1:sys.Nvir_beta
                for i = 1:sys.Nocc_beta
                    for j = i+1:sys.Nocc_beta
                        for k = j+1:sys.Nocc_beta
                            
                            % (1)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(a,b,c,i,j,k) = X3D_abcijk(a,b,c,i,j,k)/...
                                (sys.fb_oo(i,i)+sys.fb_oo(j,j)+sys.fb_oo(k,k)-sys.fb_vv(a,a)-sys.fb_vv(b,b)-sys.fb_vv(c,c)-shift);    
                            t3d(a,b,c,k,i,j) = t3d(a,b,c,i,j,k);
                            t3d(a,b,c,j,k,i) = t3d(a,b,c,i,j,k);
                            t3d(a,b,c,i,k,j) = -t3d(a,b,c,i,j,k);
                            t3d(a,b,c,j,i,k) = -t3d(a,b,c,i,j,k);
                            t3d(a,b,c,k,j,i) = -t3d(a,b,c,i,j,k);
                            
                            % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(b,a,c,i,j,k) = -t3d(a,b,c,i,j,k);
                            t3d(b,a,c,k,i,j) = -t3d(a,b,c,i,j,k);
                            t3d(b,a,c,j,k,i) = -t3d(a,b,c,i,j,k);
                            t3d(b,a,c,i,k,j) = t3d(a,b,c,i,j,k);
                            t3d(b,a,c,j,i,k) = t3d(a,b,c,i,j,k);
                            t3d(b,a,c,k,j,i) = t3d(a,b,c,i,j,k);
                            
                            % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(a,c,b,i,j,k) = -t3d(a,b,c,i,j,k);
                            t3d(a,c,b,k,i,j) = -t3d(a,b,c,i,j,k);
                            t3d(a,c,b,j,k,i) = -t3d(a,b,c,i,j,k);
                            t3d(a,c,b,i,k,j) = t3d(a,b,c,i,j,k);
                            t3d(a,c,b,j,i,k) = t3d(a,b,c,i,j,k);
                            t3d(a,c,b,k,j,i) = t3d(a,b,c,i,j,k);
                            
                            % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(c,b,a,i,j,k) = -t3d(a,b,c,i,j,k);
                            t3d(c,b,a,k,i,j) = -t3d(a,b,c,i,j,k);
                            t3d(c,b,a,j,k,i) = -t3d(a,b,c,i,j,k);
                            t3d(c,b,a,i,k,j) = t3d(a,b,c,i,j,k);
                            t3d(c,b,a,j,i,k) = t3d(a,b,c,i,j,k);
                            t3d(c,b,a,k,j,i) = t3d(a,b,c,i,j,k);
                            
                            % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(b,c,a,i,j,k) = t3d(a,b,c,i,j,k);
                            t3d(b,c,a,k,i,j) = t3d(a,b,c,i,j,k);
                            t3d(b,c,a,j,k,i) = t3d(a,b,c,i,j,k);
                            t3d(b,c,a,i,k,j) = -t3d(a,b,c,i,j,k);
                            t3d(b,c,a,j,i,k) = -t3d(a,b,c,i,j,k);
                            t3d(b,c,a,k,j,i) = -t3d(a,b,c,i,j,k);
                            
                            % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(c,a,b,i,j,k) = t3d(a,b,c,i,j,k);
                            t3d(c,a,b,k,i,j) = t3d(a,b,c,i,j,k);
                            t3d(c,a,b,j,k,i) = t3d(a,b,c,i,j,k);
                            t3d(c,a,b,i,k,j) = -t3d(a,b,c,i,j,k);
                            t3d(c,a,b,j,i,k) = -t3d(a,b,c,i,j,k);
                            t3d(c,a,b,k,j,i) = -t3d(a,b,c,i,j,k);
                        end
                    end
                end
            end
        end
    end
       

end